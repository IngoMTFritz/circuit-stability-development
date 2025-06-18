import pymaid
import re 
import numpy as np
from navis.core import TreeNeuron, NeuronList
from navis import config, graph
from typing import Union
from typing_extensions import Literal
# Set up logging
logger = config.get_logger(__name__)

NeuronObject = Union[TreeNeuron, NeuronList]
import scipy
import pandas as pd
from tqdm import tqdm

# Determine significance level based on p-value
def get_significance_star(p_value):
    if p_value <= 0.001:
        return '***'
    elif p_value <= 0.01:
        return '**'
    elif p_value <= 0.05:
        return '*'
    else:
        return 'ns'  # not significant
    
def try_with_retries(func, max_retries=5, *args, **kwargs):
    """
    Try to execute a function with retries. This is useful for functions that use remote servers.

    Parameters:
    - func: the function to call
    - max_retries: how many times to retry
    - *args, **kwargs: arguments to pass to the function

    Returns:
    - the result of func(*args, **kwargs) if successful

    Raises:
    - the last exception if all retries fail
    """
    for attempt in range(max_retries):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            print(f"Attempt {attempt+1}/{max_retries} failed with error: {e}")
            if attempt == max_retries - 1:
                print("All retries failed.")
                raise


def filter_neurons(neurons, name_pattern, print_flag=False):
    # Name will be match pattern "Uniglomerular {tract} DA1 {lineage}"
    prog = re.compile(name_pattern)
    # Match all neuron names in the paper against that pattern
    is_name = list(map(lambda x: prog.match(x) != None, neurons.name))
    # Subset list 
    filtered_neurons = neurons[is_name]
    if print_flag:
        print('Filtering neurons with pattern: ', name_pattern)
        print('Neurons filtered from ', len(neurons), ' to ', len(filtered_neurons), '\n')
    return filtered_neurons


def get_cable_overlap_and_ids(pre: NeuronObject,
                  post: NeuronObject,
                  dist: Union[float, str] = 2,
                  method: Union[Literal['min'], Literal['max'], Literal['mean'],
                                Literal['forward'], Literal['reverse']] = 'min'
                  ) -> pd.DataFrame:
    """Calculate the amount of cable of neuron A within distance of neuron B.
    (Not correct, actually it is only the amount of the clostest nodes/cable of neuron A within distance of neuron B.)

    Parameters
    ----------
    pre,post :       TreeNeuron | NeuronList
                Neuron(s) for which to compute cable within distance. It is
                highly recommended to resample neurons to guarantee an even
                sampling rate.
    dist :      int | float, optional
                Maximum distance. If the neurons have their `.units` set, you
                can also provides this as a string such as "2 microns".
    method :    'min' | 'max' | 'mean' | 'forward' | 'reverse'
                Method by which to calculate the overlapping cable between
                two cables::

                  Assuming that neurons A and B have 300 and 150 um of cable
                  within given distances, respectively:

                    1. 'min' returns 150
                    2. 'max' returns 300
                    3. 'mean' returns 225
                    4. 'forward' returns 300 (i.e. A->B)
                    5. 'reverse' returns 150 (i.e. B->A)

    Returns
    -------
    pandas.DataFrame
            Matrix in which neurons A are rows, neurons B are columns. Cable
            within distance is given in the neuron's native units::

                          neuronD  neuronE   neuronF  ...
                neuronA         5        1         0
                neuronB        10       20         5
                neuronC         4        3        15
                ...

    See Also
    --------
    :func:`navis.resample_skeleton`
                Use to resample neurons before calculating overlap.

    Examples
    --------
    >>> import navis
    >>> nl = navis.example_neurons(4)
    >>> # Cable overlap is given in the neurons' units
    >>> # Converting the example neurons from 8x8x8 voxel space into microns
    >>> # make the results easier to interpret
    >>> nl = nl.convert_units('um')
    >>> # Resample to half a micron
    >>> nl_res = nl.resample('.5 micron', inplace=False)
    >>> # Get overlapping cable within 2 microns
    >>> ol = navis.cable_overlap(nl_res[:2], nl_res[2:], dist='2 microns')

    """

    if not isinstance(pre, (TreeNeuron, NeuronList)) \
       or not isinstance(post, (TreeNeuron, NeuronList)):
        raise TypeError(f'Expected `TreeNeurons`, got "{type(pre)}" and "{type(post)}"')

    if not isinstance(pre, NeuronList):
        pre = NeuronList(pre)

    if not isinstance(post, NeuronList):
        post = NeuronList(post)

    # Make sure neurons have the same units
    # Do not use np.unique here because unit_str can be `None`
    units = set(np.append(pre._unit_str, post._unit_str))
    units = np.array(list(units)).astype(str)
    if len(units) > 1:
        logger.warning('Neurons appear to have different units: '
                       f'{", ".join(units)}. If that is the case, cable '
                       'matrix overlap results will be garbage.')

    allowed_methods = ['min', 'max', 'mean', 'forward', 'reverse']
    if method not in allowed_methods:
        raise ValueError(f'Unknown method "{method}". Allowed methods: '
                         f'"{", ".join(allowed_methods)}"')

    dist = pre[0].map_units(dist, on_error='raise')

    matrix = pd.DataFrame(np.zeros((pre.shape[0], post.shape[0])),
                          index=pre.id, columns=post.id)

    # Initialize an empty list to store data for each neuron
    overlap_loc = pd.DataFrame()

    # Compute required props
    treesPre = []
    lengthsPre = []
    for nPre in pre:
        points, vect, length = graph.neuron2tangents(nPre)
        treesPre.append(scipy.spatial.cKDTree(points))
        lengthsPre.append(length)

    treesPost = []
    lengthsPost = []
    for nPost in post:
        points, vect, length = graph.neuron2tangents(nPost)
        treesPost.append(scipy.spatial.cKDTree(points))
        lengthsPost.append(length)


    with config.tqdm(total=len(pre), desc='Calc. overlap',
                     disable=config.pbar_hide,
                     leave=config.pbar_leave) as pbar:
        for i, nPre in enumerate(pre):
            # Get cKDTree for nA
            tPre = treesPre[i]

            for k, nPost in enumerate(post):
                # Get cKDTree for nB
                tPost = treesPost[k]

                # # Note: ixPre and ixPost can be redundant (closest to datapoints), while valid_tPre and valid_tPost are always unique
                # # Query nB -> nA
                # # Gives all the closest distances from B to A, and all the indices of corresponding A-nodes (can be redundant)
                dist_to_Pre_loc, ixPre_loc = tPre.query(nPost.nodes[['x', 'y', 'z']].values,
                                        k=1,
                                        distance_upper_bound=dist,
                                        workers=1
                                        )
                valid_Post_loc = np.where(dist_to_Pre_loc != float('inf'))[0]  # all post indices that have overlap in proximity 

                # Query nA -> nB
                # Gives all the closest distances from A to B, and all the indices of corresponding B-nodes (can be redundant)
                dist_to_Post_loc, ixPost_loc = tPost.query(nPre.nodes[['x', 'y', 'z']].values,
                                        k=1,
                                        distance_upper_bound=dist,
                                        workers=1
                                        )
                valid_Pre_loc = np.where(dist_to_Post_loc != float('inf'))[0]  # all pre indices that have overlap in proximity


                # Calculate length of each node, by calculating the distance to the parent and child and dividing by two
                if method != 'reverse':
                    nPost_node_lengths = np.array([])
                    for _, node in nPost.nodes.iloc[valid_Post_loc].iterrows():
                        if node.parent_id == -1:
                            nPost_node_lengths = np.append(nPost_node_lengths, 0)
                        else:
                            parent_loc = nPost.nodes[nPost.nodes.node_id == node.parent_id][['x', 'y', 'z']].values
                            node_vect = node[['x', 'y', 'z']].values - parent_loc
                            nPost_node_lengths = np.append(nPost_node_lengths, np.sqrt(np.sum(node_vect ** 2)))

                if method != 'forward':
                    nPre_node_lengths = np.array([])
                    for _, node in nPre.nodes.iloc[valid_Pre_loc].iterrows():
                        if node.parent_id == -1:
                            nPre_node_lengths = np.append(nPre_node_lengths, 0)
                        else:
                            parent_loc = nPre.nodes[nPre.nodes.node_id == node.parent_id][['x', 'y', 'z']].values
                            node_vect = node[['x', 'y', 'z']].values - parent_loc
                            nPre_node_lengths = np.append(nPre_node_lengths, np.sqrt(np.sum(node_vect ** 2)))

                if method == 'mean':
                    overlap = (nPre_node_lengths.sum() + nPost_node_lengths.sum()) / 2
                elif method == 'max':
                    overlap = max(nPre_node_lengths.sum(), nPost_node_lengths.sum())
                elif method == 'min':
                    overlap = min(nPre_node_lengths.sum(), nPost_node_lengths.sum())
                elif method == 'forward':
                    overlap = nPost_node_lengths.sum()
                elif method == 'reverse':
                    overlap = nPre_node_lengths.sum()

                matrix.iloc[i, k] = overlap

                if method == 'forward':
                    # Extract overlap locations and unique indices
                    overlap_nodes_loc = nPost.nodes[['x', 'y', 'z']].values[valid_Post_loc]
                    overlap_nodes_id =  nPost.nodes.node_id[valid_Post_loc]
                    distances = dist_to_Pre_loc[valid_Post_loc]
                    nodes_length = nPost_node_lengths

                    # Create a DataFrame with neuron metadata and potential locations
                    overlap_df = pd.DataFrame({
                        'source_name': [nPre.name] * len(valid_Post_loc),
                        'target_name': [nPost.name] * len(valid_Post_loc),
                        'skeleton_id': [nPre.skeleton_id] * len(valid_Post_loc),
                        'x': overlap_nodes_loc[:, 0],
                        'y': overlap_nodes_loc[:, 1],
                        'z': overlap_nodes_loc[:, 2],
                        'distance': distances,
                        'length': nodes_length,
                        'target_node_index': overlap_nodes_id
                        })
                    # locations
                    overlap_loc = pd.concat([overlap_loc, overlap_df])
                elif method == 'reverse':
                    # Extract overlap locations and unique indices
                    overlap_nodes_loc = nPre.nodes[['x', 'y', 'z']].values[valid_Pre_loc]
                    overlap_nodes_id =  nPre.nodes.node_id[valid_Pre_loc]
                    distances = dist_to_Post_loc[valid_Pre_loc]
                    nodes_length = nPre_node_lengths

                    # Create a DataFrame with neuron metadata and potential locations
                    overlap_df = pd.DataFrame({
                        'source_name': [nPre.name] * len(valid_Pre_loc),
                        'target_name': [nPost.name] * len(valid_Pre_loc),
                        'skeleton_id': [nPre.skeleton_id] * len(valid_Pre_loc),
                        'x': overlap_nodes_loc[:, 0],
                        'y': overlap_nodes_loc[:, 1],
                        'z': overlap_nodes_loc[:, 2],
                        'distance': distances,
                        'length': nodes_length,
                        'target_node_index': overlap_nodes_id
                        })
                    # locations
                    overlap_loc = pd.concat([overlap_loc, overlap_df])

            pbar.update(1)

    return matrix, overlap_loc