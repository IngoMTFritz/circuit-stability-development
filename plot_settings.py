import matplotlib as plt

params = {
        #   "font.family": "Arial",
          "legend.fontsize": 14,
          "legend.handlelength": 1.5,
          "legend.edgecolor": 'inherit',
          "legend.columnspacing": 0.8,
          "legend.handletextpad": 0.5,
          "axes.labelsize": 18,
          "axes.spines.right": False,
          "axes.spines.top": False,
          "axes.edgecolor": "#d3d3d3",
          "axes.linewidth": 2,
          "xtick.color": '#d3d3d3',
          "xtick.labelcolor": 'black',
          "xtick.labelsize": 14,
          "xtick.major.size": 6,
          "xtick.major.width": 2,
          "xtick.minor.size": 4,
          "xtick.minor.width": 1,
          "ytick.color": '#d3d3d3',
          "ytick.labelcolor": 'black',
          "ytick.labelsize": 14,
          "ytick.major.size": 6,
          "ytick.major.width": 2,
          "ytick.minor.size": 4,
          "ytick.minor.width": 1,
          'mathtext.default': 'regular',
          'lines.markersize': 6,
          'lines.linewidth': 2,
          'grid.linestyle': ":",
          'grid.color': "#d3d3d3",
          'svg.fonttype': 'none',
}

default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# l1_color = '#E6BD52'
l1_color = '#97c7dd'
# l3_color = '#94416A'
l3_color = '#645aa3'

color_dict = {'L1': l1_color, 'L3': l3_color}

axon_color = '#DF5B5F'
dendrite_color = '#20B49C'

# Neuron colors
neuron_dict = {
    'ddaC': '#cd5241ff',
    'v\'ada': '#ee8329ff',
    'vdaB': '#eede7bff',
    'A02n': 'cornflowerblue',
    'A02m': 'lightblue',
    'A09a': 'seagreen',
    'A09c': 'lightgreen',
    'A09l': 'mediumorchid',
    'A10a': 'pink'
}

# left and right different colors
# neuron_color_dict = {
#     'ddaC_a1l': '#681883ff',
#     'ddaC_a1r': '#ad28d8ff',
#     'v\'ada_a1l': '#2471ddff',
#     'v\'ada_a1r': '#85afebff',
#     'vdaB_a1l': '#00a774ff',
#     'vdaB_a1r': '#23ffbcff',
#     'ddaC_a3l': '#681883ff',
#     'ddaC_a3r': '#ad28d8ff',
#     'v\'ada_a3l': '#2471ddff',
#     'v\'ada_a3r': '#85afebff',
#     'vdaB_a3l': '#00a774ff',
#     'vdaB_a3r': '#23ffbcff',
# }

# left and right same colors
neuron_color_dict = {
    'ddaC_a1l': '#ad28d8ff',
    'ddaC_a1r': '#ad28d8ff',
    'v\'ada_a1l': '#85afebff',
    'v\'ada_a1r': '#85afebff',
    'vdaB_a1l': '#00eaa2ff',
    'vdaB_a1r': '#00eaa2ff',
    'ddaC_a3l': '#ad28d8ff',
    'ddaC_a3r': '#ad28d8ff',
    'v\'ada_a3l': '#85afebff',
    'v\'ada_a3r': '#85afebff',
    'vdaB_a3l': '#00eaa2ff',
    'vdaB_a3r': '#00eaa2ff',
}

## Camera settings
l1_scene_camera=dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=0, y=1.7, z=0))
l3_scene_camera=dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=0, y=-1.7, z=0))

## Figure layout
# Define your default layout parameters
default_layout_L1 = {
    'autosize': False,
    'width': 1500,
    'height': 1000,
    'scene_camera': l1_scene_camera,
    'scene': {
        'xaxis': {'visible': False},
        'yaxis': {'visible': False},
        'zaxis': {'visible': False}
    }
}
default_layout_L3 = {
    'autosize': False,
    'width': 1500,
    'height': 1000,
    'scene_camera': l3_scene_camera,
    'scene': {
        'yaxis_autorange':'reversed', 'zaxis_autorange': 'reversed',
        'xaxis': {'visible': False},
        'yaxis': {'visible': False},
        'zaxis': {'visible': False}
    }
}