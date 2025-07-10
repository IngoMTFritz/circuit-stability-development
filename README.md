# Synaptic density and relative connectivity conservation maintain circuit stability across development

## Overview of the Repository

This repository contains the code and data used to generate the figures for the publication. Each figure is associated with a specific script or Jupyter notebook, which contains the analysis and visualization steps. The repository is organized as follows:

## Repository Structure

- **Scripts and Jupyter Notebooks**:
  - `Fig2_S1.ipynb`: Analysis and plots for Figure 2 and Supplementary Figure 1.
  - `Fig3.ipynb`: Analysis and plots for Figure 3.
  - `Fig4_S2.ipynb`: Analysis and plots for Figure 4 and Supplementary Figure 2.
  - `Fig5_S3.ipynb`: Analysis and plots for Figure 5 and Supplementary Figure 3.
  - `Fig6_S4.ipynb`: Analysis and plots for Figure 6 and Supplementary Figure 4.
  - `MATLAB model/Figure_7AB_S5_mean_levels.m`: MATLAB script for analysis and plots related to Figure 7AB and Supplementary Figure 5.
  - `MATLAB model/Figure_7CDE_real_rel_weights.m`: MATLAB script for analysis and plots related to Figure 7CDE.

- **Data Files**:
  - `l1_PSD_areas.pkl` and `l3_PSD_areas.pkl`: Preprocessed data files used in the analysis (for Python notebooks).
  - `MATLAB model/data_swc_processed/`: Folder containing SWC files iwht the neuron morphology data used by the MATLAB scripts.
  - `MATLAB model/data/density_l1_l3_dendrite.csv`: Data file for synapse densities used by the MATLAB scripts.
  - `MATLAB model/data/t.mat`: Output data file from the MATLAB script `Figure_6B_S6_mean_levels.m`.

- **Plots Directory**:
  - `plots/`: Contains subdirectories for each figure (`Fig2/`, `Fig3/`, etc.) where generated plots are saved.
  - `MATLAB model/output_plots/`: Plots generate by the MATLAB scripts.

- **Source Code**:
  - `helper.py`: Contains utility functions for data processing and analysis (for Python notebooks).
  - `plot_settings.py`: Contains matplotlib style settings for consistent figure formatting (for Python notebooks).
  - `MATLAB model/compute_r2.m`: MATLAB function used for R-squared calculations.
  - `MATLAB model/compute_overshoot.m`: MATLAB function used for calculating the necessary overshoot voltage to account for synaptic conductance effects at a distal dendritic site.
  - `MATLAB model/plotBinnedAverages.m`: MATLAB function used for plotting binned averages.
  - `MATLAB model/process_swc_synapses.m`: MATLAB script for processing SWC synapse data.
  - `MATLAB model/rel_weights.m`: MATLAB script used to get the relative weights of a neuron.

## Required Toolboxes and Libraries

### Python
To run the Python notebooks, the following Python libraries are required:

- **Specialized Libraries**:
  - [`pymaid`](https://github.com/navis-org/pymaid): For interacting with CATMAID instances. Install using:
    ```bash
    pip3 install python-catmaid
    ```
  - [`navis`](https://github.com/navis-org/navis): For neuron analysis and visualization. Install with all extras using:
    ```bash
    pip3 install "navis[all]"
    ```
- **General Libraries**: Add any other common libraries like `numpy`, `pandas`, `matplotlib`, `scipy` if they are used in the notebooks.

### MATLAB
To run the MATLAB scripts, you will need:

- **TREES Toolbox**: Most MATLAB scripts rely on functions from the [TREES Toolbox](https://www.treestoolbox.org/).

## Instructions for Use

1. Clone the repository:
   ```bash
   git clone https://github.com/IngoMTFritz/circuit-stability-development.git
