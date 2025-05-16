# Synaptic density and relative connectivity conservation maintain circuit stability across development

# Overview of the Repository

This repository contains the code and data used to generate the figures for the publication. Each figure is associated with a specific Jupyter notebook, which contains the analysis and visualization steps. The repository is organized as follows:

## Repository Structure

- **Jupyter Notebooks**:  
  Each figure has a corresponding notebook:
  - `Fig2_S1.ipynb`: Analysis and plots for Figure 2 and Supplementary Figure 1.
  - `Fig3_S2.ipynb`: Analysis and plots for Figure 3 and Supplementary Figure 2.
  - `Fig4_S3.ipynb`: Analysis and plots for Figure 4 and Supplementary Figure 3.
  - `Fig5_S4_S5.ipynb`: Analysis and plots for Figure 5 and Supplementary Figures 4 and 5.

- **Data Files**:  
  - `l1_PSD_areas.pkl` and `l3_PSD_areas.pkl`: Preprocessed data files used in the analysis.

- **Plots Directory**:  
  - Contains subdirectories for each figure (`Fig2/`, `Fig3/`, etc.) where generated plots are saved.

- **Source Code**:  
  - `src/helper.py`: Contains utility functions for data processing and analysis.
  - `src/plot_settings.py`: Contains matplotlib style settings for consistent figure formatting.

## Required Toolboxes and Libraries

To run the notebooks, the following Python libraries are required:

- **Specialized Libraries**:
  - `pymaid`: For interacting with CATMAID instances.
  - `navis`: For neuron analysis and visualization.

## Instructions for Use

1. Clone the repository:
   ```sh
   git clone https://github.com/IngoMTFritz/circuit-stability-development.git
   cd https://github.com/IngoMTFritz/circuit-stability-development.git