# Genome scale metabolic modeling with proteomics and 13C-MFA
This repository contains information and a MATLAB script for genome scale metabolic modeling using empirical data from bulk proteomics and 13C-metabolic flux analysis (13C-MFA).

Here, we implement protein expression data from bulk proteomics and empirical flux data from 13C-metabolic flux analysis to infer fluxes of various pathways in cancer cells at the genome scale. We reduced the Homo sapien Recon3D metabolic network model [1] to only contain reactions reflected in our 13C-MFA and bulk proteomics datasets. Our approach is based off of the E-Flux method [2], where protein expression is used to set upper and lower bounds of reaction fluxes of the reduced Recon3D model, with the assumption that higher expression leads to higher flux. Flux data from 13C-MFA was used to set constraints for corresponding reactions in the reduced model and Coordinate-Hit-and-Run with Rounding (CHRR) flux sampling [3] was performed to sample flux distributions from the constrined metabolic model.

## Requirements
The following optimization solver and MATLAB toolboxes are required to run this code:
* Constraint-Based Optimization and Reconstruction Analysis (COBRA) Toolbox
* Optimization Toolbox
* Parallel Computing Toolbox
* Statistics and Machine Learning Toolbox
* Gurobi solver (obtain license and install using these instructions: https://support.gurobi.com/hc/en-us/articles/4533938303505-How-do-I-install-Gurobi-for-Matlab

## Usage

### Filtering proteomic data to include only metabolism-related proteins based on Recon3D model 
To use script for filtering the bulk proteomics dataset to include only metabolism-related proteins, a log2-transformed protein expression .csv file containing ENTREZ IDs for each protein is needed. Example column name and data are shown below:

| Protein.ID | log2_Group1 | log2_Group 2 |
| -----------|-------------|--------------|
| FBP1       | 20.2087649  | 21.3422368   |
| SOD2       | 23.8277987  | 26.2988924   |

In addition, the Recon3D.mat file should be downloaded from the BiGG Model Database [4]. This can be found at: http://bigg.ucsd.edu/models/Recon3D.

Open the metabolic_protein_filtering_proteomics.m file in MATLAB. 





## References
1. Brunk, E. et al. Recon3D enables a three-dimensional view of gene variation in human metabolism. Nature Biotechnology 36, 272-281 (2018). 
2. Colijn, C. et al. Interpreting Expression Data with Metabolic Flux Models: Predicting Mycobacterium tuberculosis Mycolic Acid Production. PLOS Computational Biology 5, e1000489 (2009).
3. Haraldsdóttir, H. S., Cousins, B., Thiele, I., Fleming, R. M. T. & Vempala, S. CHRR: coordinate hit-and-run with rounding for uniform sampling of constraint-based models. Bioinformatics 33, 1741-1743 (2017).
4. Charles J Norsigian, Neha Pusarla, John Luke McConn, James T Yurkovich, Andreas Dräger, Bernhard O Palsson, Zachary King, BiGG Models 2020: multi-strain genome-scale models and expansion across the phylogenetic tree, Nucleic Acids Research (2020).

