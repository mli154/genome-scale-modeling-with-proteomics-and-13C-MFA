# Genome scale metabolic modeling with proteomics and 13C-MFA
This repository contains information and a MATLAB script for genome scale metabolic modeling using empirical data from bulk proteomics and 13C-metabolic flux analysis (13C-MFA).

Here, we implement protein expression data from bulk proteomics and empirical flux data from 13C-metabolic flux analysis to infer fluxes of various pathways in cancer cells at the genome scale. We reduced the Homo sapien Recon3D metabolic network model [1] to only contain reactions reflected in our 13C-MFA and bulk proteomics datasets. Our approach is based off of the E-Flux method [2], where protein expression is used to set upper and lower bounds of reaction fluxes of the reduced Recon3D model, and with the assumption that higher expression leads to higher flux. Flux data from 13C-MFA was used to set constraints for corresponding reactions in the reduced model and Coordinate-Hit-and-Run with Rounding (CHRR) flux sampling [3] was performed to sample flux distributions from the constrined metabolic model.

## Requirements
The following optimization solver and MATLAB toolboxes are required to run this code:
* Constraint-Based Optimization and Reconstruction Analysis (COBRA) Toolbox
* Optimization Toolbox
* Parallel Computing Toolbox
* Statistics and Machine Learning Toolbox
* Gurobi solver (obtain license and install using these instructions: https://support.gurobi.com/hc/en-us/articles/4533938303505-How-do-I-install-Gurobi-for-Matlab

## Usage



















## References
1. Brunk, E. et al. Recon3D enables a three-dimensional view of gene variation in human metabolism. Nature Biotechnology 36, 272-281 (2018). 
2. Colijn, C. et al. Interpreting Expression Data with Metabolic Flux Models: Predicting Mycobacterium tuberculosis Mycolic Acid Production. PLOS Computational Biology 5, e1000489 (2009).
3. Haraldsdóttir, H. S., Cousins, B., Thiele, I., Fleming, R. M. T. & Vempala, S. CHRR: coordinate hit-and-run with rounding for uniform sampling of constraint-based models. Bioinformatics 33, 1741-1743 (2017).

