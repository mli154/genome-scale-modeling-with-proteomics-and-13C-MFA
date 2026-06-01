# Genome scale metabolic modeling with proteomics and 13C-MFA
This repository contains information and MATLAB scripts for genome scale metabolic modeling using empirical data from bulk proteomics and 13C-metabolic flux analysis (13C-MFA). These files are part of the manuscript: Li M, Priem B, Loftus LV, Betenbaugh MJ, Pienta KJ, Amend SR. Polyploid cancer cells surviving cisplatin reallocate central carbon sources to fuel antioxidant metabolism for survival. Mol Metab. 2026 Jun;108:102370. doi: 10.1016/j.molmet.2026.102370. Epub 2026 Apr 18. PMID: 42009164; PMCID: PMC13158425.

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
To use the script for filtering the bulk proteomics dataset to include only metabolism-related proteins, a log2-transformed protein expression .csv file containing ENTREZ IDs for each protein is needed. Example column name and data are shown below:

| Protein.ID | log2_Group1 | log2_Group 2 | ENTREZID |
| -----------|-------------|--------------|----------|
| FBP1       | 20.2087649  | 21.3422368   |2203      |
| SOD2       | 23.8277987  | 26.2988924   |6648      |

In addition, the Recon3D.mat file should be downloaded from the BiGG Model Database [4]. This can be found at: http://bigg.ucsd.edu/models/Recon3D.

Run the code in the metabolic_protein_filtering_proteomics.m file in MATLAB and export the filtered proteomics data as a .csv file.

### Reducing the Recon3D model to only include reactions that correspond to proteomic and 13C-MFA datasets, nutrient exchange, and nutrient transport.
To use the script for reduce the Recon3D model to only include reactions that correspond to proteomic and 13C-MFA datasets, the original Recon3D.mat file, the filtered metabolic proteomic data .csv file, and the 13C-MFA flux result .csv file are needed. 

Run the code in reducing_recon3d_model.m in MATLAB and export the reduced model as a .mat file.

### Genome scale metabolic modeling (GEM)

We provide the script for GEM in the gem_flux_sampling.m file. Open the file in MATLAB. Here is a summary of the steps below:
1. Start up the COBRA Toolbox.
2. Set the solver to the gurobi solver.
3. Load the reduced Recon3D model.
4. Load the filtered proteomics data containing only metabolism-related proteins (log2-transformed) with ENTREZ IDs.
5. Perform normalization of the log2-transformed data for proteomic integration. This step performs a global min-max normalization on the log2-transformed protein expression data. The lowest and highest expression values are identified across all proteins in both experimental groups and all values are scaled from 0 to 1 based on those global minimum and maximum values. These normalized values are used as "expression coefficients" that act as scaling factors for reaction bounds, thus integrating protein expression data into the genome scale metabolic model.
   
```
   % Cross-group normalization of log2 data
all_log2 = [prot_Par; prot_DTP];
global_min = min(all_log2, [], 'omitnan');
global_max = max(all_log2, [], 'omitnan');
norm_Parental = (prot_Par - global_min) / (global_max - global_min);
norm_DTP = (prot_DTP - global_min) / (global_max - global_min);

 ```

6. Scale the reaction bounds based on normalized proteomic expression across groups.
7. Apply 13C-MFA constraints. Load the corresponding 13C-MFA flux results .csv file and set the upper and lower bounds of reactions based on the empirical flux results.
8. Ensure that the objective function is set to 0, since we are interested in objective-independent sampling.
9. Run objective-independent sampling with the Coordinate-Hit-and-Run with Rounding (CHRR) sampler and the Parallel Computing Toolbox. Set the number of samples to 1000, with 100 skips per sample. Save the sample results, model, and roundedPolytope that was generated based on the specific flux constraints and proteomic integration set on the metabolic model.

```
% Run objective-independent flux sampling (chrrSampler)
% using multiple cores
changeCobraSolver('gurobi','all')
numSamples = 1000; numSkip = 100;
parpool('local', 24, 'IdleTimeout', Inf)
[samples_Parental, roundedPolytope_Parental] = chrrSampler(model, numSkip, numSamples, [], [], false); % Use 'false' for useFastFVA

% Save samples, model, and roundedPolytope for specific condition 
save('Parental_samples_model_roundedPolytope.mat', 'samples_Parental', 'model', 'roundedPolytope_Parental')

```

10. Plot results with various types of graphs, such as density plots and violin plots comparing different groups. The code for extracting certain data and plotting are also provided in the gem_flux_sampling.m file.


## References
1. Brunk, E. et al. Recon3D enables a three-dimensional view of gene variation in human metabolism. Nature Biotechnology 36, 272-281 (2018). 
2. Colijn, C. et al. Interpreting Expression Data with Metabolic Flux Models: Predicting Mycobacterium tuberculosis Mycolic Acid Production. PLOS Computational Biology 5, e1000489 (2009).
3. Haraldsdóttir, H. S., Cousins, B., Thiele, I., Fleming, R. M. T. & Vempala, S. CHRR: coordinate hit-and-run with rounding for uniform sampling of constraint-based models. Bioinformatics 33, 1741-1743 (2017).
4. Charles J Norsigian, Neha Pusarla, John Luke McConn, James T Yurkovich, Andreas Dräger, Bernhard O Palsson, Zachary King, BiGG Models 2020: multi-strain genome-scale models and expansion across the phylogenetic tree, Nucleic Acids Research (2020).

