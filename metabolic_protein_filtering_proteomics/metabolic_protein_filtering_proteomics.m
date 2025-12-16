% Startup
addpath(genpath('C:\COBRA'))
savepath
initCobraToolbox(false)

% Load the COBRA metabolic model
model = readCbModel('Recon3D.mat');

% Extract unique metabolic gene IDs from the model
model_gene_ids = cellfun(@(x) strtok(x, '_'), model.genes, 'UniformOutput', false);
model_gene_ids = unique(strtrim(model_gene_ids));  % <--- Remove spaces

fprintf('Number of unique metabolic genes: %d\n', numel(model_gene_ids));

% Read the cleaned proteomics dataset (must include ENTREZID column)
proteomics_tbl = readtable('log2_transformed_proteomics_data_average_cleaned_WITH_ENTREZID.csv');

% Ensure ENTREZID column is formatted as strings with no decimals or spaces
% Convert to cellstr if numeric
if isnumeric(proteomics_tbl.ENTREZID)
    proteomics_tbl.ENTREZID = cellstr(num2str(proteomics_tbl.ENTREZID));
end

% Remove leading/trailing spaces
proteomics_tbl.ENTREZID = strtrim(proteomics_tbl.ENTREZID);

% Filter the proteomics data to include only metabolic genes present in the model
is_in_model = ismember(proteomics_tbl.ENTREZID, model_gene_ids);
filtered_proteomics = proteomics_tbl(is_in_model, :);

% Display and report the filtered results
disp('First few filtered rows:');
disp(filtered_proteomics(1:min(10,height(filtered_proteomics)), :));
fprintf('Filtered to metabolic proteins: %d\n', height(filtered_proteomics));

% Export the filtered proteomics dataset to a new CSV file
writetable(filtered_proteomics, 'metabolic_only_proteomics.csv');