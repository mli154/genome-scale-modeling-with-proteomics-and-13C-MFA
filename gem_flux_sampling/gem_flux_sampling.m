% GEM Flux Sampling Script

% Startup
addpath(genpath('C:\COBRA'))
savepath
initCobraToolbox(false)
changeCobraSolver('gurobi','all')

% Load the model
model = readCbModel('Recon3D_reduced.mat');  % Use your SBML (or .mat if you have one)

% Load proteomics data (metabolic proteins CSV)
protT = readtable('metabolic_only_proteomics.csv');
entrezID = string(protT.ENTREZID);
prot_Par = protT.log2_Parental;
prot_DTP = protT.log2_DTP;

% Cross-group normalization of log2 data
all_log2 = [prot_Par; prot_DTP];
global_min = min(all_log2, [], 'omitnan');
global_max = max(all_log2, [], 'omitnan');
norm_Parental = (prot_Par - global_min) / (global_max - global_min);
norm_DTP = (prot_DTP - global_min) / (global_max - global_min);

% Flexible: apply to parental or resistant
group = 'DTP'; % 'DTP' or 'parental'
if strcmp(group, 'parental')
    norm_vec = norm_parental;
else
    norm_vec = norm_DTP;
end

%extract EntrezIDs out of model.genes
modelEntrez = cellfun(@(g) strtok(g, '_'), model.genes, 'UniformOutput', false); % returns cell array of EntrezIDs as strings

% 3. Proteomics scaling (reactxbounds)
maxFlux = 10000;
minScale = 0.05;
for i = 1:numel(model.rxns)
    rule_str = model.rules{i};
    if isempty(rule_str)
        continue
    end
    idx_nums = regexp(rule_str, 'x\((\d+)\)', 'tokens');
    if isempty(idx_nums), continue, end
    idxs = unique(str2double([idx_nums{:}]));
    entrez_ids = string(modelEntrez(idxs));
    [~, locs] = ismember(entrez_ids, entrezID); % indices into protT
    abund = norm_vec(locs(locs>0));             % nonzero locs found
    if ~isempty(abund) && all(~isnan(abund))
        scale = max(mean(abund), minScale);
        if model.lb(i) < 0
            model.lb(i) = -maxFlux * scale;
        else
            model.lb(i) = 0;
        end
        model.ub(i) = maxFlux * scale;
    else
        model.ub(i) = 1;
        if model.lb(i) < 0
            model.lb(i) = -1;
        else
            model.lb(i) = 0;
        end
    end
end

% 4. Apply MFA constraints
mfaT = readtable('Parental_13C-MFA_flux_results_for_GEM.csv');
minStd = 0.05;
for k = 1:height(mfaT)
    rxnid = mfaT.reaction_id{k};
    idx = find(strcmpi(model.rxns, rxnid));
    if ~isempty(idx)
        flux = mfaT.flux(k);
        std = max(mfaT.std(k), minStd);
        lb = flux - std;
        ub = flux + std;
        model.lb(idx) = lb;
        model.ub(idx) = ub;
    end
end

%% for oxPPP KO ONLY
to_KO = {'G6PDH2r', 'G6PDH2c', 'G6PDHer', 'GNDer', 'PGLer', 'GND', 'GNDc', 'PGL', 'PGLc'}; % Include all present in your rxn list
for i = 1:numel(to_KO)
    rxnIdx = find(strcmpi(model.rxns, to_KO{i}));
    if ~isempty(rxnIdx)
        model.lb(rxnIdx) = 0;
        model.ub(rxnIdx) = 0;
    end
end

% Check if KO worked
for i = 1:numel(to_KO)
    rxnIdx = find(strcmpi(model.rxns, to_KO{i}));
    if ~isempty(rxnIdx)
        disp([to_KO{i} ': lb=' num2str(model.lb(rxnIdx)) ', ub=' num2str(model.ub(rxnIdx))])
    end
end

%% Ensure that objective is cleared for objective-independent sampling
model.c(:) = 0;

idx_obj = find(model.c);
if isempty(idx_obj)
    disp('No objective is set (model.c is all zeros). Sampling will be objective-independent.');
else
    disp('Objective(s) set for:');
    disp(model.rxns(idx_obj));
    disp(['Objective coefficients: ', num2str(model.c(idx_obj)')]);
end

% Run objective-independent flux sampling (chrrSampler)
% using multiple cores
changeCobraSolver('gurobi','all')
numSamples = 1000; numSkip = 100;
parpool('local', 24, 'IdleTimeout', Inf)
[samples_Parental, roundedPolytope_Parental] = chrrSampler(model, numSkip, numSamples, [], [], false); % Use 'false' for useFastFVA

% Save samples, model, and roundedPolytope for specific condition 
save('Parental_samples_model_roundedPolytope.mat', 'samples_Parental', 'model', 'roundedPolytope_Parental')

%% To plot the distribution for a particular reaction:
rxnID = 'GTHOr';  % or your reaction of interest
rxnIdx = find(strcmpi(model.rxns, rxnID));

% Get the sampled fluxes for this reaction (a 1 × nSamples vector)
rxnSamples_Parental = samples_Parental(rxnIdx, :);  % Still a row; if needed, use rxnSamples(:) for column

[x_density, y_density] = ksdensity(rxnSamples_Parental(:));   % Estimates the density

figure;
plot(y_density, x_density, 'LineWidth', 2)
xlabel([rxnID ' Flux']);
ylabel('Probability Density');
title(['CHRR sampled flux distribution: ' rxnID]);

%% Plotting density plots of fluxes to compare between groups
% Note: must run flux sampling on another group for comparison
% Modify the script to run flux sampling for another group (ie. DTP, DTP_oxPPP KO)

rxnID = 'GTHOr';  % or your reaction of interest
rxnIdx = find(strcmpi(model.rxns, rxnID));

% Get sampled fluxes for your reaction in each group
rxnSamples_Parental = samples_Parental(rxnIdx, :);
rxnSamples_DTP = samples_DTP(rxnIdx, :);

% Estimate density for both groups (swap output if needed for your MATLAB version)
[x_density_Parental, y_density_Parental] = ksdensity(rxnSamples_Parental(:));
[x_density_DTP, y_density_DTP] = ksdensity(rxnSamples_DTP(:));

figure;
% Plot Parental group
plot(y_density_Parental, x_density_Parental, 'k-', 'LineWidth', 2)
hold on

% Plot DTP group
plot(y_density_DTP, x_density_DTP, 'r-', 'LineWidth', 2)
hold off

xlabel([rxnID ' Flux'])
ylabel('Probability Density')
title([rxnID ' Flux Distribution'])
legend({'Parental', 'DTP'}, 'Location', 'Best')

%% Statistics and plotting violin plots
% Prepare data
fluxes = [samples_Parental(rxnIdx, :)'; samples_Parental(rxnIdx, :)'];
groups = [repmat({'Parental'}, size(samples_Parental,2), 1); 
          repmat({'DTP'}, size(samples_DTP,2), 1)];
% Ensure fluxes are numeric
if ~isnumeric(fluxes)
    error('fluxes must be a numeric array.');
end

% Bootstrapping Mann-Whitney p-values
nBoots = 100; % Number of bootstrap trials
pvals = zeros(nBoots,1);

for b = 1:nBoots
    idxA = randsample(numel(rxnSamples_Parental), 100);
    idxB = randsample(numel(rxnSamples_DTP), 100);
    pvals(b) = ranksum(rxnSamples_Parental(idxA), rxnSamples_DTP(idxB));
end

% Report summary statistics
pval_median = median(pvals);
pval_mean   = mean(pvals);

fprintf('Bootstrapped Mann-Whitney p-value: median = %.3g, mean = %.3g\n', pval_median, pval_mean);


% Violinplot (if using Ben Hart's or similar violinplot.m)
group_labels = categorical(groups, {'Parental', 'DTP'}, 'Ordinal', true);
ax = gca;
ax.ColorOrder = [0 0 0; 1 0 0];     % Black for group 1, red for group 2
ax.ColorOrderIndex = 1;

violinplot(fluxes, group_labels, 'ShowData', false);

ylimVals = ylim;
text(1.5, ylimVals(2), ['p = ' num2str(pval_median, 2)], ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold');

%% Exporting sampled data and reaction IDs from reduced model

% Labeling reactions in the data output (wasn't included in the export)
flux_table_Parental = array2table(samples_Parental); % shape: [6032 x 1000]
flux_table_Parental.Reaction = model.rxns; % add a column for reaction names

% Move reaction names to the first column
flux_table_Parental = [flux_table_Parental(:,end), flux_table_Parental(:,1:end-1)];

% Export as csv
writetable(flux_table_Parental, 'reaction_flux_samples_Parental.csv');

%% Identifying, extracting, and plotting NADPH consuming reactions
% Assuming that flux sampling was performed on DTP and DTP_oxPPP_KO groups
% and that flux_table_DTP and flux_table_DTP_oxPPP_KO are generated (view
% previous section)

% Identifying NADPH consuming reactions
nadph_met_idx = find(strcmp(model.mets, 'nadph_c'));    % index of NADPH
nadph_stoich_row = model.S(nadph_met_idx, :);            % stoichiometry vector for NADPH across all rxns

% NADPH-consuming reactions: negative coefficient (NADPH appears as a reactant)
nadph_consuming_rxn_idxs = find(nadph_stoich_row < -1e-8);   % threshold to avoid rounding errors

nadph_consuming_rxns = model.rxns(nadph_consuming_rxn_idxs);

% Extract NADPH consuming fluxes from labeled dataframe 

% Find rows (reactions) matching NADPH consumers for DTP
rows_to_extract_DTP = ismember(flux_table_DTP.Reaction, nadph_consuming_rxns);

% Extract subtable of NADPH-consuming reactions for DTP
nadph_flux_table_DTP = flux_table_DTP(rows_to_extract_DTP, :);  % All samples, only NADPH-consuming reactions


writetable(nadph_flux_table_DTP, 'nadph_consuming_flux_samples_DTP.csv');

% Find rows (reactions) matching NADPH consumers for DTP oxPPP KO
rows_to_extract_DTP_oxPPP_KO = ismember(flux_table_DTP_oxPPP_KO.Reaction, nadph_consuming_rxns);

% Extract subtable of NADPH-consuming reactions for oxPPP KO
nadph_flux_table_DTP_oxPPP_KO = flux_table_PC3_DTP_oxPPP_KO(rows_to_extract_DTP_oxPPP_KO, :);  % All samples, only NADPH-consuming reactions

writetable(nadph_flux_table_DTP_oxPPP_KO, 'nadph_consuming_flux_samples_DTP_oxPPP_KO.csv');


% Plotting NADPH Consuming Fluxes for DTP WT and KO 

common_rxns = intersect(nadph_flux_table_DTP.Reaction, ...
                        nadph_flux_table_DTP_oxPPP_KO.Reaction, 'stable');

% Align both tables to the same reaction order
[~, WT_lookup] = ismember(common_rxns, nadph_flux_table_DTP.Reaction);
[~, KO_lookup] = ismember(common_rxns, nadph_flux_table_DTP_oxPPP_KO.Reaction);

WT_flux = table2array(nadph_flux_table_DTP(WT_lookup, 2:end));   % nRxns x nSamples
KO_flux = table2array(nadph_flux_table_DTP_oxPPP_KO(KO_lookup, 2:end));
rxnLabels = nadph_flux_table_DTP.Reaction(WT_lookup);            % newly aligned

% For plotting, transpose so MATLAB violins plot samples as rows:
WT_flux = WT_flux';   % [nSamples x nRxns]
KO_flux = KO_flux';   % [nSamples x nRxns]

% Top 10 by mean WT flux
meanWT = mean(WT_flux, 1);
[~, sortIdx] = sort(meanWT, 'descend');
topN = 10;
topIdx = sortIdx(1:topN);

WT_top = WT_flux(:, topIdx);
KO_top = KO_flux(:, topIdx);
labels_top = rxnLabels(topIdx);

% Exclude GTHOr from the top 10 nadph consuming fluxes
% Since we were interested in non-antioxidant nadph-consuming fluxes
exclude_rxns = {'GTHOr'};
exclude_idx = ismember(labels_top, exclude_rxns);

labels_top_noGTHOr = labels_top(~exclude_idx);
WT_top_noGTHOr = WT_top(:, ~exclude_idx);
KO_top_noGTHOr = KO_top(:, ~exclude_idx);

topN = size(WT_top_noGTHOr, 2);

% Build interleaved matrix and color/style arrays
interleaved_data = zeros(size(WT_top_noGTHOr,1), 2*topN);
group_labels = cell(1, 2*topN);
colors      = cell(1, 2*topN);

for i = 1:topN
    interleaved_data(:, 2*i-1) = WT_top_noGTHOr(:, i);           % WT
    interleaved_data(:, 2*i)   = KO_top_noGTHOr(:, i);           % KO
    group_labels{2*i-1} = [labels_top_noGTHOr{i} ' WT'];
    group_labels{2*i}   = [labels_top_noGTHOr{i} ' KO'];
    colors{2*i-1} = [0 1 1];        % Cyan for DTP WT
    colors{2*i}   = [1 0.65 0];     % Orange for DTP oxPPP KO
end

% Remove constant violins (all identical values)
goodCols = arrayfun(@(i) numel(unique(interleaved_data(:,i))) > 1, 1:size(interleaved_data,2));
interleaved_data_clean = interleaved_data(:, goodCols);
group_labels_clean     = group_labels(goodCols);
colors_clean           = colors(goodCols);

% Plot
figure;
v = violinplot(interleaved_data_clean, group_labels_clean, 'ShowData', false);

for i = 1:length(colors_clean)
    if isprop(v(i),'ViolinColor')
        v(i).ViolinColor = colors_clean{i};
    end
end
set(gca, 'YScale', 'log');
xtickangle(45);
ylabel('Flux (nmol/mg/hr)');
title('Top NADPH-Consuming Reaction Fluxes');
set(gca, 'FontSize', 14);

% Only one tick per reaction, placed between each pair (centered)
nReactions = numel(labels_top_noGTHOr);         % or num labels after exclusions
xticks = 1.5:2:2*nReactions;
set(gca, 'XTick', xticks, 'XTickLabel', labels_top_noGTHOr);

hold on;
hWT = plot(NaN,NaN,'o','MarkerFaceColor',[0 1 1],'MarkerEdgeColor','none','DisplayName','DTP WT');
hKO = plot(NaN,NaN,'o','MarkerFaceColor',[1 0.65 0],'MarkerEdgeColor','none','DisplayName','DTP oxPPP KO');
legend([hWT, hKO],'DTP WT','DTP oxPPP KO','Location','best');
hold off;


%% Identifying, extracting, and plotting NADPH producing reactions
% Assuming that flux sampling was performed on DTP and DTP_oxPPP_KO groups
% and that flux_table_DTP and flux_table_DTP_oxPPP_KO are generated (view
% previous sections)

% Identify nadph producing reactions
nadph_met_idx = find(strcmp(model.mets, 'nadph_c'));    % index of NADPH
nadph_stoich_row = model.S(nadph_met_idx, :);            % stoichiometry vector for NADPH across all rxns

% NADPH-producing reactions: positive coefficient (NADPH appears as a product)
nadph_producing_rxn_idxs = find(nadph_stoich_row > 1e-8);   % threshold to avoid rounding errors

nadph_producing_rxns = model.rxns(nadph_producing_rxn_idxs);

% Extract nadph-producing fluxes from labeled dataframe

% Find rows (reactions) matching NADPH producers
rows_to_extract_producer_DTP = ismember(flux_table_DTP.Reaction, nadph_producing_rxns);

nadph_prod_flux_table_DTP = flux_table_DTP(rows_to_extract_producer_DTP, :);

writetable(nadph_prod_flux_table_DTP, 'nadph_producing_flux_samples_DTP.csv');

% repeat for other group

rows_to_extract_producer_DTP_oxPPP_KO = ismember(flux_table_DTP_oxPPP_KO.Reaction, nadph_producing_rxns);

nadph_prod_flux_table_DTP_oxPPP_KO = flux_table_DTP_oxPPP_KO(rows_to_extract_producer_DTP_oxPPP_KO, :);

writetable(nadph_prod_flux_table_DTP_oxPPP_KO, 'nadph_producing_flux_samples_DTP_oxPPP_KO.csv');


% Plotting NADPH Producing Fluxes for DTP WT and KO 

% Find intersection/order for both tables
common_rxns = intersect(nadph_prod_flux_table_DTP.Reaction, ...
                        nadph_prod_flux_table_DTP_oxPPP_KO.Reaction, 'stable');

% Align both tables to the same reaction order
[~, WT_lookup] = ismember(common_rxns, nadph_prod_flux_table_DTP.Reaction);
[~, KO_lookup] = ismember(common_rxns, nadph_prod_flux_table_DTP_oxPPP_KO.Reaction);

WT_flux = table2array(nadph_prod_flux_table_DTP(WT_lookup, 2:end));   % nRxns x nSamples
KO_flux = table2array(nadph_prod_flux_table_DTP_oxPPP_KO(KO_lookup, 2:end));
rxnLabels = nadph_prod_flux_table_DTP.Reaction(WT_lookup);            % newly aligned

% For plotting, transpose so MATLAB violins plot samples as rows:
WT_flux = WT_flux';   % [nSamples x nRxns]
KO_flux = KO_flux';

% Top 10 by mean WT flux
meanWT = mean(WT_flux, 1);
[~, sortIdx] = sort(meanWT, 'descend');
topN = 10;
topIdx = sortIdx(1:topN);

WT_top = WT_flux(:, topIdx);
KO_top = KO_flux(:, topIdx);
labels_top = rxnLabels(topIdx);

% exclude G6PDH2r, GND, and G3PD2 from the top 10 nadph producing fluxes
% we were interested in nadph producing fluxes other than the pentose
% phosphate pathway and G3PD2
exclude_rxns = {'G6PDH2r', 'GND', 'G3PD2'};
exclude_idx = ismember(labels_top, exclude_rxns);

labels_top_noGG = labels_top(~exclude_idx);
WT_top_noGG = WT_top(:, ~exclude_idx);
KO_top_noGG = KO_top(:, ~exclude_idx);

topN = size(WT_top_noGG, 2);

% Build interleaved matrix and color/style arrays
interleaved_data = zeros(size(WT_top_noGG, 1), 2 * topN);
group_labels = cell(1, 2 * topN);
colors      = cell(1, 2 * topN);

for i = 1:topN
    interleaved_data(:, 2*i-1) = WT_top_noGG(:, i);     % WT
    interleaved_data(:, 2*i)   = KO_top_noGG(:, i);     % KO
    group_labels{2*i-1} = [labels_top_noGG{i} ' WT'];
    group_labels{2*i}   = [labels_top_noGG{i} ' KO'];
    colors{2*i-1} = [0 1 1];        % Cyan for WT
    colors{2*i}   = [1 0.65 0];     % Orange for KO
end

% Remove constant violins (all identical values)
goodCols = arrayfun(@(i) numel(unique(interleaved_data(:, i))) > 1, 1:size(interleaved_data, 2));
interleaved_data_clean = interleaved_data(:, goodCols);
group_labels_clean     = group_labels(goodCols);
colors_clean           = colors(goodCols);

% cap negative values at 1 to plot on log scale
interleaved_data_capped = interleaved_data_clean;
interleaved_data_capped(interleaved_data_capped < 1) = 1;

% Plot
figure;
ax = gca;
ax.ColorOrder = repmat([0 1 1; 1 0.65 0], ceil(topN / 2), 1);
ax.ColorOrderIndex = 1;

v = violinplot(interleaved_data_capped, group_labels_clean, 'showData', false);

for i = 1:length(colors_clean)
    if isprop(v(i), 'ViolinColor')
        v(i).ViolinColor = colors_clean{i};
    end
end

set(gca, 'YScale', 'log');
xtickangle(45);
ylabel('Flux (nmol/mg/hr)');
title('Top NADPH-Producing Reaction Fluxes');
set(gca, 'FontSize', 14);

% ----- Single reaction labels: One per pair, centered -----
Nxticks = numel(labels_top_noGG);
xtick_positions = 1.5:2:(2*Nxticks);   % centers of each WT/KO pair
set(gca, 'XTick', xtick_positions, 'XTickLabel', labels_top_noGG);

% ----- Manual legend -----
hold on;
hWT = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0 1 1], 'MarkerEdgeColor', 'none', 'DisplayName', 'DTP WT');
hKO = plot(NaN, NaN, 'o', 'MarkerFaceColor', [1 0.65 0], 'MarkerEdgeColor', 'none', 'DisplayName', 'DTP oxPPP KO');
legend([hWT, hKO], 'DTP WT', 'DTP oxPPP KO', 'Location', 'best');
hold off;


%% Plotting fluxes of multiple reactions in a grouped violin plot 

rxnID_list = {'GTHOr', 'GTHOm','GTHPi','CATp','GTHRDt', 'GTHPm','TRDR'}; % your selected reactions
nRxns = numel(rxnID_list);
FluxList = [];
GroupList = [];
ReactionList = [];

for i = 1:nRxns
    rxnIdx = find(strcmpi(model.rxns, rxnID_list{i}));
    Parental_flux = samples_Parental(rxnIdx, :)';
    DTP_flux = samples_DTP(rxnIdx, :)';
    FluxList = [FluxList; Parental_flux;DTP_flux ];
    GroupList = [GroupList; repmat({'Parental'}, length(Parental_flux), 1); ...
                           repmat({'DTP'}, length(DTP_flux), 1)];
    ReactionList = [ReactionList; repmat(rxnID_list(i), length(Parental_flux) + length(DTP_flux), 1)];
end

GroupLabelList = strcat(ReactionList, " | ", GroupList); % e.g. 'GTHOr | Parental'
% Ensure group_labels categorical is in the order of rxnID_list:
ordered_labels = {};
for i = 1:numel(rxnID_list)
    ordered_labels{end+1} = [rxnID_list{i} ' | Parental'];
    ordered_labels{end+1} = [rxnID_list{i} ' | DTP'];
end
ordered_labels = string(ordered_labels);
group_labels = categorical(GroupLabelList, ordered_labels, 'Ordinal', true);

% Set up alternating colors: black for Parental, red for DTP
figure;
ax = gca;
ax.ColorOrder = repmat([0 0 0; 1 0 0], ceil(nRxns/2), 1);
ax.ColorOrderIndex = 1;
violinplot(FluxList, group_labels, 'ShowData', false);

set(gca, 'YScale', 'log');
xtickangle(45);
ylabel('Flux (nmol/mg/hr)');
title('Antioxidant Reaction Fluxes');
set(gca, 'FontSize', 14);

% Single reaction labels
xticks = 1.5:2:2*nRxns;
set(gca, 'XTick', xticks, 'XTickLabel', rxnID_list);

hold on;
hParental = plot(NaN,NaN,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','DisplayName','Parental');
hDTP = plot(NaN,NaN,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','none','DisplayName','cPACC10');
legend([hParental, hDTP],'Parental','DTP','Location','best');
hold off;


%% Generating mean, SEM, 95 percent CI statistics 

% Ensure orientation
if size(samples_DTP,2) ~= numel(model.rxns)
    samples_DTP = samples_DTP'; % Should become [nSamples x nRxns]
end
if size(samples_DTP_oxPPP_KO,2) ~= numel(model.rxns)
    samples_DTP_oxPPP_KO = samples_DTP_oxPPP_KO';
end

rxnNames = model.rxns;

% Statistics 
nRxns = numel(rxnNames);

% WT group
mean_WT  = mean(samples_DTP,1);
sem_WT   = std(samples_DTP,0,1) ./ sqrt(size(samples_DTP,1));
ci_low_WT  = mean_WT - 1.96 * sem_WT;
ci_high_WT = mean_WT + 1.96 * sem_WT;
percentile_low_WT  = prctile(samples_DTP, 2.5, 1);
percentile_high_WT = prctile(samples_DTP, 97.5,1);

% KO group
mean_KO  = mean(samples_DTP_oxPPP_KO,1);
sem_KO   = std(samples_DTP_oxPPP_KO,0,1) ./ sqrt(size(samples_DTP_oxPPP_KO,1));
ci_low_KO  = mean_KO - 1.96 * sem_KO;
ci_high_KO = mean_KO + 1.96 * sem_KO;
percentile_low_KO  = prctile(samples_DTP_oxPPP_KO, 2.5, 1);
percentile_high_KO = prctile(samples_DTP_oxPPP_KO, 97.5,1);

% Combine into Table 
T = table(rxnNames(:), ...
    mean_WT(:), sem_WT(:), ci_low_WT(:), ci_high_WT(:), percentile_low_WT(:), percentile_high_WT(:), ...
    mean_KO(:), sem_KO(:), ci_low_KO(:), ci_high_KO(:), percentile_low_KO(:), percentile_high_KO(:), ...
    'VariableNames', {'Reaction', ...
    'WT_Mean', 'WT_SEM', 'WT_Norm_CI_Low', 'WT_Norm_CI_High', 'WT_Perc_CI_Low', 'WT_Perc_CI_High', ...
    'KO_Mean', 'KO_SEM', 'KO_Norm_CI_Low', 'KO_Norm_CI_High', 'KO_Perc_CI_Low', 'KO_Perc_CI_High'});

% Export as CSV 
writetable(T, 'flux_statistics_DTP_WT_vs_KO.csv');





