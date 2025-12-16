% Startup
addpath(genpath('C:\COBRA'))
savepath
initCobraToolbox(false)
changeCobraSolver('gurobi','all')

% Load full model (after initializing your COBRA session and solver)
model = readCbModel('Recon3D.mat');

% Load proteomics data
protT = readtable('metabolic_only_proteomics_normalized_abundances.csv');
expressedEntrezIDs = string(protT.ENTREZID);

% List MFA/13C flux reaction IDs (from your flux file)
mfaT = readtable('Parental_flux_results_for_GEM.csv');
mfa_rxns = unique(string(mfaT.reaction_id));  % string array of MFA reaction IDs

% Extract EntrezIDs from model.genes
modelEntrez = cellfun(@(g) strtok(g, '_'), model.genes, 'UniformOutput', false); % returns cell array of EntrezIDs as strings

% Identify reactions to keep (initialize once!)
nRxns = numel(model.rxns);
keepRxns = false(nRxns, 1);

for i = 1:nRxns
    r_id = string(model.rxns{i});
    if any(r_id == mfa_rxns)
        keepRxns(i) = true;
        continue
    end
    rule_str = model.rules{i};
    if isempty(rule_str)
        continue
    end
    idx_nums = regexp(rule_str, 'x\((\d+)\)', 'tokens');
    if ~isempty(idx_nums)
        idxs = unique(str2double([idx_nums{:}]));
        % Instead of model.genes, use stripped Entrez IDs:
        geneEntrez = string(modelEntrez(idxs));
        if any(ismember(geneEntrez, expressedEntrezIDs))
            keepRxns(i) = true;
        end
    end
end

% Keep exchange, transport, and objectives 
isExchange  = reshape(startsWith(model.rxns, "EX_"), [], 1);
isTransport = reshape(contains(lower(model.rxnNames), 'transport'), [], 1);
isObj       = reshape(strcmpi(model.rxns, 'GTHOr') | strcmpi(model.rxns, 'TRDR'), [], 1);

% Generate reduced model
modelReduced = removeRxns(model, model.rxns(~keepRxns));

% Clean up dead ends (optional but recommended)
if exist('findBlockedRxns','file')
    deadRxns = findBlockedRxns(modelReduced);  % Should be cell or string array
elseif exist('detectDeadEnds','file')
    deadMets = detectDeadEnds(modelReduced);
    if isnumeric(deadMets)
        metsToLoop = modelReduced.mets(deadMets); % get metabolite IDs
    elseif iscell(deadMets)
        metsToLoop = deadMets;
    elseif isstring(deadMets)
        metsToLoop = cellstr(deadMets);
    else
        error('Unexpected deadMets type');
    end
    deadRxns = [];
    for j = 1:length(metsToLoop)
        rxnIDs = findRxnsFromMets(modelReduced, metsToLoop{j});  % use {} for cell array of char
        deadRxns = [deadRxns; rxnIDs(:)];
    end
    deadRxns = unique(deadRxns);
else
    error('No dead-end or blocked reaction finder available in COBRA Toolbox.');
end

modelReduced = removeRxns(modelReduced, deadRxns);

% Set the objective
modelReduced = changeObjective(modelReduced, 'GTHOr'); % or other objective rxn ID

% Test if reduced model is solvable
FBAsolution = optimizeCbModel(modelReduced);
disp(FBAsolution.f);

% Save reduced model
save('Recon3D_reduced.mat', 'modelReduced');