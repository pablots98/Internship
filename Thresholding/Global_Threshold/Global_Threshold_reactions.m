%%%%%%%%%%%%%%%%%%%% Global thresholding reactions %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize the COBRA Toolbox and set the solver to Gurobi
gurobi_setup;
initCobraToolbox(false);
setenv('GUROBI_PATH', 'C:\gurobi1001\win64\matlab');
changeCobraSolver('gurobi', 'all');

%%                              Load data                                %%
% Load transcriptomics data, housekeeping genes, and metabolic model
data = readtable('Merged_data.xlsx'); % Transcriptomics data
h_k_g = readtable('NM2ENSG.xlsx');    % Housekeeping genes with Ensembl IDs
model = load('Human-GEM_Cobra_v1.01.mat'); % Human1 metabolic model
model = model.model;

%%          Ensure the model does not contain blocked reactions          %%
% [fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool, ~, fluxConsistModel] = findFluxConsistentSubset(model);
% model = fluxConsistModel; % load the new model with no blocked reactions
% save('model', 'model');
model = load('model.mat');
model = model.model;
model_genes = model.genes; % ENSEMBL_IDs of the genes in the model
%%                              Preprocessing                            %%
% Colect the data needed to create the table
Ensembl_id = data(:, 1);
sampleNames = data.Properties.VariableNames(7:end); 
data_to_log = table2array(data(:, 7:end));
logged_data = log10(data_to_log + 1); % Normalize the data (+1 to avoid having 0 values) 

%Create the new table with the data obtained before
log_data = [Ensembl_id, array2table(logged_data)]; % Load the data on the new table  
log_data.Properties.VariableNames(2:end) = sampleNames; % Variables names 
log_data.Properties.VariableNames{1} = 'gene'; % Change ENSEML_ID name to gene (To run findUsedGenesLevels)

%%              Processing the gene expression dataset                   %% 
% Define new cells and new Matrixes
allGenes = {};
geneExpressionMatrix = [];


for i = 1: width(sampleNames)
    log_data.value = log_data{:, i + 1}; %- min(log_data{:, i + 1}); % shift minimum to 0 ¡¡ ASK MARIAN HOW HE WANT TO NORMALIZE IT!! ACTUALLY IT SAYS FPKM, SO I ALREADY HAVE IT
    [geneList, geneExpression] = findUsedGenesLevels(model, log_data); % Athough is just for the last sample, the genes present in the models are always the same
    geneExpressionMatrix = [geneExpressionMatrix, geneExpression'];
    log_data.value = [];
end

metabolic_genes = array2table(geneExpressionMatrix, 'RowNames', geneList, 'VariableNames', sampleNames);

% Delete the NaN values
metabolic_genes = rmmissing(metabolic_genes, 'MinNumMissing', size(metabolic_genes, 2)); % So we can use the local thresholding 

%%      Global thresholding to determine core and non-core genes         %%
results = table;

for i = 1:length(sampleNames)
    col_name = sampleNames{i};
    expData = metabolic_genes(:, 1);
    expData(:, 2) = metabolic_genes(:, col_name);
    expData.Properties.VariableNames{1} = 'gene';
    expData.Properties.VariableNames{2} = 'ExpressionValue';
    expData.LogExpressionValue = log10(expData{:, 2} + 1);
    expData.Value = expData.LogExpressionValue - min(expData.LogExpressionValue);
    threshold = prctile(expData.Value, 75);
    coreGenesIdx = expData.Value >= threshold;
    coreGenes = expData(coreGenesIdx, :);
    coreGeneIDs = coreGenes.Properties.RowNames;
    results{i, 'Sample'} = {col_name};
    results{i, 'CoreGeneIDs'} = {coreGeneIDs};
end

results_CoreGenes = results;
results_CoreGenes(1, :) = [];
disp(results_CoreGenes);

%%          Map core genes to core reactions for each sample             %%
numSamples = height(results_CoreGenes);
expressionDataSamples = cell(numSamples, 1);

for i = 1:numSamples
    coreGenes = results_CoreGenes.CoreGeneIDs{i};
    [~, idx] = ismember(coreGenes, metabolic_genes.Properties.RowNames);
    expressionData = struct();
    expressionData.gene = metabolic_genes.Properties.RowNames(idx);
    expressionData.value = table2array(metabolic_genes(idx, 2:end));
    expressionDataSamples{i} = expressionData;
end

mappedReactions = cell(numSamples, 1);

for i = 1:numSamples
    expressionData = expressionDataSamples{i};
    [expressionRxns, ~, ~] = mapExpressionToReactions(model, expressionData);
    mappedReactions{i} = expressionRxns;
end

% Check for repeated reaction names
react_names = model.rxns;
numSamp = numel(mappedReactions);
reactionNamesPerSample = cell(numSamp, 1);

for i = 1:numSamp
    currentReactions = mappedReactions{i};
    currentReactionNames = {};
    for j = 1:numel(currentReactions)
        if ~isnan(currentReactions(j))
            currentReactionNames{end+1} = react_names{j};
        end
    end
    reactionNamesPerSample{i} = currentReactionNames;
end

% Check and display unique reaction names
allNames = {};
for i = 1:length(reactionNamesPerSample)
    allNames = [allNames; reactionNamesPerSample{i}(:)];
end

uniqueNames = unique(allNames);
if length(uniqueNames) < length(allNames)
    disp('There are repeated names');
else
    disp('No repeated names');
end

%%                  Housekeeping genes analysis                          %% (CHECK, I THINK THIS IS WRONG)
gene_names = metabolic_genes.Properties.RowNames;
index_names = ismember(gene_names, h_k_g.converted_alias);
hkg_met = metabolic_genes(index_names, :);
hkg_met_ens = hkg_met.Properties.RowNames;
%geneIDs = table2cell(hkg_met_ens);
[results] = findRxnsFromGenes(model, hkg_met_ens);

% Extract unique housekeeping reaction names
fields = fieldnames(results);
housekeep_react = {};
for i = 1:length(fields)
    field = fields{i};
    cellArray = results.(field);
    for j = 1:size(cellArray, 1)
        firstcol = cellArray{j, 1};
        housekeep_react{end+1, 1} = firstcol;
    end
end
housekeep_react_unique = unique(housekeep_react);
disp(housekeep_react_unique);

%% Compare the number of housekeeping core genes
numSample_genes = numel(results_CoreGenes.CoreGeneIDs);
housekep_core_gene_G = struct("numHousekeepingCoreGenes", [], 'housekeepingCoreGenes', [], 'percentage', []);
totalHousekeepingGenes = numel(hkg_met_ens);

for i = 1:numSample_genes
    coreGenes = results_CoreGenes.CoreGeneIDs{i};
    housekeepingGeneNames = hkg_met_ens;
    housekeepingInCore_g = ismember(housekeepingGeneNames, coreGenes);
    housekep_core_gene_G(i).numHousekeepingCoreGenes = sum(housekeepingInCore_g);
    housekep_core_gene_G(i).housekeepingCoreGenes = housekeepingGeneNames(housekeepingInCore_g);
    housekep_core_gene_G(i).percentage = (housekep_core_gene_G(i).numHousekeepingCoreGenes / totalHousekeepingGenes) * 100;
end

%% Compare the number of housekeeping core reactions
numSample = numel(reactionNamesPerSample);
housekep_core_react_G = struct('numHousekeepingCoreReactions', [], 'housekeepingCoreReactions', [], 'percentage', []);
totalHousekeepingReactions = numel(housekeep_react_unique);

for i = 1:numSample
    coreReactions = reactionNamesPerSample{i};
    housekeepingInCore = ismember(housekeep_react_unique, coreReactions);
    housekep_core_react_G(i).numHousekeepingCoreReactions = sum(housekeepingInCore);
    housekep_core_react_G(i).housekeepingCoreReactions = housekeep_react_unique(housekeepingInCore);
    housekep_core_react_G(i).percentage = (housekep_core_react_G(i).numHousekeepingCoreReactions / totalHousekeepingReactions) * 100;
end

