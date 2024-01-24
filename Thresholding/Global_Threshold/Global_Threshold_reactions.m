%%%%%%%%%%%%%%%%%%%% Global thresholding reactions %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize the COBRA Toolbox and set the solver to Gurobi
gurobi_setup;
initCobraToolbox(false);
setenv('GUROBI_PATH', 'C:\gurobi1001\win64\matlab');
changeCobraSolver('gurobi', 'all');

%% Load necessary data
data = readtable('Merged_data.xlsx'); % Transcriptomics data
h_k_g = readtable('NM2ENSG.xlsx'); % Housekeeping genes with Ensembl IDs
model = load('Human-GEM_Cobra_v1.01.mat'); % Human1 metabolic model
model = model.model;
%First, in all model extraction it is good to pre_process the model to
%ensure it does not contain blocked reactions

[fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool, ~, fluxConsistModel] = findFluxConsistentSubset(model);
model = fluxConsistModel
model
%% Preprocess the data
genes = data(:, 2); % Keep the gene column
Ensembl_id = data(:, 1); % Keep the Ensembl ID column

% Normalize the expression data (log10(x + 1))
data_to_log = table2array(data(:, 7:end));
logged_data = log10(data_to_log + 1);
log_data = [Ensembl_id, genes, array2table(logged_data)];
col_names = data.Properties.VariableNames(7:end);
log_data.Properties.VariableNames(3:end) = col_names;
log_data.Properties.VariableNames{1} = 'gene';
log_data.Properties.VariableNames{2} = 'Ensembl_GeneID';

%% Process gene expression dataset
num_samples = width(log_data(:, 3:end));
allGenes = {};
sampleNames = log_data.Properties.VariableNames(3:end);
geneExpressionMatrix = [];

for i = 1:num_samples
    log_data.value = log_data{:, i + 2} - min(log_data{:, i + 2});
    [geneList, geneExpression] = findUsedGenesLevels(model, log_data);
    allGenes = unique([allGenes; geneList]);
    tempMatrix = zeros(length(allGenes), 1);
    [~, geneIdx] = ismember(geneList, allGenes);
    tempMatrix(geneIdx) = geneExpression;
    geneExpressionMatrix = [geneExpressionMatrix, tempMatrix];
    log_data.value = [];
end

metabolic_genes = array2table(geneExpressionMatrix, 'RowNames', allGenes, 'VariableNames', sampleNames);

% Take the metabolic genes expression data
met_genes = metabolic_genes.Properties.RowNames
% Find indexes of the genes that are present in the metabolic model. 
index_names = ismember(log_data.Ensembl_GeneID, met_genes);
% Take the rows of the dataset that match the indexes obtained before, with all the data 
data_met = log_data(index_names, :);

%% Global thresholding to determine core and non-core genes
results = table;

for i = 1:length(col_names)
    col_name = col_names{i};
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

%% Map core genes to core reactions for each sample
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

%% Housekeeping genes analysis (CHECK, I THINK THIS IS WRONG)
genes_table = data_met.Ensembl_GeneID;
index_names = ismember(genes_table, h_k_g.converted_alias);
hkg_met = data_met(index_names, :);
hkg_met_ens = hkg_met(:, "Ensembl_GeneID");
geneIDs = table2cell(hkg_met_ens);
[results] = findRxnsFromGenes(model, geneIDs);

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
housekep_core_gene = struct("numHousekeepingCoreGenes", [], 'housekeepingCoreGenes', [], 'percentage', []);
totalHousekeepingGenes = numel(hkg_met.Ensembl_GeneID);

for i = 1:numSample_genes
    coreGenes = results_CoreGenes.CoreGeneIDs{i};
    housekeepingGeneNames = hkg_met.Ensembl_GeneID;
    housekeepingInCore_g = ismember(housekeepingGeneNames, coreGenes);
    housekep_core_gene(i).numHousekeepingCoreGenes = sum(housekeepingInCore_g);
    housekep_core_gene(i).housekeepingCoreGenes = housekeepingGeneNames(housekeepingInCore_g);
    housekep_core_gene(i).percentage = (housekep_core_gene(i).numHousekeepingCoreGenes / totalHousekeepingGenes) * 100;
end

%% Compare the number of housekeeping core reactions
numSample = numel(reactionNamesPerSample);
housekep_core_react = struct('numHousekeepingCoreReactions', [], 'housekeepingCoreReactions', [], 'percentage', []);
totalHousekeepingReactions = numel(housekeep_react_unique);

for i = 1:numSample
    coreReactions = reactionNamesPerSample{i};
    housekeepingInCore = ismember(housekeep_react_unique, coreReactions);
    housekep_core_react(i).numHousekeepingCoreReactions = sum(housekeepingInCore);
    housekep_core_react(i).housekeepingCoreReactions = housekeep_react_unique(housekeepingInCore);
    housekep_core_react(i).percentage = (housekep_core_react(i).numHousekeepingCoreReactions / totalHousekeepingReactions) * 100;
end

