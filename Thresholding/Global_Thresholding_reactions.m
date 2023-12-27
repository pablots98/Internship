%%%%%%%%%%%%%%%%%%%% Global thresholding reactions %%%%%%%%%%%%%%%%%%%%%%%%

%% Setup and initialization
gurobi_setup;
initCobraToolbox(false);
setenv('GUROBI_PATH', 'C:\gurobi1001\win64\matlab');
changeCobraSolver('gurobi', 'all');

%% Load the data
% Load the transcriptomics data, housekeeping genes, and the metabolic model
data = readtable('Merged_data.xlsx');
h_k_g = readtable('NM2ENSG.xlsx');
model = load('Human-GEM_Cobra_v1.01.mat');
model = model.model;
col_names = data.Properties.VariableNames(7:end);

%% Pre processing
% Original dataset
data_log = data(:, 1:end);
genes = data_log(:, 2);
Ensembl_id = data_log(:, 1);
data_to_log = table2array(data_log(:, 7:end));
logged_data = log10(data_to_log + 1);
log_data = [Ensembl_id, genes, array2table(logged_data)];
log_data.Properties.VariableNames(3:end) = col_names;
log_data.Properties.VariableNames{2} = 'gene';
log_data.Properties.VariableNames{1} = 'Ensembl_GeneID';

% Create the values columns
numSamples = width(log_data(:, 2:end)) - 1;
allGenes = {};
sampleNames = log_data.Properties.VariableNames(2:end);
geneExpressionMatrix = [];

for i = 1:numSamples
    log_data.value = log_data{:, i + 1} - min(log_data{:, i + 1});
    [geneList, geneExpression] = findUsedGenesLevels(model, log_data);
    allGenes = unique([allGenes; geneList]);
    tempMatrix = zeros(length(allGenes), 1);
    [~, geneIdx] = ismember(geneList, allGenes);
    tempMatrix(geneIdx) = geneExpression;
    geneExpressionMatrix = [geneExpressionMatrix, tempMatrix];
    log_data.value = [];
end

metabolic_genes = array2table(geneExpressionMatrix, 'RowNames', allGenes, 'VariableNames', sampleNames);

%% Calculate the mean and standard deviation form different samples
stat_results = calculateGeneExpressionStats(metabolic_genes);

%% Take the gene expression from the metabolic genes
met_genes = metabolic_genes.Properties.RowNames;
index_names = ismember(log_data.Ensembl_GeneID, met_genes);
data_met = log_data(index_names, :);

%% Global thresholding
results = table;
expression_col = metabolic_genes.Properties.VariableNames(1:end);

for i = 1:length(expression_col)
    col_name = expression_col{i};
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

%% From core-genes to core-reactions for each sample

% STEP 1: Prepare Expression Data for Each Sample
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

% STEP 2: Map Expression to Reactions for Each Sample
mappedReactions = cell(numSamples, 1);

for i = 1:numSamples
    expressionData = expressionDataSamples{i};
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model, expressionData);
    mappedReactions{i} = expressionRxns;
end

% STEP 3: Identify Reactions Present in Each Sample
react_names = model.rxns;
reactionNamesPerSample = cell(numSamples, 1);

for i = 1:numSamples
    currentReactions = mappedReactions{i};
    currentReactionNames = {};
    for j = 1:numel(currentReactions)
        if ~isnan(currentReactions(j))
            currentReactionNames{end+1} = react_names{j};
        end
    end
    reactionNamesPerSample{i} = currentReactionNames;
end

%% Housekeeping genes
genes_table = data_met.Ensembl_GeneID;
index_names = ismember(genes_table, h_k_g.converted_alias);
hkg_met = data_met(index_names, :);
hkg_met_ens = hkg_met(:, "Ensembl_GeneID");
geneIDs = table2cell(hkg_met_ens);

% Obtain the core reactions from the core genes
[results] = findRxnsFromGenes(model, geneIDs);

% Obtain just the housekeeping reactions names
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

% Delet the repeated names
housekeep_react_unique = unique(housekeep_react);
disp(housekeep_react_unique);

%% Compare the amount of housekeeping core reactions

% STEP 1: Prepare the data structure
numSamples = numel(reactionNamesPerSample);
housekep_core_react = struct('numHousekeepingCoreReactions', [], 'housekeepingCoreReactions', [], 'percentage', []);
totalHousekeepingReactions = numel(housekeep_react_unique);

% STEP 2: Do the comparison and the calculus
for i = 1:numSamples
    coreReactions = reactionNamesPerSample{i};
    housekeepingInCore = ismember(housekeep_react_unique, coreReactions);
    housekep_core_react(i).numHousekeepingCoreReactions = sum(housekeepingInCore);
    housekep_core_react(i).housekeepingCoreReactions = housekeep_react_unique(housekeepingInCore);
    housekep_core_react(i).percentage = (housekep_core_react(i).numHousekeepingCoreReactions / totalHousekeepingReactions) * 100;
end

