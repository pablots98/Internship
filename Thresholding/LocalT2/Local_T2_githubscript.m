%%%%%%%%%%%%%%%%% Local T2 thresholding reactions %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Gurobi setup and Cobra Toolbox initialization
gurobi_setup;
initCobraToolbox(false);
setenv('GUROBI_PATH', 'C:\gurobi1001\win64\matlab');
changeCobraSolver('gurobi', 'all');

%% Load data
% Load transcriptomics data, housekeeping genes, and metabolic model
data = readtable('Merged_data.xlsx'); % Transcriptomics data
h_k_g = readtable('NM2ENSG.xlsx');    % Housekeeping genes with Ensembl IDs
model = load('Human-GEM_Cobra_v1.01.mat'); % Human1 metabolic model
model = model.model;

%ensure the model does not contain blocked reactions
[fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool, ~, fluxConsistModel] = findFluxConsistentSubset(model);
model = fluxConsistModel
model

%% Preprocessing
genes = data(:, 2);
Ensembl_id = data(:, 1);
data_to_log = table2array(data(:, 7:end));
logged_data = log10(data_to_log + 1);
log_data = [Ensembl_id, genes, array2table(logged_data)];
col_names = data.Properties.VariableNames(7:end);
log_data.Properties.VariableNames(3:end) = col_names;
log_data.Properties.VariableNames{2} = 'gene';
log_data.Properties.VariableNames{1} = 'Ensembl_GeneID';

%% Processing the gene expression dataset
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

%% Taking gene expression from the metabolic genes
met_genes = metabolic_genes.Properties.RowNames;
index_names = ismember(log_data.Ensembl_GeneID, met_genes);
data_met = log_data(index_names, :);
logdata = data_met(:, 3:end);
gene_names = data_met{:, 1};

%% Local thresholding
up_percentage = 75;
low_percentage = 25;
[expression_scoreMatrix, coreMat] = localT2_new(logdata, low_percentage, up_percentage);

%% Map to Expression

% Adjusting the structure to what is likely expected by the function
expressionData.gene = gene_names; % Assuming gene_names is a cell array of gene identifiers
expressionData.value = expression_scoreMatrix; % Assuming this is a numeric matrix of expression values

Rxns_local25_75 = [];
geneUsed_local25_75 = {};

% Iterate over each sample
for i = 1:length(sampleNames) % Assuming sampleNames is a list of your samples
    % Extract the expression data for the current sample
    expressionDataSample = struct();
    expressionDataSample.gene = expressionData.gene;
    expressionDataSample.value = expressionData.value(:, i); % Selecting the column for the current sample

    % Map expression data to reactions for the current sample
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model, expressionDataSample, 'false');
    Rxns_local25_75 = [Rxns_local25_75, expressionRxns]; % Store the mapped reactions
    geneUsed_local25_75{i} = gene_used; % Store the genes used in the mapping
end

%% Check core reactions

for i = 1:length(sampleNames)
    ActiveRxns_local25_75{i} = (find(Rxns_local25_75(:, i)>=1));
end
disp(ActiveRxns_local25_75)