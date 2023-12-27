% Local T2 thresholding reactions

% Gurobi setup and Cobra Toolbox initialization
gurobi_setup;
initCobraToolbox(false);
setenv('GUROBI_PATH', 'C:\gurobi1001\win64\matlab');
changeCobraSolver('gurobi', 'all');

% Load data
% Load transcriptomics data, housekeeping genes, and metabolic model
data = readtable('Merged_data.xlsx'); % Transcriptomics data
h_k_g = readtable('NM2ENSG.xlsx');    % Housekeeping genes with Ensembl IDs
model = load('Human-GEM_Cobra_v1.01.mat'); % Human1 metabolic model
model = model.model;

% Preprocessing
genes = data(:, 2);
Ensembl_id = data(:, 1);
data_to_log = table2array(data(:, 7:end));
logged_data = log10(data_to_log + 1);
log_data = [Ensembl_id, genes, array2table(logged_data)];
col_names = data.Properties.VariableNames(7:end);
log_data.Properties.VariableNames(3:end) = col_names;
log_data.Properties.VariableNames{2} = 'gene';
log_data.Properties.VariableNames{1} = 'Ensembl_GeneID';

% Processing the gene expression dataset
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

% Taking gene expression from the metabolic genes
met_genes = metabolic_genes.Properties.RowNames;
index_names = ismember(log_data.Ensembl_GeneID, met_genes);
data_met = log_data(index_names, :);
logdata = data_met(:, 3:end);
gene_names = data_met{:, 1};

% Local thresholding
up_percentage = 75;
low_percentage = 25;
coreMat = localT2_new(logdata, low_percentage, up_percentage);

% Mapping reactions with gene expression
[numGenes, numSamples] = size(coreMat);
mappedResults = cell(1, numSamples);

for sampleIdx = 1:numSamples
    coreGeneIndices = find(coreMat(:, sampleIdx) == 1);
    coreGeneNames = gene_names(coreGeneIndices);
    expressionValues = table2array(logdata(coreGeneIndices, sampleIdx));
    expressionData = struct('gene', coreGeneNames, 'value', expressionValues);
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model, expressionData);
    mappedResults{sampleIdx} = expressionRxns;
end

% Processing housekeeping genes
genes_table = data_met.Ensembl_GeneID;
index_names = ismember(genes_table, h_k_g.converted_alias);
hkg_met = data_met(index_names, :);
hkg_met_ens = hkg_met(:, "Ensembl_GeneID");
geneIDs = table2cell(hkg_met_ens);

% Obtaining core reactions from core genes
[results] = findRxnsFromGenes(model, geneIDs);

% Obtaining just the housekeeping reaction names
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

% Deleting duplicate names
housekeep_react_unique = unique(housekeep_react);

% Compare the amount of housekeeping core reactions
numSamples = numel(mappedResults);
housekep_core_react = struct('numHousekeepingCoreReactions', [], 'housekeepingCoreReactions', [], 'percentage', []);

% Total number of housekeeping reactions
totalHousekeepingReactions = numel(housekeep_react_unique);

% Comparison and calculation
for i = 1:numSamples
    % Core reactions of the current sample
    coreReactions = mappedResults{i};

    % Find housekeeping reactions in the core reactions
    housekeepingInCore = ismember(housekeep_react_unique, coreReactions);

    % Count and list housekeeping reactions in the core reactions
    housekep_core_react(i).numHousekeepingCoreReactions = sum(housekeepingInCore);
    housekep_core_react(i).housekeepingCoreReactions = housekeep_react_unique(housekeepingInCore);

    % Calculate the percentage
    housekep_core_react(i).percentage = (housekep_core_react(i).numHousekeepingCoreReactions / totalHousekeepingReactions) * 100;
end
