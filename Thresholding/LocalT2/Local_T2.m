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

%% Obtain core-genes
% Core gene matrix
Coregene_Matrix = expression_scoreMatrix >= 1;

% Link with gene names
coreGenesStructure = cell(1, size(Coregene_Matrix, 2));

for i = 1:size(Coregene_Matrix, 2)
    activeGeneIndices = find(Coregene_Matrix(:, i)); % Gene active indexes
    coreGenesStructure{i} = gene_names(activeGeneIndices); % Gene actie names
end
%% Map to Expression
% AND rules model
% Obtén las reglas GPR originales del modelo
originalGPRs = model.grRules;

% Procesa las reglas para eliminar las partes "OR"
processedGPRs = cell(size(originalGPRs));
for i = 1:length(originalGPRs)
    originalRule = originalGPRs{i};
    
    % Elimina cualquier parte "OR" usando expresiones regulares
    processedRule = regexprep(originalRule, '\s*and\s*', ''); % Elimina "or" y espacios cambiar and por or
    
    processedGPRs{i} = processedRule;
end

% Crea un nuevo modelo con las reglas GPR procesadas
newModel = model;
newModel.grRules = processedGPRs;


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
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(newModel, expressionDataSample, 'false');
    Rxns_local25_75 = [Rxns_local25_75, expressionRxns]; % Store the mapped reactions
    geneUsed_local25_75{i} = gene_used; % Store the genes used in the mapping
end

% Define reactionNamesPerSample
reactionNamesPerSample = cell(size(sampleNames));

% Obtén los nombres de las reacciones activas para cada muestra
for i = 1:length(sampleNames)
    activeReactions = find(Rxns_local25_75(:, i) >= 1);  % Encuentra reacciones activas
    reactionNamesPerSample{i} = model.rxns(activeReactions);  % Guarda los nombres de las reacciones activas
end
%% Check core reactions

for i = 1:length(sampleNames)
    ActiveRxns_local25_75{i} = (find(Rxns_local25_75(:, i)>=1));
end
disp(ActiveRxns_local25_75)

%% Obtain the core genes from core reactions
geneData = cell(size(ActiveRxns_local25_75));

% Turn reactions to genes
for i = 1:length(ActiveRxns_local25_75)
    activeReactionsIndices = ActiveRxns_local25_75{i};
    activeReactionsNames = model.rxns(activeReactionsIndices);

    geneList = findGenesFromRxns(newModel, activeReactionsNames);

    geneData{i} = geneList;
end
% Remove the duplicated genes
numColumns = size(geneData, 2);

uniqueGenesArray = cell(1, numColumns);

for i = 1:numColumns
    genesInColumn = geneData{i};

    % If the genes are in nested cells, flatten them out
    if iscell(genesInColumn{1})
        genesInColumn = vertcat(genesInColumn{:});
    end

    % Remove duplicate gene names
    uniqueGenes = unique(genesInColumn);

    uniqueGenesArray{i} = uniqueGenes;
end

%% Compare both and get the names of matching genes
numColumns = length(uniqueGenesArray);
matchedGenesNames = cell(1, numColumns); 

for i = 1:numColumns
    uniqueGenes = uniqueGenesArray{i};    
    activeGenes = coreGenesStructure{i};

    geneMatches = ismember(uniqueGenes, activeGenes);
    matchedGenes = uniqueGenes(geneMatches); 
    matchedGenesNames{i} = matchedGenes;  
end
%% 
%%%%%%%%%%%%%%%%% Compare housekeeping genes and reactions with 
% core genes and reactions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping genes analysis
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

