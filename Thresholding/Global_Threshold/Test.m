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

results_CoreGene = results;
results_CoreGene(1, :) = [];
disp(results_CoreGene);

%% Test
%%      Global thresholding to determine core and non-core genes         %%
% Inicializar results_CoreGenes con nombres de muestra como nombres de columnas

results_CoreGene = array2table(cell(1, length(sampleNames)), 'VariableNames', sampleNames);

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
    
    % Asignar los IDs de los genes centrales a la celda correspondiente
    results_CoreGene{1, col_name} = {coreGeneIDs};
end

% Mostrar los resultados
disp(results_CoreGene);

%% l
% Create expressionData structure to use mapExpressionToReactions function
expressionData.gene = gene_names; 
expressionData.value = expression_scoreMatrix; 

%Initialize variables
Rxns_local25_75 = [];
geneUsed_local25_75 = {};
parsedGPR_local25_75 = {};

%Iterate over each sample
for i = 1:width(sampleNames)
    % Extract the expression data for the current sample
    expressionDataSample = struct();
    expressionDataSample.gene = expressionData.gene;
    expressionDataSample.value = expressionData.value(:, i); % Selecting the column for the current sample

    % Map expression data to reactions for the current sample
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model, expressionDataSample, 'false');
    Rxns_local25_75 = [Rxns_local25_75, expressionRxns]; % Store the mapped reactions
    geneUsed_local25_75{i} = gene_used; % Store the genes used in the mapping
    parsedGPR_local25_75{i} = parsedGPR; % Store the genes used in the mapping GO THROUGH IT!!!!!
end

%Save the results
save('Rxns_local25_75', "Rxns_local25_75");
save('geneUsed_local25_75', 'geneUsed_local25_75');
save('parsedGPR_local25_75', 'parsedGPR_local25_75');

% Load them
Rxns_local25_75 = load('Rxns_local25_75.mat');
Rxns_local25_75 = Rxns_local25_75.Rxns_local25_75;
geneUsed_local25_75 = load('geneUsed_local25_75.mat');
geneUsed_local25_75 = geneUsed_local25_75.geneUsed_local25_75;
parsedGPR_local25_75 = load('parsedGPR_local25_75.mat');
parsedGPR_local25_75 = parsedGPR_local25_75.parsedGPR_local25_75;
