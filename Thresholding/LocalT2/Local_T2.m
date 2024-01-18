%%%%%%%%%%%%%%%%% Local T2 thresholding reactions %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%             Gurobi setup and Cobra Toolbox initialization             %%
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
model_genes = model.genes; % ENSEMBL_IDs of the genes in the model

%%          Ensure the model does not contain blocked reactions          %%
[fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool, ~, fluxConsistModel] = findFluxConsistentSubset(model);
model = fluxConsistModel; % load the new model with no blocked reactions

%%                              Preprocessing                            %%
% Colect the data needed to create the table
genes = data(:, 2); % take the Entrez_ID (needed for findUsedGenesLevels)
Ensembl_id = data(:, 1);
sampleNames = data.Properties.VariableNames(7:end); 
data_to_log = table2array(data(:, 7:end));
logged_data = log10(data_to_log + 1); % Normalize the data (+1 to avoid having 0 values) 

%Create the new table with the data obtained before
log_data = [Ensembl_id, genes, array2table(logged_data)]; % Load the data on the new table  
log_data.Properties.VariableNames(3:end) = sampleNames; % Variables names 
log_data.Properties.VariableNames{2} = 'gene'; % Change Entrez_ID name to gene (To run findUsedGenesLevels)

%%              Processing the gene expression dataset                   %% REVIEW IT
% Define new cells and new Matrixes
allGenes = {};
geneExpressionMatrix = [];

for i = 1: width(sampleNames)
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

%%          Taking gene expression from the metabolic genes              %%
% Obtained the metabolic genes, also take their normalized expression
met_genes_names = metabolic_genes.Properties.RowNames; % 2495
index_names = ismember(log_data.Ensembl_GeneID, met_genes_names); % here 2471 genes
data_met = log_data(index_names, :);
logdata = data_met(:, 3:end);

% Create a vector with metabolic genes ENSEMBL_IDs
gene_names = data_met{:, 1};

%%                          Local thresholding                           %%
% Set the upper and lower threshold
up_percentage = 75;
low_percentage = 25;

% Use the localT2 function
[expression_scoreMatrix, coreMat] = localT2_new(logdata, low_percentage, up_percentage);

%%                          Obtain core-genes                            %%
% Core gene matrix
Coregene_Matrix = expression_scoreMatrix >= 1;

% Link with gene names
coreGenesStructure = cell(1, size(Coregene_Matrix, 2));

for i = 1:size(Coregene_Matrix, 2)
    activeGeneIndices = find(Coregene_Matrix(:, i)); % Gene active indexes
    coreGenesStructure{i} = gene_names(activeGeneIndices); % Gene active names
end


%% AND rules CHECK IT SERIOUSLY
% % Get the original GPR rules of the model
% originalGPRs = model.grRules;
% 
% % Process the rules to remove the "AND" parts.
% processedGPRs = cell(size(originalGPRs));
% for i = 1:length(originalGPRs)
%     originalRule = originalGPRs{i};
% 
%     % Remove any "OR" parts using regular expressions
%     processedRule = regexprep(originalRule, '\s*or\s*', ''); % Remove "or" and spaces change and to "or".
% 
%     processedGPRs{i} = processedRule;
% end
% 
% % Creates a new model with the processed GPR rules
% newModel = model;
% newModel.grRules = processedGPRs;

%%                          Map to Expression                            %%
% Create expressionData structure to use mapExpressionToReactions function
expressionData.gene = gene_names; 
expressionData.value = expression_scoreMatrix; 

% Initialize variables
Rxns_local25_75 = [];
geneUsed_local25_75 = {};
parsedGPR_local25_75 = {};

% Iterate over each sample
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
%%

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
% creo que solo me da los números, ssacar los nombres

%% Obtain the core genes from core reactions
% Supongamos que 'data' es tu cell array original
numColumns = size(parsedGPR_local25_75, 2);
result = cell(1, numColumns); % Inicializa la cell array de resultados

for i = 1:numColumns
    geneNames = {}; % Para almacenar temporalmente los nombres de los genes
    for j = 1:size(parsedGPR_local25_75{i}, 1)
        % Revisa si el contenido de la cell es otra cell array
        if iscell(parsedGPR_local25_75{i}{j})
            % Si es otra cell array, recórrela
            for k = 1:size(parsedGPR_local25_75{i}{j}, 1)
                % Continúa con este proceso si te encuentras con más cell arrays anidadas
                if iscell(parsedGPR_local25_75{i}{j}{k})
                    for l = 1:size(parsedGPR_local25_75{i}{j}{k}, 1)
                        % Aquí finalmente asumimos que has llegado a los nombres de los genes
                        geneNames = [geneNames; parsedGPR_local25_75{i}{j}{k}{l}];
                    end
                else
                    geneNames = [geneNames; parsedGPR_local25_75{i}{j}{k}];
                end
            end
        else
            geneNames = [geneNames; parsedGPR_local25_75{i}{j}];
        end
    end
    % Almacena los nombres únicos de los genes en la cell array de resultados
    result{i} = unique(geneNames);
end

%% Compare both and get the names of matching genes
numColumns = length(result);
matchedGenesNames = cell(1, numColumns); 

for i = 1:numColumns
    uniqueGenes = result{i};    
    activeGenes = coreGenesStructure{i};

    geneMatches = ismember(uniqueGenes, activeGenes);
    matchedGenes = uniqueGenes(geneMatches); 
    matchedGenesNames{i} = matchedGenes;  
end
%% 
%%%%%%%%%%%%%%%%% Compare housekeeping genes and reactions with 
% core genes and reactions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping genes analysis
ens_hkg = h_k_g.converted_alias;
% Check which ones are metabolic
% Find indexes of the genes that are present in the metabolic model. 
index_names = ismember(ens_hkg, model_genes);
% take the rows of the dataset that match the indexes obtained before, with all the data 
met_hkg = ens_hkg(index_names);

results = findRxnsFromGenes(model, met_hkg);


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
numSample_genes = numel(coreGenesStructure);
housekep_core_gene = struct("numHousekeepingCoreGenes", [], 'housekeepingCoreGenes', [], 'percentage', []);
totalHousekeepingGenes = numel(met_hkg);

for i = 1:numSample_genes
    coreGenes = coreGenesStructure{i};
    housekeepingGeneNames = met_hkg;
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

