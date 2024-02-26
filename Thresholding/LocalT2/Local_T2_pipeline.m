%%%%%%%%%%%%%%%%% Local T2 thresholding reactions %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%             Gurobi setup and Cobra Toolbox initialization             %%
gurobi_setup;
initCobraToolbox(false);
setenv('GUROBI_PATH', 'C:\gurobi1001\win64\matlab');
changeCobraSolver('gurobi', 'all');

%%                              Load data                                %%
% Load transcriptomics data, housekeeping genes, and metabolic model
%data = readtable('Merged_data.xlsx'); % Transcriptomics data
data = readtable("data_TPM.xlsx");
h_k_g = readtable('NM2ENSG.xlsx'); 
% h_k_g = readtable('ENS_ID_HKG.xlsx');    % Housekeeping genes with Ensembl IDs
model = load('Human-GEM_Cobra_v1.01.mat'); % Human1 metabolic model
model = model.model;

%%          Ensure the model does not contain blocked reactions          %%
% [fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool, ~, fluxConsistModel] = findFluxConsistentSubset(model);
% model = fluxConsistModel; % load the new model with no blocked reactions
% save('model', 'model');
model = load('model.mat');
model = model.model;
% model_genes = model.genes; % ENSEMBL_IDs of the genes in the model
%%                              Preprocessing                            %%
% Colect the data needed to create the table
Ensembl_id = data(:, 1);
sampleNames = data.Properties.VariableNames(2:end); 
data_to_log = table2array(data(:, 2:end));
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
model_genes = metabolic_genes.Properties.RowNames;
%%                          Local thresholding                           %%
% Set the upper and lower threshold
up_percentage = 75;
low_percentage = 25;

% Use the localT2 function
[adjusted_matrix, expression_scoreMatrix] = localT2_function(metabolic_genes, low_percentage, up_percentage);

%%                          Obtain core-genes                            %% 
gene_names = metabolic_genes.Properties.RowNames;

% Link with gene names
coreGenesStructure = cell(1, size(adjusted_matrix, 2));

for i = 1:size(adjusted_matrix, 2)
    activeGeneIndices = adjusted_matrix(:, i); % Gene active indexes
    coreGenesStructure{i} = gene_names(activeGeneIndices); % Gene active names
end

%%                          Map to Expression                            %%  
% Create expressionData structure to use mapExpressionToReactions function
% expressionData.gene = gene_names; 
% expressionData.value = expression_scoreMatrix; 
% 
% %Initialize variables
% Rxns_local25_75 = [];
% geneUsed_local25_75 = {};
% parsedGPR_local25_75 = {};
% 
% %Iterate over each sample
% for i = 1:width(sampleNames)
%     % Extract the expression data for the current sample
%     expressionDataSample = struct();
%     expressionDataSample.gene = expressionData.gene;
%     expressionDataSample.value = expressionData.value(:, i); % Selecting the column for the current sample
% 
%     % Map expression data to reactions for the current sample
%     [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model, expressionDataSample, 'false');
%     Rxns_local25_75 = [Rxns_local25_75, expressionRxns]; % Store the mapped reactions
%     geneUsed_local25_75{i} = gene_used; % Store the genes used in the mapping
%     parsedGPR_local25_75{i} = parsedGPR; % Store the genes used in the mapping GO THROUGH IT!!!!!
% end
% 
% %Save the results
% save('Rxns_local25_75', "Rxns_local25_75");
% save('geneUsed_local25_75', 'geneUsed_local25_75');
% save('parsedGPR_local25_75', 'parsedGPR_local25_75');

% Load them
Rxns_local25_75 = load('Rxns_local25_75.mat');
Rxns_local25_75 = Rxns_local25_75.Rxns_local25_75;
geneUsed_local25_75 = load('geneUsed_local25_75.mat');
geneUsed_local25_75 = geneUsed_local25_75.geneUsed_local25_75;
parsedGPR_local25_75 = load('parsedGPR_local25_75.mat');
parsedGPR_local25_75 = parsedGPR_local25_75.parsedGPR_local25_75;
% Check there are not repeated reactions

%%                      Set the core reactions                           %%
% Define variable
CoreReact_Sample = cell(size(sampleNames));

% Obtain core reactions names
for i = 1:length(sampleNames)
    activeReactions = find(Rxns_local25_75(:, i) >= 1);  % Active reactions indexes
    CoreReact_Sample{i} = model.rxns(activeReactions);  % Keep active reactions names for each sample
end


%%              Obtain the core genes from core reactions                %% 

numSamples = size(geneUsed_local25_75, 2);

geneNamesBySample = cell(1, numSamples);

for i = 1:numSamples
    numGenes = size(geneUsed_local25_75{1, i}, 2);
    geneNamesList = {};
    for j = 1:numGenes
        currentElement = geneUsed_local25_75{1, i}{1, j};
        
        if iscell(currentElement) && numel(currentElement) == 1 && ischar(currentElement{1}) && ~isempty(strtrim(currentElement{1})) % This part was done by chatGPT but seems to be right
            geneNamesList{end+1} = strtrim(currentElement{1});
        elseif ischar(currentElement) && ~isempty(strtrim(currentElement))
            geneNamesList{end+1} = strtrim(currentElement);
        end
    end

    geneNamesBySample{i} = unique(geneNamesList);
end

%% Compare the core genes obtained with the previous core genes

numSamples = length(geneNamesBySample);
difGenesSCheck = cell(1, numSamples);

for i = 1:numSamples
    genesSample = geneNamesBySample{i};
    coreGenes = coreGenesStructure{i};

    genesPresent = ismember(genesSample, coreGenes);

    if all(genesPresent)
        fprintf('All genes from sample %d are presente in coreGenesStructure.\n', i);
    else
        fprintf('Some genes from the sample %d are not present in coreGenesStructure.\n', i);
        missingGenes = genesSample(~genesPresent);
        difGenesSCheck{i} = missingGenes;
    end
end

% The results are horrible / Maybe woks better establishing the core-genes,
% and from them going to core-reactions.
%%                                                                       %%
%%%%%%%%%%%%%%%%% Compare housekeeping genes and reactions with %%%%%%%%%%%
%%%%%%%%%%%%%%%%%             core genes and reactions          %%%%%%%%%%%

%%                    Housekeeping genes analysis                        %%
ens_hkg = h_k_g.converted_alias;
index_names = ismember(ens_hkg, model_genes);
met_hkg = ens_hkg(index_names);

% Use this function to find the housekeeping reactions names
HKG_react = findRxnsFromGenes(model, met_hkg); % Gives a structure, so the for loop to extract the names

% Extract unique housekeeping reaction names
fields = fieldnames(HKG_react);
housekeep_react = {};
for i = 1:length(fields)
    field = fields{i};
    cellArray = HKG_react.(field);
    for j = 1:size(cellArray, 1)
        firstcol = cellArray{j, 1};
        housekeep_react{end+1, 1} = firstcol;
    end
end

% Just want the name one time
housekeep_react_unique = unique(housekeep_react);

%%           Compare the number of housekeeping core genes               %%
housekep_core_gene_L = struct("numHousekeepingCoreGenes", [], 'housekeepingCoreGenes', [], 'percentage', []);
totalHousekeepingGenes = numel(met_hkg);

for i = 1:width(sampleNames)
    coreGenes = coreGenesStructure{i};
    housekeepingGeneNames = met_hkg;
    housekeepingInCore_g = ismember(housekeepingGeneNames, coreGenes);
    housekep_core_gene_L(i).numHousekeepingCoreGenes = sum(housekeepingInCore_g);
    housekep_core_gene_L(i).housekeepingCoreGenes = housekeepingGeneNames(housekeepingInCore_g);
    housekep_core_gene_L(i).percentage = (housekep_core_gene_L(i).numHousekeepingCoreGenes / totalHousekeepingGenes) * 100;
end

HK_G_acc_LT = [housekep_core_gene_L.percentage];
HK_G_acc_LT_col = HK_G_acc_LT(:);
T_Genes = table(HK_G_acc_LT_col, 'VariableNames', {'HK_G_acc_LT'});
writetable(T_Genes, "HK_G_acc_LT.xlsx", 'WriteRowNames', false);

%%          Compare the number of housekeeping core reactions            %%
housekep_core_react_L = struct('numHousekeepingCoreReactions', [], 'housekeepingCoreReactions', [], 'percentage', []);
totalHousekeepingReactions = numel(housekeep_react_unique);

for i = 1:width(sampleNames)
    coreReactions = CoreReact_Sample{i};
    housekeepingInCore = ismember(housekeep_react_unique, coreReactions);
    housekep_core_react_L(i).numHousekeepingCoreReactions = sum(housekeepingInCore);
    housekep_core_react_L(i).housekeepingCoreReactions = housekeep_react_unique(housekeepingInCore);
    housekep_core_react_L(i).percentage = (housekep_core_react_L(i).numHousekeepingCoreReactions / totalHousekeepingReactions) * 100;
end

HK_R_acc_LT = [housekep_core_react_L.percentage];
HK_R_acc_LT_col = HK_R_acc_LT(:);
T_React = table(HK_R_acc_LT_col, 'VariableNames', {'HK_R_acc_LT'});
writetable(T_React, "HK_R_acc_LT.xlsx", 'WriteRowNames', false);

%% Check the samples with the same accuracy with genes
% find the indexes of the same values
 equalIndex = {};

 for i = height(HK_G_acc_LT_col)
     index = find(HK_G_acc_LT_col == HK_G_acc_LT_col(i));

     if length(index) > 1 && ~any(cellfun(@(x) isequal(x, index), equalIndex))
         equalIndex{end+1} = index
     end
 end

% Check the HK core genes names 
genesArray = housekep_core_gene_L.housekeepingCoreGenes;

for i = 1:length(equalIndex)
    indicesGrupo = equalIndex{i};
    genesGrupo = genesArray{indicesGrupo(1)};

    sonIguales = true;

    for j = 2:length(indicesGrupo)
        genesComparar = genesArray{indicesGrupo(j)};

        if length(genesGrupo) ~= length(genesComparar) || ~all(strcmp(genesGrupo, genesComparar))
            sonIguales = false;
            break;
        end
    end

    resultadosComparacion{i} = sonIguales;
end

%% Check the samples with the same accuracy with reactions
% find the indexes of the same values
 equalIndex = {};

 for i = height(HK_R_acc_LT_col)
     index = find(HK_R_acc_LT_col == HK_R_acc_LT_col(i));

     if length(index) > 1 && ~any(cellfun(@(x) isequal(x, index), equalIndex))
         equalIndex{end+1} = index
     end
 end

% Check the HK core genes names 
genesArray = housekep_core_react_L.housekeepingCoreReactions;

for i = 1:length(equalIndex)
    indicesGrupo = equalIndex{i};
    genesGrupo = genesArray{indicesGrupo(1)};

    sonIguales = true;

    for j = 2:length(indicesGrupo)
        genesComparar = genesArray{indicesGrupo(j)};

        if length(genesGrupo) ~= length(genesComparar) || ~all(strcmp(genesGrupo, genesComparar))
            sonIguales = false;
            break;
        end
    end

    resultadosComparacion{i} = sonIguales;
end