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

%%          Ensure the model does not contain blocked reactions          %%
% [fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool, ~, fluxConsistModel] = findFluxConsistentSubset(model);
% model = fluxConsistModel; % load the new model with no blocked reactions
% save('model', 'model');
model = load('model.mat');
model = model.model;
model_genes = model.genes; % ENSEMBL_IDs of the genes in the model
%%                              Preprocessing                            %%
% Colect the data needed to create the table
% genes = data(:, 2); % take the Entrez_ID (needed for findUsedGenesLevels) WE CAN DELETE IT
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

% Delete the NaN values
metabolic_genes = rmmissing(metabolic_genes, 'MinNumMissing', size(metabolic_genes, 2));

%%                          Local thresholding                           %%
% Set the upper and lower threshold
up_percentage = 75;
low_percentage = 25;

% Use the localT2 function
[expression_scoreMatrix, coreMat] = localT2_new(metabolic_genes, low_percentage, up_percentage);

%%                          Obtain core-genes                            %%
% Core gene matrix
Coregene_Matrix = expression_scoreMatrix >= 1;
gene_names = metabolic_genes.Properties.RowNames;

% Link with gene names
coreGenesStructure = cell(1, size(Coregene_Matrix, 2));

for i = 1:size(Coregene_Matrix, 2)
    activeGeneIndices = find(Coregene_Matrix(:, i)); % Gene active indexes
    coreGenesStructure{i} = gene_names(activeGeneIndices); % Gene active names
end


%%                          Map to Expression                            %% WHY i CAN'T SEE GENE_USED 
% % Create expressionData structure to use mapExpressionToReactions function
% expressionData.gene = gene_names; 
% expressionData.value = expression_scoreMatrix; 
% 
% % Initialize variables
% Rxns_local25_75 = [];
% geneUsed_local25_75 = {};
% parsedGPR_local25_75 = {};
% 
% 
% % I COMMENTED BECAUSE IT TAKES A LOT OF TIME, SO I SAVE THE RESULTS AND
% % JUST LOAD THE,
% 
% % Iterate over each sample
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
% % Save the results
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

%%                      Set the core reactions                           %%
% Define variable
reactionNamesPerSample = cell(size(sampleNames));

% Obtain core reactions names
for i = 1:length(sampleNames)
    activeReactions = find(Rxns_local25_75(:, i) >= 1);  % Active reactions indexes
    reactionNamesPerSample{i} = model.rxns(activeReactions);  % Keep active reactions names for each sample
end

% %%              Obtain the core genes from core reactions                %% THERE IS SOMETHING WRONG HERE, I HAVE MORE GENES THAT I WOULD BE EXPECTED
% % Initialize variables
% numColumns = size(parsedGPR_local25_75, 2);
% CoreG_from_CoreR = cell(1, numColumns); 
% 
% for i = 1:numColumns
%     geneNames = {}; % Keep temporarly the gene's names
%     for j = 1:size(parsedGPR_local25_75{i}, 1)
%         % Check if the data is inside other array
%         if iscell(parsedGPR_local25_75{i}{j})
%             % If there is another cell, it go through it
%             for k = 1:size(parsedGPR_local25_75{i}{j}, 1)
%                 % Continue trying to find new nested cells
%                 if iscell(parsedGPR_local25_75{i}{j}{k})
%                     for l = 1:size(parsedGPR_local25_75{i}{j}{k}, 1)
%                         % Arrived to gene names
%                         geneNames = [geneNames; parsedGPR_local25_75{i}{j}{k}{l}];
%                     end
%                 else
%                     geneNames = [geneNames; parsedGPR_local25_75{i}{j}{k}];
%                 end
%             end
%         else
%             geneNames = [geneNames; parsedGPR_local25_75{i}{j}];
%         end
%     end
%     % Keep just the unique names
%     CoreG_from_CoreR{i} = unique(geneNames);
% end
% 
% %%          Compare both and get the names of matching genes             %% THIS ONE WAS A SANITY CHECK TO SEE IF ALL OF THE GENES WERE INCLUDED, APPARENTLY NO
% numColumns = length(CoreG_from_CoreR);
% matchedGenesNames = cell(1, numColumns); 
% 
% for i = 1:numColumns
%     uniqueGenes = CoreG_from_CoreR{i}; 
%     activeGenes = coreGenesStructure{i};
% 
%     % Check if uniqueGenes and activeGenes are cell arrays of strings
%     if isnumeric(uniqueGenes)
%         uniqueGenes = cellstr(num2str(uniqueGenes));
%     end
% 
%     if isnumeric(activeGenes)
%         activeGenes = cellstr(num2str(activeGenes));
%     end
% 
%     geneMatches = ismember(uniqueGenes, activeGenes);
%     matchedGenes = uniqueGenes(geneMatches); 
%     matchedGenesNames{i} = matchedGenes;  
% end

%%                                                                       %%
%%%%%%%%%%%%%%%%% Compare housekeeping genes and reactions with %%%%%%%%%%%
%%%%%%%%%%%%%%%%%             core genes and reactions          %%%%%%%%%%%

%%                    Housekeeping genes analysis                        %%
ens_hkg = h_k_g.converted_alias;
index_names = ismember(ens_hkg, model_genes);
met_hkg = ens_hkg(index_names);

% Use this function to find the housekeeping reactions names
results = findRxnsFromGenes(model, met_hkg); % Gives a structure, so for loop to extract the names


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

% Just want the name one time
housekeep_react_unique = unique(housekeep_react);
disp(housekeep_react_unique);

%%           Compare the number of housekeeping core genes               %%
housekep_core_gene_L = struct("numHousekeepingCoreGenes", [], 'housekeepingCoreGenes', [], 'percentage', []);
totalHousekeepingGenes = numel(met_hkg);

for i = 1:numColumns
    coreGenes = coreGenesStructure{i};
    housekeepingGeneNames = met_hkg;
    if isnumeric(coreGenes)
        coreGenes = cellstr(num2str(coreGenes));
    end
    housekeepingInCore_g = ismember(housekeepingGeneNames, coreGenes);
    housekep_core_gene_L(i).numHousekeepingCoreGenes = sum(housekeepingInCore_g);
    housekep_core_gene_L(i).housekeepingCoreGenes = housekeepingGeneNames(housekeepingInCore_g);
    housekep_core_gene_L(i).percentage = (housekep_core_gene_L(i).numHousekeepingCoreGenes / totalHousekeepingGenes) * 100;
end

%%          Compare the number of housekeeping core reactions            %%
housekep_core_react_L = struct('numHousekeepingCoreReactions', [], 'housekeepingCoreReactions', [], 'percentage', []);
totalHousekeepingReactions = numel(housekeep_react_unique);

for i = 1:numColumns
    coreReactions = reactionNamesPerSample{i};
    housekeepingInCore = ismember(housekeep_react_unique, coreReactions);
    housekep_core_react_L(i).numHousekeepingCoreReactions = sum(housekeepingInCore);
    housekep_core_react_L(i).housekeepingCoreReactions = housekeep_react_unique(housekeepingInCore);
    housekep_core_react_L(i).percentage = (housekep_core_react_L(i).numHousekeepingCoreReactions / totalHousekeepingReactions) * 100;
end

