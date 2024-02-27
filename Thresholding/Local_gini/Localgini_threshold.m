%%%%%%%%%%%%%%%%%%% Localgini Thresholding Method %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%                          Set up                                       %%
gurobi_setup;
initCobraToolbox(false);
setenv('GUROBI_PATH', 'C:\gurobi1001\win64\matlab')
changeCobraSolver('gurobi', 'all')
%%                      Load the data                                    %%
% Transcriptomics data
data = readtable('Merged_data.xlsx');
% Housekeeping genes with the ensembl ids
%h_k_g = readtable('NM2ENSG.xlsx');
h_k_g = readtable('ENS_ID_HKG.xlsx');
% Human1 metabolic model.
% model = load('Human-GEM_Cobra_v1.01.mat')
%model = load('SysBio_COBRA_v1.13.mat');
model = load('SysBio_COBRA_v1.17_consensus.mat'); % Human1 metabolic model
model = model.myModel;
%model = model.model

%%         Ensure the model does not contain blocked reactions           %%
% [fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool, ~, fluxConsistModel] = findFluxConsistentSubset(model);
% model = fluxConsistModel; % load the new model with no blocked reactions
% save('model', 'model');
model = load('model.mat');
model = model.model;
model_genes = model.genes; % ENSEMBL_IDs of the genes in the model

%%                           Preporcessing                               %%
%Colect the data needed to create the table
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

%%                         Prepare the data                              %%
% Separating gene names and data
genes = metabolic_genes.Properties.RowNames; % First column for gene names

% Get sample names (row names)
context = metabolic_genes.Properties.VariableNames % Assuming that 'rownames' is a function or variable containing the row names

% Calculate the log10 of the fitted expression data
value = metabolic_genes; % Normalize it

% Creating the structure
geneExpression = struct;
geneExpression.value = value;
geneExpression.context = context;
geneExpression.genes = genes;

%%                          Repo pipeline                                %%                           
% Choice of model extraction method (has to be anyone of ['FASTCORE','iMAT','MBA','GIMME','INIT','mCADRE'])
MeM = 'FASTCORE';

% The context or samples that you are interested, in this case, all of them
contexts = context;

% Set upper and lower thresholds
ut = 75;
lt = 25;

% Specifying where threshold has to be implied (1==> implied to genes; 2==> implied to enzymes; 3==> implied to reactions)
ThS = 1; % impliying at gene level

% Reactions that are manually given higher importance
biomass_id = find(strcmp(model.rxns,'biomass_reaction'));
atp_demand_id = find(strcmp(model.rxns,'DM_atp_c_'));
coreRxn=[biomass_id atp_demand_id];

% Tolerance level above which reactions are considered as expressed
tol = 1e-4;

% Folder path to save the built models Here put your path to save the models, 
% this way (\\) is for windows, for Linux, and I think for mac it should be (/)
filename = 'C:\\Users\\PC\\OneDrive\\Documentos\\Systems_Biology_master\\Internship\\Internship\\Thresholding\\Local_gini\\';

%Building context-specific models
% [Models,RxnImp] = buildContextmodels(geneExpression,model,MeM,contexts,ut,lt,ThS,coreRxn,filename,tol);
% save('Models.mat', 'Models');
% save('RxnImp.mat', 'RxnImp');
RxI_check = load('RxnImp.mat');
RxnImp = RxI_check.RxnImp;
Model_check = load('Models.mat');
Models = Model_check.Models;

%%                          Housekeeping genes                           %%
ens_hkg = h_k_g.ENSEMBL_ID; 
index_names = ismember(ens_hkg, genes);
met_hkg = ens_hkg(index_names);

% Find reactions from the housekeeping genes
reactions_hkg = findRxnsFromGenes(model, met_hkg);

%Create an array with all the housekeeping reactions
fields = fieldnames(reactions_hkg); 
housekeep_react = {}; 

for i = 1:length(fields)
    field = fields{i}; 
    cellArray = reactions_hkg.(field); 

    for j = 1:size(cellArray, 1)
        firstcol = cellArray{j, 1}; 
        housekeep_react{end+1, 1} = firstcol;
    end
end
housekeep_react;

% Finding unique elements and their frequency
[uniques, ~, index] = unique(housekeep_react);

%%                  Check the housekeeping core genes                    %%

num_models = length(Models);
results = zeros(num_models, 1); 

for i = 1:num_models
    current_model = Models{i};
    
    if isfield(current_model, 'genes')
        genes_model = current_model.genes;
        
        num_coincidences = sum(ismember(genes_model, met_hkg));
        
        percentage_coincidences = (num_coincidences / length(met_hkg)) * 100;
        
        results(i) = percentage_coincidences;
    else
        results(i) = NaN;
    end
end
HK_G_accuracy = array2table(results, 'VariableNames', {'PercentageCoincidences'});

HK_G_acc_LG = [HK_G_accuracy.PercentageCoincidences];
HK_G_acc_LG_col = HK_G_acc_LG(:);
T_Gene = table(HK_G_acc_LG_col, 'VariableNames', {'HK_G_acc_LG'});
writetable(T_Gene, "HK_G_acc_LG.xlsx", 'WriteRowNames', false);

%%              Check the housekeeping core reactions                    %%

num_models = length(Models);
results = zeros(num_models, 1); 
for i = 1:num_models
    current_model = Models{i};
    
    if isfield(current_model, 'rxns')
        reactions_model = current_model.rxns;
        
        num_coincidences = sum(ismember(reactions_model, uniques));
       
        percentage_coincidences = (num_coincidences / length(uniques)) * 100;
        
        results(i) = percentage_coincidences;
    else
        results(i) = NaN;
    end
end
HK_R_accuracy = array2table(results, 'VariableNames', {'PercentageCoincidencesReactions'});

HK_R_acc_LG = [HK_R_accuracy.PercentageCoincidencesReactions];
HK_R_acc_LG_col = HK_R_acc_LG(:);
T_React = table(HK_R_acc_LG_col, 'VariableNames', {'HK_R_acc_LG'});
writetable(T_React, "HK_R_acc_LG.xlsx", 'WriteRowNames', false);