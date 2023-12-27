%%%%%%%%%%% Localgini Thresholding Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gurobi_setup;
initCobraToolbox(false);
setenv('GUROBI_PATH', 'C:\gurobi1001\win64\matlab')
changeCobraSolver('gurobi', 'all')
%% Load the data
% Load the transcriptomics data, the housekeeping genes obtained in the 
% previous script, and the metabolic model (Human1 version).
% Transcriptomics data
data = readtable('Mod_data.xlsx');
% Housekeeping genes with the ensembl ids
h_k_g = readtable('NM2ENSG.xlsx');
% Human1 metabolic model.
%model = readCbModel('Human-GEM_Cobra_v1.01.mat')
model = load('SysBio_COBRA_v1.13.mat');
model = model.model

%% Metabolic genes
%Pick the genes related to metabolism
%In the transcriptomics dataset there are a lot of genes, but we are 
% interested just in those ones that take part in the metabolism, so, we 
% are going to compare the ensembl id of both, the transcriptomics dataset
% and the Human1 metabolic model, and pick those ones that are present.
% Create a vector just with the ensembl ids of the genes from the model
model_genes = model.genes;
% Find indexes of the genes that are present in the metabolic model. 
index_names = ismember(data.Ensembl_GeneID, model_genes);
% take the rows of the dataset that match the indexes obtained before, with all the data 
data_met = data(index_names, :);
%Select just the transcriptomics data, and normalize it, with log10, and adding a 1, to avoid having -inf values
% Data just with transcriptomics data, not explanatory columns
data_f = data_met(:, 2:end);
% Normalize the data
logdata = log10(data_f + 1)
logdata

%% Prepare the data
% Separating gene names and data
genes = data_met(:, 1); % First column for gene names
% Converts the table column to a cell array
genes = table2cell(genes);

% Get sample names (row names)
context = logdata.Properties.VariableNames % Assuming that 'rownames' is 
% a function or variable containing the row names

% Calculate the log10 of the fitted expression data
expressionData = data_met(:, 2:end); % Remaining columns for expression data
value = log10(expressionData + 1); % Normalize it

% Creating the structure
geneExpression.value = value;
geneExpression.context = context;
geneExpression.genes = genes;

%% FOLLOW THE STEPS OF THE REPO
% Choice of model extraction method (has to be anyone of ['FASTCORE',
% 'iMAT','MBA','GIMME','INIT','mCADRE'])
MeM = 'FASTCORE';
%Specifying contexts for which models have to be built (has to be in the 
% same format as geneExpression.context) QUE COÑO ES ESTO? MIRARLO MEJOR
contexts = context; % The context or samples that you are interested, in this case, all of them
%Set the upper or lower threshold
ut = 75;
lt = 25;
%Specifying where threshold has to be implied (1==> implied to genes; 2==> 
% implied to enzymes; 3==> implied to reactions)
ThS = 1; % impliying at gene level
%Reactions that are manually given higher importance
biomass_id = find(strcmp(model.rxns,'biomass_reaction'));
atp_demand_id = find(strcmp(model.rxns,'DM_atp_c_'));
coreRxn=[biomass_id atp_demand_id];
%Tolerance level above which reactions are considered as expressed
tol = 1e-4;
%Folder path to save the built models
filename = 'C:\\Users\\PC\\OneDrive\\Documentos\\Systems_Biology_master\\Internship\\Code\\Thresholding\\';
% si no poner una doble y ver que pasa (pista no pasa nada)
%Building context-specific models
%[Models,RxnImp] = buildContextmodels(geneExpression,model,MeM,contexts,ut,lt,ThS,coreRxn,filename,tol)
%save('Models.mat', 'Models');
%save('RxnImp.mat', 'RxnImp');
RxI_check = load('RxnImp.mat')
RxnImp = RxI_check.RxnImp
Model_check = load('Models.mat')
Models = Model_check.Models

%% Housekeeping genes 
ens_hkg = h_k_g.converted_alias;
% Check which ones are metabolic
% Find indexes of the genes that are present in the metabolic model. 
index_names = ismember(ens_hkg, model_genes);
% take the rows of the dataset that match the indexes obtained before, with all the data 
met_hkg = ens_hkg(index_names)

reactions_hkg = findRxnsFromGenes(model, met_hkg)

save('reactions_hkg.mat', 'reactions_hkg')