%%%%%%%%%%%%%%%%%% StanDep thresholding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%%% Set up 
gurobi_setup;
initCobraToolbox(false);


%%% Load the data
% Load the transcriptomics data, the housekeeping genes obtained in the 
% previous script, and the metabolic model (Human1 version).
% Transcriptomics data
data = readtable('Mod_data.xlsx');
% Housekeeping genes with the ensembl ids
h_k_g = readtable('housekeeping_ens.csv');
% Human1 metabolic model.
model = readCbModel('Human-GEM_Cobra_v1.01.mat');

%%% Metabolic genes
% Pick the genes related to metabolism
% In the transcriptomics dataset there are a lot of genes, but we are 
% interested just in those ones that take part in the metabolism, so, we 
% are going to compare the ensembl id of both, the transcriptomics dataset 
% and the Human1 metabolic model, and pick those ones that are present.
% Create a vector just with the ensembl ids of the genes from the model
model_genes = model.genes
% Find indexes of the genes that are present in the metabolic model. 
index_names = ismember(data.Ensembl_GeneID, model_genes);
% take the rows of the dataset that match the indexes obtained before, with all the data 
data_met = data(index_names, :)


                                % StanDep 
% First, in all model extraction it is good to pre_process the model to
% ensure it does not contain blocked reactions

[fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool, ~, fluxConsistModel] = findFluxConsistentSubset(model);
model = fluxConsistModel

%Ensuring that the first colum is names 'gene' 
data_met.Properties.VariableNames{1} = 'gene';
%Although this is not documented clearly, standep will perform a log10 on
%the expression data. To ensure that this does not introduce negative
%numbers we add 1 everywhere.
datalog10 = data_met
datalog10{:,2:end} = log10(datalog10{:,2:end}+1)

%%% Define the limits of the bins:
%Standep requires us to define the bin boundaries. IMPORTANTLY these
%boundaries have to be in a log scale...


maxlog10 = max(max(table2array(datalog10(:,2:end))));
edgeX = linspace(0,maxlog10,11); %If we want 10 bins we need 11 
edgeX = round(edgeX,1); %These are our bin bounds!

%%% Pre-processing of the data for Standep
% Standep requires a very specific pre-processing for the inputs!
% 1st, the RNA data structure
rnaData = struct();
rnaData.gene = data_met.gene;
rnaData.value = table2array(data_met(:,2:end));
rnaData.valuebyTissue = table2array(data_met(:,2:end));
rnaData.Tissue = data.Properties.VariableNames(2:end)';
% Pre-processing of the metabolic model
% 2nd the model data structure
modelData = getModelData(rnaData,model);
%%% Pre-processing of the enzyme data 
% 3rd the enzyme data structure
spec = getSpecialistEnzymes(model);  
prom = getPromEnzymes(model);
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);
%%% Set the parameters for the clustering analysis
% We then set our parameters...
distMethod = 'euclidean'; % distance method  maybe think one more robust than euclidean.
linkageMethod = 'complete'; % linkage metric for hierarchical clustering
k = 10;
%%% Clustering analysis
close all % This is important to ensure that all standep plots have the right number to be saved if we wish
clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,k,distMethod,linkageMethod);
%%% Identify core reactions and calculate ubiquity score
[coreRxnMat,enzTis,cutOff,thr] = models4mClusters1(clustObj,enzymeData.Tissue,model,edgeX,[],[],false,0,[1 1]); %CoreRxnMat defines the core set of reactions in a general sense
[ubiScore,uScore] = getUbiquityScore_2022(clustObj,edgeX,model); %ubiScore is the ubiquity score per rxns to use in mCadre!