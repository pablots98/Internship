%%%%%%%%%%%%%%%%%% StanDep Thresholding Method %%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%%% Set up %%%
% Initialize Gurobi and COBRA Toolbox
gurobi_setup;
initCobraToolbox(false);

%%% Load Data %%%
% Load transcriptomics data, housekeeping genes, and metabolic model (Human1 version).
data = readtable('Mod_data.xlsx');  % Transcriptomics data
h_k_g = readtable('housekeeping_ens.csv');  % Housekeeping genes with ensembl IDs
model = readCbModel('Human-GEM_Cobra_v1.01.mat');  % Human1 metabolic model

%%% Metabolic Genes %%%
% Extract genes related to metabolism from the transcriptomics dataset by comparing 
% ensembl IDs with the metabolic model.
model_genes = model.genes;  % Ensembl IDs of genes in the model
index_names = ismember(data.Ensembl_GeneID, model_genes);  % Indices of metabolic genes
data_met = data(index_names, :);  % Rows of dataset matching metabolic genes

%%% StanDep %%%
% Preprocess the metabolic model to ensure it doesn't contain blocked reactions
[~, ~, ~, ~, ~, fluxConsistModel] = findFluxConsistentSubset(model);
model = fluxConsistModel;

% Ensure first column is named 'gene'
data_met.Properties.VariableNames{1} = 'gene';

% Add 1 to avoid negative numbers after log10 transformation
datalog10 = data_met;
datalog10{:, 2:end} = log10(datalog10{:, 2:end} + 1);

%%% Define Bin Limits %%%
% Define bin boundaries for StanDep
maxlog10 = max(max(table2array(datalog10(:, 2:end))));
edgeX = linspace(0, maxlog10, 11);  % 10 bins require 11 boundaries
edgeX = round(edgeX, 1);  % Round to 1 decimal place

%%% Pre-processing for StanDep %%%
% 1. RNA data structure
rnaData = struct();
rnaData.gene = data_met.gene;
rnaData.value = table2array(data_met(:, 2:end));
rnaData.valuebyTissue = table2array(data_met(:, 2:end));
rnaData.Tissue = data.Properties.VariableNames(2:end)';

% 2. Model data structure
modelData = getModelData(rnaData, model);

% 3. Enzyme data structure
spec = getSpecialistEnzymes(model);
prom = getPromEnzymes(model);
enzymeData = comparePromiscuousSpecific(spec, prom, modelData);

%%% Clustering Analysis Parameters %%%
distMethod = 'euclidean';  % Distance method
linkageMethod = 'complete';  % Linkage metric for hierarchical clustering
k = 10;  % Number of clusters

%%% Clustering Analysis %%%
close all;
clustObj = geneExprDist_hierarchy(enzymeData, [], edgeX, k, distMethod, linkageMethod);

%%% Identify Core Reactions and Calculate Ubiquity Score %%%
[coreRxnMat, ~, ~, ~] = models4mClusters1(clustObj, enzymeData.Tissue, model, edgeX, [], [], false, 0, [1 1]);
[ubiScore, ~] = getUbiquityScore_2022(clustObj, edgeX, model);
