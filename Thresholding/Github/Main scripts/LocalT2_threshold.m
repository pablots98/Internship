%%%%%%%%%%%%%%%%%%%% LocalT2 Thresholding Method %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;

%%% Load Data %%%
% Load transcriptomics data, housekeeping genes list, and Human1 metabolic model
data = readtable('Mod_data.xlsx');  % Transcriptomics data
h_k_g = readtable('housekeeping_ens.csv');  % Housekeeping genes with Ensembl IDs
model = readCbModel('Human-GEM_Cobra_v1.01.mat');  % Human1 metabolic model

%%% Metabolic Genes %%%
% Extract genes related to metabolism from transcriptomics dataset
model_genes = model.genes;  % Extract Ensembl IDs of genes from the metabolic model
index_names = ismember(data.Ensembl_GeneID, model_genes);  % Find indices of metabolic genes in transcriptomics data
data_met = data(index_names, :);  % Extract metabolic genes from transcriptomics data

%%% Normalize Transcriptomics Data %%%
% Select transcriptomics data and normalize using log10 transformation
data_f = data_met(:, 2:end);  % Extract transcriptomics data
logdata = log10(data_f + 1);  % Normalize data with log10 transformation

%%% LocalT2 Core Analysis %%%
% Set lower and upper percentiles for threshold calculation in LocalT2 method
low_percentage = 25;
up_percentage = 75;
% Calculate core matrix using LocalT2 method
coreMat = localT2_core(logdata, low_percentage, up_percentage);

%%% Core Metabolic Genes %%%
% Identify core metabolic genes based on LocalT2 results
core_indices = find(mean(coreMat, 2) > 0.5);
gene_names = data_met{:, 1};
core_genes = gene_names(core_indices);
core_genes_tab = table(core_genes);

%%% Housekeeping Genes Analysis %%%
% Identify metabolic housekeeping genes
genes_table = data_met.Ensembl_GeneID;
index_names = ismember(genes_table, h_k_g.converted_alias);
hkg_met = data_met(index_names, :);
hkg_met_ens = hkg_met(:, "Ensembl_GeneID");
% Identify core housekeeping genes
correctly_identified = intersect(core_genes_tab.core_genes, hkg_met_ens.Ensembl_GeneID);
hkg_met_ens = table2array(hkg_met(:, 'Ensembl_GeneID'));
percentage_correct = (length(correctly_identified) / length(hkg_met_ens)) * 100;
disp(['Correctly identified housekeeping genes: ', num2str(percentage_correct), '%']);

%%% Metabolic Housekeeping Genes %%%
% Analyze coverage of metabolic housekeeping genes in each sample
sel_rows = coreMat(index_names == 1, :);
num_ones = sum(sel_rows);
num_rows = size(sel_rows, 1);
prop_genes = num_ones / num_rows;
expression_col = data_met.Properties.VariableNames(2:end);
results_col_FPKM = array2table(prop_genes, 'VariableNames', expression_col)

%------------------------FPKM to TPM Conversion---------------------------%
% Load data again for TPM conversion
data = readtable('Mod_data.xlsx');
h_k_g = readtable('housekeeping_ens.csv');
model = readCbModel('Human-GEM_Cobra_v1.01.mat');

% Convert FPKM values to TPM
data_matrix = data{:, 2:end};
column_sums = sum(data_matrix, 1);
normalized_matrix = (data_matrix ./ column_sums) * 1e6;
normalized_table = array2table(normalized_matrix, 'VariableNames', data.Properties.VariableNames(2:end));
data(:, 2:end) = normalized_table;

%%% Metabolic Genes (TPM) %%%
% Repeat metabolic gene extraction for TPM data
model_genes = model.genes;
index_names = ismember(data.Ensembl_GeneID, model_genes);
data_met = data(index_names, :);

%%% Normalize Transcriptomics Data (TPM) %%%
data_f = data_met(:, 2:end);
logdata = log10(data_f + 1);

%%% LocalT2 Core Analysis (TPM) %%%
coreMat = localT2_core(logdata, low_percentage, up_percentage);

%%% Core Metabolic Genes (TPM) %%%
core_indices = find(mean(coreMat, 2) > 0.5);
gene_names = data_met{:, 1};
core_genes = gene_names(core_indices);
core_genes_tab = table(core_genes);

%%% Housekeeping Genes Analysis (TPM) %%%
genes_table = data_met.Ensembl_GeneID;
index_names = ismember(genes_table, h_k_g.converted_alias);
hkg_met = data_met(index_names, :);
hkg_met_ens = hkg_met(:, "Ensembl_GeneID");
correctly_identified = intersect(core_genes_tab.core_genes, hkg_met_ens.Ensembl_GeneID);
hkg_met_ens = table2array(hkg_met(:, 'Ensembl_GeneID'));
percentage_correct = (length(correctly_identified) / length(hkg_met_ens)) * 100;
disp(['Correctly identified housekeeping genes (TPM): ', num2str(percentage_correct), '%']);

%%% Metabolic Housekeeping Genes (TPM) %%%
sel_rows = coreMat(index_names == 1, :);
num_ones = sum(sel_rows);
num_rows = size(sel_rows, 1);
prop_genes = num_ones / num_rows;
results_col_TPM = array2table(prop_genes, 'VariableNames', expression_col)
