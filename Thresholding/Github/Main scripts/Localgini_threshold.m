%%%%%%%%%%%%%%%%%%% Localgini Thresholding Method %%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%%% Load the data %%%

% Load transcriptomics data
data = readtable('Mod_data.xlsx');

% Load housekeeping genes with the Ensembl IDs
h_k_g = readtable('housekeeping_ens.csv');

% Load Human1 metabolic model
model = readCbModel('Human-GEM_Cobra_v1.01.mat');

%%% Metabolic genes %%%

% Extract genes related to metabolism from the transcriptomics dataset
model_genes = model.genes; % Ensembl IDs of genes from the metabolic model
index_names = ismember(data.Ensembl_GeneID, model_genes); % Find indices of metabolic genes in the dataset
data_met = data(index_names, :); % Extract rows corresponding to metabolic genes

% Normalize transcriptomics data
data_f = data_met(:, 2:end); % Extract transcriptomics data (excluding explanatory columns)
logdata = log10(data_f + 1); % Normalize data with log10 transformation

%%% Gini coefficient calculation %%%

% Calculate the Gini coefficient using the localgini function
gini = localgini_core(logdata);

% Identify core genes of metabolism based on Gini coefficients
core_indices = find(mean(gini, 2) > 0.5); % Indices of core genes
gene_names = data_met{:, 1}; % Names of genes
core_genes = gene_names(core_indices); % Names of core genes
core_genes_tab = table(core_genes); % Convert to table

%%% Housekeeping genes %%%

% Identify housekeeping genes with metabolic functions
genes_table = data_met.Ensembl_GeneID;
index_names = ismember(genes_table, h_k_g.converted_alias); % Indices of housekeeping genes in the metabolic gene list
hkg_met = data_met(index_names, :); % Extract rows corresponding to housekeeping genes with metabolic functions

% Identify core housekeeping genes
correctly_identified = intersect(core_genes_tab.core_genes, hkg_met.Ensembl_GeneID); % Core housekeeping genes
percentage_correct = (length(correctly_identified) / length(hkg_met.Ensembl_GeneID)) * 100; % Percentage of correctly identified core housekeeping genes
disp(['Correctly identified housekeeping genes: ', num2str(percentage_correct), '%']);

%%% Metabolic housekeeping genes %%%

% Calculate coverage of metabolic housekeeping genes for each sample
sel_rows = gini(index_names == 1, :); % Gini coefficients for metabolic housekeeping genes
num_ones = sum(sel_rows, 1); % Count of ones in each column
num_rows = size(sel_rows, 1); % Total number of rows
prop_genes = num_ones / num_rows; % Proportion of ones in each column

% Create a table with the results
expression_col = data_met.Properties.VariableNames(2:end);
results_col_FPKM = array2table(prop_genes, 'VariableNames', expression_col)

%----------------------- Conversion from FPKM to TPM -----------------------

% Reload the data
data = readtable('Mod_data.xlsx');
h_k_g = readtable('housekeeping_ens.csv');
model = readCbModel('Human-GEM_Cobra_v1.01.mat');

% Convert FPKM to TPM
data_matrix = data{:, 2:end};
column_sums = sum(data_matrix, 1);
normalized_matrix = (data_matrix ./ column_sums) * 1e6;
normalized_table = array2table(normalized_matrix, 'VariableNames', data.Properties.VariableNames(2:end));
data(:, 2:end) = normalized_table;

% Repeat metabolic gene extraction, normalization, Gini coefficient calculation, and core gene identification for TPM data
% (same steps as before, with TPM data instead of FPKM data)
index_names = ismember(data.Ensembl_GeneID, model_genes);
data_met = data(index_names, :);
data_f = data_met(:, 2:end);
logdata = log10(data_f + 1);
gini = localgini_core(logdata);
core_indices = find(mean(gini, 2) > 0.5);
gene_names = data_met{:, 1};
core_genes = gene_names(core_indices);
core_genes_tab = table(core_genes);
index_names = ismember(genes_table, h_k_g.converted_alias);
hkg_met = data_met(index_names, :);
correctly_identified = intersect(core_genes_tab.core_genes, hkg_met.Ensembl_GeneID);
percentage_correct = (length(correctly_identified) / length(hkg_met.Ensembl_GeneID)) * 100;
disp(['Correctly identified housekeeping genes: ', num2str(percentage_correct), '%']);
sel_rows = gini(index_names == 1, :);
num_ones = sum(sel_rows, 1);
num_rows = size(sel_rows, 1);
prop_genes = num_ones / num_rows;
expression_col = data_met.Properties.VariableNames(2:end);
results_col_TPM = array2table(prop_genes, 'VariableNames', expression_col)
