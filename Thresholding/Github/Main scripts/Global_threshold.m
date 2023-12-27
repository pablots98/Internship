%%%%%%%%%%%%%%%%%% Global Thresholding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%%% Load Data %%%
% Load transcriptomics data, housekeeping genes, and metabolic model
data = readtable('Mod_data.xlsx');
h_k_g = readtable('housekeeping_ens.csv');
model = readCbModel('Human-GEM_Cobra_v1.01.mat');

%%% Metabolic Genes %%%
% Extract genes related to metabolism
model_genes = model.genes;
index_names = ismember(data.Ensembl_GeneID, model_genes);
data_met = data(index_names, :);

%%% Calculate Percentage for FPKM %%%
results = table;
expression_col = data_met.Properties.VariableNames(2:end);

for i = 1:length(expression_col)
    col_name = expression_col{i};
    expfibroblast = data_met(:,1);
    expfibroblast(:,2) = data_met(:, col_name);
    expfibroblast.Properties.VariableNames{1} = 'gene';
    expfibroblast.Properties.VariableNames{2} = 'FPKMvalue';

    % Log-transform and normalize expression data
    expfibroblast.logFPKMvalue = log10(expfibroblast{:,2} + 1);
    expfibroblast.value = expfibroblast.logFPKMvalue - min(expfibroblast.logFPKMvalue);

    % Calculate threshold and filter data
    percentage = 70;
    up_threshold_fib = prctile(expfibroblast.value, percentage);
    idx = expfibroblast.value >= up_threshold_fib;
    filtered_data = expfibroblast(idx, :);

    % Identify housekeeping genes
    index_names = ismember(data_met.Ensembl_GeneID, h_k_g.converted_alias);
    hkg_met = data_met(index_names, :);
    hkg_met_ens = hkg_met(:, "Ensembl_GeneID");

    matching_names_idx = ismember(filtered_data.gene, hkg_met_ens.Ensembl_GeneID);
    matching_names = filtered_data(matching_names_idx, 1);

    % Calculate percentage of maintenance genes
    percentage_hk = (length(matching_names.gene) / length(filtered_data.gene)) * 100;

    % Save results
    results{i, 'Column'} = {col_name};
    results{i, 'PercentageHK'} = percentage_hk;
end
results_FPKM = results
disp(results_FPKM);

%%% Convert FPKM to TPM %%%
data_matrix = data{:, 2:end};
column_sums = sum(data_matrix, 1);
normalized_matrix = (data_matrix ./ column_sums) * 1e6;
normalized_table = array2table(normalized_matrix, 'VariableNames', data.Properties.VariableNames(2:end));
data(:, 2:end) = normalized_table;

%%% Recalculate Percentage for TPM %%%
results = table;
expression_col = data_met.Properties.VariableNames(2:end);

for i = 1:length(expression_col)
    col_name = expression_col{i};
    expfibroblast = data_met(:,1);
    expfibroblast(:,2) = data_met(:, col_name);
    expfibroblast.Properties.VariableNames{1} = 'gene';
    expfibroblast.Properties.VariableNames{2} = 'FPKMvalue';

    % Log-transform and normalize expression data
    expfibroblast.logFPKMvalue = log10(expfibroblast{:,2} + 1);
    expfibroblast.value = expfibroblast.logFPKMvalue - min(expfibroblast.logFPKMvalue);

    % Calculate threshold and filter data
    percentage = 70;
    up_threshold_fib = prctile(expfibroblast.value, percentage);
    idx = expfibroblast.value >= up_threshold_fib;
    filtered_data = expfibroblast(idx, :);

    % Identify housekeeping genes
    index_names = ismember(data_met.Ensembl_GeneID, h_k_g.converted_alias);
    hkg_met = data_met(index_names, :);
    hkg_met_ens = hkg_met(:, "Ensembl_GeneID");

    matching_names_idx = ismember(filtered_data.gene, hkg_met_ens.Ensembl_GeneID);
    matching_names = filtered_data(matching_names_idx, 1);

    % Calculate percentage of maintenance genes
    percentage_hk = (length(matching_names.gene) / length(filtered_data.gene)) * 100;

    % Save results
    results{i, 'Column'} = {col_name};
    results{i, 'PercentageHK'} = percentage_hk;
end
results_TPM = results
disp(results_TPM);
