function [coreRxns_binary, reactionFractions] = ubiquityScore_Thresholding(expressionRxns_Sample, sampleNames, model_p)
    % Create the initial matrix with appropriate dimensions and zeros
    row_names = model_p.rxns;
    coreRxns_binary = zeros(length(row_names), length(sampleNames));
    coreRxns_binary = array2table(coreRxns_binary, 'RowNames', row_names, 'VariableNames', sampleNames);

    % Populate the table with data from expressionRxns_Sample
    for i = 1:length(sampleNames)
        % Check that the length of the vector in expressionRxns_Sample{i} is correct
        if length(expressionRxns_Sample{i}) == size(coreRxns_binary, 1)
            coreRxns_binary{:, i} = expressionRxns_Sample{i};
        else
            error('The size of the expression vector does not match the number of reactions');
        end
    end

    % Initialize the matrix for coreRxns_binary
    for i = 1:width(coreRxns_binary)
        current_column = coreRxns_binary{:, i};
        
        % Values greater than 1 are set to 1
        current_column(current_column > 1) = 1;
        
        % Values equal to 0 are set to -1e-6
        current_column(current_column == 0) = -1e-6;
        
        % Convert the modified column back to the table
        coreRxns_binary{:, i} = current_column;
    end

    % Calculate the median ubiquity score for non-zero non-core values for NaN replacement
    non_core_values = table2array(coreRxns_binary);
    non_core_values = non_core_values(non_core_values > 0 & non_core_values < 1);
    median_non_core = median(non_core_values, 'omitnan'); % Omitting NaN values for the median calculation

    % Replace NaN values with the calculated median ubiquity score
    for i = 1:width(coreRxns_binary)
        current_column = coreRxns_binary{:, i};
        current_column(isnan(current_column)) = median_non_core;
        coreRxns_binary{:, i} = current_column;
    end

    % Calculate the fraction of 1s for each reaction across all samples
    reactionFractions = sum(table2array(coreRxns_binary) == 1, 2) / width(coreRxns_binary);

    % Return the updated matrix and reaction fractions
    return
end

