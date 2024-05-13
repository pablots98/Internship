function [coreRxns_binary, reactionFractions] = ubiquityScore_Thresolding(expressionRxns_Sample, sampleNames, model_p)
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

    % Convert values to 1 or 0 based on the specified condition
    for i = 1:width(coreRxns_binary)
        coreRxns_binary{:, i} = coreRxns_binary{:, i} >= 1;
    end

    % Calculate the fraction of 1s for each reaction across all samples
    reactionFractions = sum(table2array(coreRxns_binary), 2) / width(coreRxns_binary);

    % Return the updated matrix
    return
end

