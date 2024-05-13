function [coreRxns_binary, reactionFractions] = ubiquityScore_Thresholding(Rxns_local25_75, sampleNames, row_names)
    % Apply threshold to the Rxns_local25_75 matrix
    % Values >= 1 are converted to 1, and values < 1 are converted to 0
    binaryData = Rxns_local25_75 >= 1;

    % Convert the binary data matrix into a table with sample and reaction names
    coreRxns_binary = array2table(binaryData, 'RowNames', row_names, 'VariableNames', sampleNames);

    % Calculate the fraction of 1s for each reaction across all samples
    reactionFractions = sum(binaryData, 2) / size(binaryData, 2);
    reactionFractions = array2table(reactionFractions, 'RowNames', row_names, 'VariableNames', {'FractionActive'});

    % Return the processed matrices as tables
    return
end


