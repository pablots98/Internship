function [adjusted_matrix, expression_scoreMatrix] = localT2_new(mappedDat, lowerThres, upperThres)
% This function adjusts gene expression data based on provided lower and 
% upper percentile thresholds. It converts the input data into a binary 
% format where values equal to or greater than 1 are set to 1, and values 
% less than 1 are set to 0. This binary adjustment is based on the 
% comparison of each gene's expression level against the calculated thresholds.
% 
% Inputs:
% - mappedDat: A table containing gene expression data. Each row represents
% a different gene, and each column represents different samples or conditions.
% - lowerThres (optional): The lower percentile threshold for adjustment. 
% If not provided, it defaults to 25.
% - upperThres (optional): The upper percentile threshold for adjustment. 
% If not provided, it defaults to 75.
% Outputs:
% - adjusted_matrix: A binary matrix of the same size as mappedDat. Values 
% are set to 1 if the original expression value is equal to or exceeds the 
% adjusted threshold, and 0 otherwise.
%
% Reference (Richelle et al. 2019, doi: 10.1371/JOURNAL.PCBI.1007185)

    
    % Check for the existence of thresholds and issue a warning if necessary
    if exist('lowerThres', 'var') && exist('upperThres', 'var') && lowerThres < 1 && upperThres < 1
        warning('The thresholds should be percentiles in the 0 to 100 range, you may be using the 0 to 1 range here.')
    end
    if ~exist('lowerThres', 'var')
        lowerThres = 25;
    end
    if ~exist('upperThres', 'var')
        upperThres = 75;
    end

    % Convert the table to a matrix
    mappedDat = table2array(mappedDat);

    linData = reshape(mappedDat, [], 1);
    lowerThres = quantile(linData, lowerThres / 100);
    upperThres = quantile(linData, upperThres / 100);

    coreMat = false(size(mappedDat));

    % Initialize the vector for adjusted thresholds
    expression_ths_local_T2_25_75 = zeros(size(mappedDat, 1), 1);

    % Calculate the adjusted thresholds for each gene
    for i = 1:size(mappedDat, 1)
        expressionValue = mappedDat(i, :);
        if mean(expressionValue) >= upperThres
            expression_ths_local_T2_25_75(i) = upperThres;
        else
            expression_ths_local_T2_25_75(i) = max(mean(expressionValue), lowerThres);
        end
    end

    % Adjust expression values
    expression_scoreLocal_T2_25_75 = zeros(size(mappedDat));
    for i = 1:size(mappedDat, 2)
        expression_scoreLocal_T2_25_75(:, i) = mappedDat(:, i) ./ expression_ths_local_T2_25_75;
    end

    % Convert the adjusted values to 1 if they are >= 1, otherwise to 0
    adjusted_matrix = expression_scoreLocal_T2_25_75 %>= 1;
    %adjusted_matrix(adjusted_matrix >= 1) = 1
    % New matrix with values between 0 and 1
    expression_scoreMatrix = expression_scoreLocal_T2_25_75 >= 1;
    %expression_scoreMatrix(expression_scoreMatrix >= 1) = 1;

end
