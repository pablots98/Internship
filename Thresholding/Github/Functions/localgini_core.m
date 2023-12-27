function [coreMat] = localgini_core(mappedDat)
%LOCALT2_CORE Uses localgini to define the core set of reactions based on Gini coefficient
%   INPUTS:     
%       mappedDat, a matrix of expression values (gene/reaction X Sample). Usually gene
%               expression values or gene expression values mapped onto a models
%               reactions using GPR rules.
%
%   OUTPUTS:    
%       coreMat, the binary core set matrix (gene/reaction X Sample).

% Convierte la tabla en una matriz
if istable(mappedDat)
    mappedDat = table2array(mappedDat);
end

% Compute the Gini coefficients for each gene/reaction
giniCoefficients = zeros(size(mappedDat, 1), 1);
for i = 1:size(mappedDat, 1)
    sortedRow = sort(mappedDat(i, :), 'ascend');
    rowSum = sum(sortedRow);
    cumulativeSum = cumsum(sortedRow);
    giniCoefficients(i) = 1 - 2 * (sum((2*(1:size(mappedDat, 2)) - size(mappedDat, 2) - 1) .* sortedRow) / (size(mappedDat, 2) * rowSum));
end

% Create coreMat based on Gini coefficients
coreMat = false(size(mappedDat));
for i = 1:size(mappedDat, 1)
    coreMat(i, :) = mappedDat(i, :) > giniCoefficients(i);
end
end
