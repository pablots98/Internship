function statisticalResults = calculateGeneExpressionStats(metabolicGenesThing)
    % Total number of samples
    numSamples = width(metabolicGenesThing);

    % Extract sample names
    sampleNames = metabolicGenesThing.Properties.VariableNames;

    % Initialize arrays to store the means and standard deviations
    means = zeros(numSamples, 1);
    stdDevs = zeros(numSamples, 1);

    % Iterate over each sample (column)
    for i = 1:numSamples
        % Get gene expression data for the current sample
        geneExpr = metabolicGenesThing{:, i};

        % Filter the data to include only those that are >= 0
        filteredGeneExpr = geneExpr(geneExpr >= 0);

        % Calculate the mean and standard deviation
        means(i) = mean(filteredGeneExpr);
        stdDevs(i) = std(filteredGeneExpr);
    end

    % Create a table with the results
    statisticalResults = table(sampleNames', means, stdDevs, 'VariableNames', {'Sample', 'Mean', 'StdDev'});
end
