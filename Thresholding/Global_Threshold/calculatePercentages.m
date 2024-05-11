function percentages = calculatePercentages(geneRulesCount)
    % Initialize total counters
    totalSingle = 0;
    totalAnd = 0;
    totalOr = 0;
    totalAndOr = 0;
    
    % Get the list of genes
    genes = fieldnames(geneRulesCount);
    
    % Sum the values for each type of rule
    for i = 1:numel(genes)
        gene = genes{i};
        totalSingle = totalSingle + geneRulesCount.(gene).single;
        totalAnd = totalAnd + geneRulesCount.(gene).and;
        totalOr = totalOr + geneRulesCount.(gene).or;
        totalAndOr = totalAndOr + geneRulesCount.(gene).and_or;
    end
    
    % Total rules
    totalRules = totalSingle + totalAnd + totalOr + totalAndOr;
    
    % Check if the total rules is zero to avoid division by zero
    if totalRules == 0
        warning('No rules found containing the specified genes.');
        percentages.single = NaN;
        percentages.and = NaN;
        percentages.or = NaN;
        percentages.and_or = NaN;
    else
        % Calculate percentages
        percentages = struct();
        percentages.single = (totalSingle / totalRules) * 100;
        percentages.and = (totalAnd / totalRules) * 100;
        percentages.or = (totalOr / totalRules) * 100;
        percentages.and_or = (totalAndOr / totalRules) * 100;
    end
end

