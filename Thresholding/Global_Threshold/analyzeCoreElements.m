function results = analyzeCoreElements(GenExp_all, GeneExp_P, met_ens, model_p)
    % Initialize storage structures
    results = struct();

    % Necessary data lengths
    IMR90_Y2_length = length(GeneExp_P);
    IMR90_Y2_RLength = length(model_p.rxns);

    % Iterate from the 100th percentile down to the 60th
    for p = 100:-1:1
        percentile_value = prctile(GenExp_all, p);
        GeneExp_P_norm = GeneExp_P ./ percentile_value;
        coreGenes = met_ens(GeneExp_P_norm >= 1);
        expressionData.gene = met_ens;
        expressionData.value = GeneExp_P_norm;
        IMR90_Y2_normCR = mapExpressionToReactions(model_p, expressionData, false);
        coreReactions = model_p.rxns(IMR90_Y2_normCR >= 1);
        ruleAnalysis = analyzeGPR(coreReactions, model_p.grRules, model_p);

        % Save results including percentages
        results.(sprintf('P%d', p)) = struct(...
            'threshold', percentile_value, ...
            'numCoreGenes', length(coreGenes), ...
            'numCoreReactions', length(coreReactions), ...
            'numOrRules', ruleAnalysis.or, ...
            'numAndRules', ruleAnalysis.and, ...
            'numSingleRules', ruleAnalysis.single, ...
            'numCombinedRules', ruleAnalysis.combined, ...
            'percentOrRules', ruleAnalysis.orPct, ...
            'percentAndRules', ruleAnalysis.andPct, ...
            'percentSingleRules', ruleAnalysis.singlePct, ...
            'percentCombinedRules', ruleAnalysis.combinedPct);
    end
end


function ruleAnalysis = analyzeGPR(coreReactions, grRules, model_p)
    ruleAnalysis = struct('or', 0, 'and', 0, 'single', 0, 'combined', 0);
    totalRules = length(coreReactions);  % Store the total number of core reactions

    for i = 1:totalRules
        reaction = coreReactions(i);
        rule = grRules{strcmp(model_p.rxns, reaction)};
        if contains(rule, 'and') && contains(rule, 'or')
            ruleAnalysis.combined = ruleAnalysis.combined + 1;
        elseif contains(rule, 'and')
            ruleAnalysis.and = ruleAnalysis.and + 1;
        elseif contains(rule, 'or')
            ruleAnalysis.or = ruleAnalysis.or + 1;
        else
            ruleAnalysis.single = ruleAnalysis.single + 1;
        end
    end

    % Calculate percentages for each type of rule
    if totalRules > 0
        ruleAnalysis.orPct = 100 * ruleAnalysis.or / totalRules;
        ruleAnalysis.andPct = 100 * ruleAnalysis.and / totalRules;
        ruleAnalysis.singlePct = 100 * ruleAnalysis.single / totalRules;
        ruleAnalysis.combinedPct = 100 * ruleAnalysis.combined / totalRules;
    else
        ruleAnalysis.orPct = 0;
        ruleAnalysis.andPct = 0;
        ruleAnalysis.singlePct = 0;
        ruleAnalysis.combinedPct = 0;
    end
end



