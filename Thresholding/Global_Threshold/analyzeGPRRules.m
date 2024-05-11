function geneRulesCount = analyzeGPRRules(ensemblIDs, grRules)
    % Initialize the output structure
    geneRulesCount = struct();
    
    % Iterate over each provided gene
    for i = 1:length(ensemblIDs)
        geneID = ensemblIDs{i};
        geneRulesCount.(geneID) = struct('single', 0, 'and', 0, 'or', 0, 'and_or', 0);
        
        % Search each rule that contains the gene
        for j = 1:length(grRules)
            rule = grRules{j};
            if contains(rule, geneID)
                % Count the type of rule
                if contains(rule, 'and') && contains(rule, 'or')
                    geneRulesCount.(geneID).and_or = geneRulesCount.(geneID).and_or + 1;
                elseif contains(rule, 'and')
                    geneRulesCount.(geneID).and = geneRulesCount.(geneID).and + 1;
                elseif contains(rule, 'or')
                    geneRulesCount.(geneID).or = geneRulesCount.(geneID).or + 1;
                else
                    geneRulesCount.(geneID).single = geneRulesCount.(geneID).single + 1;
                end
            end
        end
    end
end

