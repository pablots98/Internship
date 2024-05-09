function geneRulesCount = analyzeGPRRules(ensemblIDs, grRules)
    % Inicializa la estructura de salida
    geneRulesCount = struct();
    
    % Recorre cada gen proporcionado
    for i = 1:length(ensemblIDs)
        geneID = ensemblIDs{i};
        geneRulesCount.(geneID) = struct('single', 0, 'and', 0, 'or', 0, 'and_or', 0);
        
        % Busca cada regla que contenga el gen
        for j = 1:length(grRules)
            rule = grRules{j};
            if contains(rule, geneID)
                % Cuenta el tipo de regla
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
