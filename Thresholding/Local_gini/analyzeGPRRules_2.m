function geneRulesCount = analyzeGPRRules(ensemblIDs, grRules)
    % Inicializa la estructura de salida
    % Initialize the output structure
    geneRulesCount = struct();

    % Recorre cada gen proporcionado
    % Iterate over each provided gene
    for i = 1:length(ensemblIDs)
        geneID = ensemblIDs{i};
        geneRulesCount.(geneID) = struct('single', 0, 'and', 0, 'or', 0, 'and_or', 0);

        % Busca cada regla que contenga el gen
        % Search each rule that contains the gene
        for j = 1:length(grRules)
            rule = grRules{j};
            if contains(rule, geneID)
                % Cuenta el tipo de regla
                % Count the type of rule
                if contains(rule, 'and') && contains(rule, 'or')
                    geneRulesCount.(geneID).and_or = geneRulesCount.(geneID).and_or + 1;
                elseif contains(rule, 'and')

        end
    end
end