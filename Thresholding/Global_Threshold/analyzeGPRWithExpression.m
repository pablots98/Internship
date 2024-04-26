function [onlyOR, onlyAND, singleGene, mixed, onlyORgenes, onlyANDgenes, singleGenes, mixedGenes] = analyzeGPRWithExpression(model, expressionData)
    onlyOR = 0;
    onlyAND = 0;
    singleGene = 0;
    mixed = 0;
    onlyORgenes = {};  
    onlyANDgenes = {};  
    singleGenes = {};
    mixedGenes = {};
    
    ensemblIDs = expressionData.Ensembl_GeneID;
    modelGenes = model.genes;  % adjust this according to your specific model
    
    % Find common genes
    commonGenes = intersect(ensemblIDs, modelGenes);
    
    for i = 1:length(model.rxns)
        rule = model.grRules{i};
        if isempty(rule)
            continue;
        end
        
        % Extract the genes mentioned in the rule
        genesInRule = extractGenesFromRule(rule);
        
        % Check if any gene in the rule is on the list of common genes
        commonGenesInRule = intersect(genesInRule, commonGenes);
        if ~isempty(commonGenesInRule)
        
            simplifiedRule = lower(strtrim(rule));
            hasOR = contains(simplifiedRule, 'or');
            hasAND = contains(simplifiedRule, 'and');
            
            if hasOR && ~hasAND
                onlyOR = onlyOR + 1;
                onlyORgenes = [onlyORgenes; commonGenesInRule];
            elseif hasAND && ~hasOR
                onlyAND = onlyAND + 1;
                onlyANDgenes = [onlyANDgenes; commonGenesInRule];
            elseif hasAND && hasOR
                mixed = mixed + 1;
                mixedGenes = [mixedGenes; commonGenesInRule];  % Sumar aquellos que contienen ambos
            elseif ~hasAND && ~hasOR
                singleGene = singleGene + 1;
                singleGenes = [singleGenes; commonGenesInRule];
            end
        end
    end
end
