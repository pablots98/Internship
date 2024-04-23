function [onlyOR, anyAND, singleGene, onlyORgenes, anyANDgenes, singleGenes] = analyzeGPRWithExpression(model, expressionData)
    onlyOR = 0;
    anyAND = 0;
    singleGene = 0;
    onlyORgenes = {};  % Almacena los ENSEMBL_ID para only OR
    anyANDgenes = {};  % Almacena los ENSEMBL_ID para any AND
    singleGenes = {};  % Almacena los ENSEMBL_ID para singleGenes
    
    % Assume that expressionData is a table with 'EnsemblID' as one of its columns
    ensemblIDs = expressionData.Ensembl_GeneID;
    
    % Assume that the model has a gene correspondence that also uses Ensembl IDs
    modelGenes = model.genes; % adjust this according to your specific model
    
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
                onlyORgenes = [onlyORgenes; commonGenesInRule];  % Agrega los genes a la lista
            end
            
            if hasAND
                anyAND = anyAND + 1;
                anyANDgenes = [anyANDgenes; commonGenesInRule];  % Agrega los genes a la lista
            end
            
            if ~hasAND && ~hasOR
                singleGene = singleGene + 1;
                singleGenes = [singleGenes; commonGenesInRule];  % Agrega los genes a la lista
            end
        end
    end
end

function genes = extractGenesFromRule(rule)
    % This function needs to be defined to correctly extract genes from a GPR rule
    % depending on the specific format of the rules in your model
    genes = regexp(rule, '[a-zA-Z0-9_]+', 'match');
end
