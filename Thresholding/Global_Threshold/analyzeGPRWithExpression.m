function [onlyOR, anyAND, singleGene] = analyzeGPRWithExpression(model, expressionData)
    onlyOR = 0;
    anyAND = 0;
    singleGene = 0;
    
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
        if ~isempty(intersect(genesInRule, commonGenes))
        
            simplifiedRule = lower(strtrim(rule));
            hasOR = contains(simplifiedRule, 'or');
            hasAND = contains(simplifiedRule, 'and');
            
            if hasOR && ~hasAND
                onlyOR = onlyOR + 1;
            end
            
            if hasAND
                anyAND = anyAND + 1;
            end
            
            if ~hasAND && ~hasOR
                singleGene = singleGene + 1;
            end
        end
    end
end

function genes = extractGenesFromRule(rule)
    % This function needs to be defined to correctly extract genes from a GPR rule
    % depending on the specific format of the rules in your model
    genes = regexp(rule, '[a-zA-Z0-9_]+', 'match');
end
