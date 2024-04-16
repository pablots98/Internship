function [enzymesInfo] = getEnzymesByGenes(model, geneList)
    % Parse and organize GPRs
    parsedGPR = GPRparser(model);
    [parsedGPR, ix] = linearization_index(parsedGPR, 'rows');
    
    % Ensure that ix does not exceed the number of genes
    ix = ix(ix <= length(model.genes));
    
    corrRxns = model.rxns(ix);
    corrSys = model.subSystems(ix);
    modelGenes = model.genes(ix);  % This should no longer give an error

    % Initialization of structures to store enzyme information
    enzymes = {};
    rxns = {};
    subSystems = {};

    % Iterate through each gene in the input list
    for i = 1:length(geneList)
        currentGene = geneList{i};
        
        % Find indices where the current gene is involved in GPR
        geneInGPR = find(contains(modelGenes, currentGene));
        
        if ~isempty(geneInGPR)
            % For each reaction found, search for the associated enzyme according to the GPR
            for j = 1:length(geneInGPR)
                enzyme = parsedGPR{geneInGPR(j)};
                if iscell(enzyme)
                    enzyme = char(enzyme{1});  % Assuming the cell is not empty
                else
                    enzyme = char(enzyme);  % Explicitly convert to char
                end
                enzymes{end+1} = enzyme;
                rxns{end+1} = corrRxns{geneInGPR(j)};
                subSystems{end+1} = corrSys{geneInGPR(j)};
            end
        end
    end
    
    % Remove duplicates and maintain only unique enzyme identifiers
    [enzymes, uniqueIdx] = unique(enzymes);
    rxns = rxns(uniqueIdx);
    subSystems = subSystems(uniqueIdx);
    
    % Create output structure
    enzymesInfo.enzymes = enzymes;
    enzymesInfo.rxns = rxns;
    enzymesInfo.subSystems = subSystems;
end






