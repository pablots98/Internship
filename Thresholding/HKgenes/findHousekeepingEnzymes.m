function [hkDetails] = findHousekeepingDetails(model, hkGenes)
    % USAGE:
    %    [hkDetails] = findHousekeepingDetails(model, hkGenes)
    %
    % INPUT:
    %    model:      A COBRA model with gene-protein-reaction (GPR) associations
    %    hkGenes:    A cell array of housekeeping gene IDs
    %
    % OUTPUT:
    %    hkDetails:  A structure containing detailed housekeeping gene associations

    % Initialize output
    hkDetails = struct('geneAssociations', {}, 'rxnID', {}, 'subSystems', {});

    % Initialize index for structure array
    index = 1;

    % Process each GPR to find housekeeping genes and their associations
    for i = 1:length(model.grRules)
        currentGPR = model.grRules{i};
        if isempty(currentGPR)
            continue;
        end

        % Check if any housekeeping gene is mentioned in the GPR
        containsHKGene = any(cellfun(@(gene) any(contains(currentGPR, gene)), hkGenes));

        if containsHKGene
            % Housekeeping gene found within GPR, add details to the list
            hkDetails(index).geneAssociations = currentGPR;  % Original GPR expression
            hkDetails(index).rxnID = model.rxns{i};          % Associated reaction
            if isfield(model, 'subSystems') && length(model.subSystems) >= i
                hkDetails(index).subSystems = model.subSystems{i};
            else
                hkDetails(index).subSystems = 'Subsystem not available';
            end
            index = index + 1;  % Increment index for next entry
        end
    end
end

