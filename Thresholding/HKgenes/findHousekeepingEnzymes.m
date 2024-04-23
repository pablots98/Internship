function [hkDetails] = findHousekeepingDetails(model, hkGenes, ruleType)
    % USAGE:
    %    [hkDetails] = findHousekeepingDetails(model, hkGenes, ruleType)
    %
    % INPUT:
    %    model:      A COBRA model with gene-protein-reaction (GPR) associations
    %    hkGenes:    A cell array of housekeeping gene IDs
    %    ruleType:   A string specifying the rule type ('AND', 'OR', 'INDIVIDUAL')
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
        containsHKGene = any(cellfun(@(gene) any(contains(currentGPR, gene, 'IgnoreCase', true)), hkGenes));

        % Continue only if a housekeeping gene is found in the GPR
        if containsHKGene
            % Define patterns based on the ruleType
            matchesRuleType = false;
            switch upper(ruleType)
                case 'AND'
                    matchesRuleType = contains(currentGPR, ' and ', 'IgnoreCase', true);
                case 'OR'
                    matchesRuleType = contains(currentGPR, ' or ', 'IgnoreCase', true);
                case 'INDIVIDUAL'
                    % This pattern checks if the GPR contains only the gene ID and possibly whitespace
                    individualPattern = strcat('^[\s]*', strjoin(hkGenes, '$|[\s]*'), '$');
                    matchesRuleType = ~isempty(regexp(currentGPR, individualPattern, 'once'));
                otherwise
                    error('Invalid ruleType specified. Choose from ''AND'', ''OR'', ''INDIVIDUAL''.');
            end

            % Add details to hkDetails if the rule matches the specified type
            if matchesRuleType
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
end
