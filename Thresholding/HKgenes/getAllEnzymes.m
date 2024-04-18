function [enzymes] = getAllEnzymes(model)
    % USAGE:
    %    [enzymes] = getAllEnzymes(model)
    %
    % INPUT:
    %    model:      A COBRA model with gene-protein-reaction (GPR) associations
    %
    % OUTPUT:
    %    enzymes:    A cell array containing enzyme information, including associated reactions and subsystems

    % Check for GPR field in the model
    if ~isfield(model, 'grRules')
        error('Model does not contain gene-protein-reaction associations (grRules).');
    end

    % Extract GPRs, reactions, and subsystems
    grRules = model.grRules;
    reactions = model.rxns;
    subsystems = model.subSystems;

    % Preallocate output structure
    enzymes = struct('names', {}, 'rxns', {}, 'subSystems', {});

    % Iterate through each GPR to construct the enzyme list
    uniqueGPRs = unique(grRules);
    for i = 1:length(uniqueGPRs)
        if isempty(uniqueGPRs{i})
            continue;
        end

        % Find indices where current GPR occurs
        idx = strcmp(grRules, uniqueGPRs{i});
        
        % Collect reactions and subsystems for this GPR
        enzymeRxns = reactions(idx);
        enzymeSubSys = subsystems(idx);

        % Store in structure
        enzymes(i).names = uniqueGPRs{i};
        enzymes(i).rxns = enzymeRxns;
        enzymes(i).subSystems = enzymeSubSys;
    end
end
