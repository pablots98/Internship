function HK_rxns = find_housekeeping_reactions(model, HKG, mode)
    %% FIND_HOUSEKEEPING_REACTIONS identifies housekeeping reactions in a metabolic model
    %
    %   HK_rxns = FIND_HOUSEKEEPING_REACTIONS(model, HKG, mode) returns a list of reactions
    %   that are considered housekeeping reactions. The definition of housekeeping
    %   reactions depends on the selected mode.
    %
    %%   Inputs:
    %       model - A metabolic model structure containing fields:
    %           model.grRules - Cell array of gene-reaction rules as strings
    %           model.rxns - Cell array of reaction identifiers
    %       HKG - Cell array of housekeeping gene identifiers
    %       mode - String specifying the mode, either 'strict' or
    %       'lenient'. 
    %           strict - at least one enzyme that catalyzes the reaction
    %                    is composed entirely of housekeeping genes.
    %           lenient - at least one enzyme that catalyzes the reaction
    %                     is composed part of housekeeping genes. (i.e any
    %                     reaction as long as at least one HK gene appears in
    %                     grRules.)
    %
    %%   Output:
    %       HK_rxns - Cell array of reaction identifiers that are housekeeping reactions
    %
    %%  Authur: ChatGPT，Shuyi
    %
    %%  Test code
    % model.grRules = {'ENSG00000123415 and ENSG00000123416', 'ENSG00000123417 or (ENSG00000123418 and ENSG00000123419)', 'ENSG00000123420'};
    % model.rxns = {'RXN1', 'RXN2', 'RXN3'};
    % HKG = {'ENSG00000123415', 'ENSG00000123416', 'ENSG00000123418', 'ENSG00000123420'};
    % HK_rxns = find_housekeeping_reactions(model, HKG,'strict');
    % disp(HK_rxns);
    %%  Usage：HK_rxns_met = find_housekeeping_reactions(model, HKG,'strict');
    %
    %% Begin
    % Validate the input mode
    if ~ismember(mode, {'strict', 'lenient'})
        error('Invalid mode. Choose either ''strict'' or ''lenient''.');
    end
    
    % Initialize the result as an empty cell array
    HK_rxns = {};
    
    % Iterate over each reaction in the model
    for i = 1:length(model.grRules)
        grRule = model.grRules{i}; % Get the gene-reaction rule for the current reaction
        
        % Parse the grRule to get individual enzymes
        enzymes = strsplit(grRule, ' or '); % Split the rule by 'or' to get different enzyme options
        
        % Flag to check if this reaction is a housekeeping reaction
        is_HK_rxn = false; % Initialize flag as false
        
        % Check each enzyme
        for j = 1:length(enzymes)
            enzyme = enzymes{j}; % Get the current enzyme
            
            % Remove parentheses and split into individual genes
            enzyme = strrep(enzyme, '(', ''); % Remove left parenthesis
            enzyme = strrep(enzyme, ')', ''); % Remove right parenthesis
            genes = strsplit(enzyme, ' and '); % Split the enzyme into individual genes by 'and'
            
            % Check the mode
            switch mode
                case 'strict'
                    % In strict mode, all genes in this enzyme must be housekeeping genes
                    if all(ismember(genes, HKG)) % If all genes are in the HKG list
                        is_HK_rxn = true; % Set the flag to true
                        break; % No need to check other enzymes for this reaction
                    end
                case 'lenient'
                    % In lenient mode, if any gene in this enzyme is a housekeeping gene
                    if any(ismember(genes, HKG)) % If any gene is in the HKG list
                        is_HK_rxn = true; % Set the flag to true
                        break; % No need to check other enzymes for this reaction
                    end
            end
        end
        
        % If this is a housekeeping reaction, add it to the result
        if is_HK_rxn
            HK_rxns{end+1} = model.rxns{i}; % Add the reaction identifier to the result list
        end
    end
end