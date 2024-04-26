function genes = extractGenesFromRule(rule)
    % Remove logical operators and parentheses
    cleanedRule = regexprep(rule, 'and|or|not|\(|\)', '');
    % This function needs to be defined to correctly extract genes from a GPR rule
    % depending on the specific format of the rules in your model
    genes = regexp(cleanedRule, '[a-zA-Z0-9_]+', 'match');
end