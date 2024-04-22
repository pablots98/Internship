function [expressionOR, expressionAND, expressionSingle] = extractExpressionByGPR(model, expressionData)
    % Carga los datos de expresión desde un archivo o una tabla ya en el entorno
    % expressionData = readtable('expression_data.csv'); % Descomentar si necesitas cargar los datos

    % Inicializa listas para almacenar expresiones según el tipo de regla
    expressionOR = [];
    expressionAND = [];
    expressionSingle = [];

    % Procesa cada reacción en el modelo
    for i = 1:length(model.rxns)
        rule = model.grRules{i};
        if isempty(rule)
            continue;
        end

        % Extrae genes de la regla
        genesInRule = extractGenesFromRule(rule);

        % Encuentra los niveles de expresión para estos genes
        idx = ismember(expressionData.Ensembl_GeneID, genesInRule);
        geneExpressions = expressionData.Expression(idx);

        % Determina el tipo de regla y almacena la expresión correspondiente
        simplifiedRule = lower(strtrim(rule));
        hasOR = contains(simplifiedRule, 'or');
        hasAND = contains(simplifiedRule, 'and');

        if hasOR && ~hasAND
            expressionOR = [expressionOR; geneExpressions];
        elseif hasAND
            expressionAND = [expressionAND; geneExpressions];
        elseif ~hasAND && ~hasOR
            expressionSingle = [expressionSingle; geneExpressions];
        end
    end
end
