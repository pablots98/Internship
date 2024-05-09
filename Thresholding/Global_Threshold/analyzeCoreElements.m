function results = analyzeCoreElements(GenExp_all, GeneExp_P, met_ens, model_p)
    % Inicializa las estructuras de almacenamiento
    results = struct();

    % Longitudes de datos necesarias
    IMR90_Y2_length = length(GeneExp_P);
    IMR90_Y2_RLength = length(model_p.rxns); % Supone una lista de todas las reacciones en el modelo

    % Itera desde el percentil 100 hasta el 60
    for p = 100:-1:60
        % Calcula el umbral de expresión
        percentile_value = prctile(GenExp_all, p);
        GeneExp_P_norm = GeneExp_P ./ percentile_value;

        % Determina los core genes
        coreGenes = met_ens(GeneExp_P_norm >= 1);

        % Prepara los datos de expresión para el mapeo
        expressionData.gene = met_ens;
        expressionData.value = GeneExp_P_norm;

        % Mapea la expresión genética a reacciones
        IMR90_Y2_normCR = mapExpressionToReactions(model_p, expressionData, false);
        coreReactions = model_p.rxns(IMR90_Y2_normCR >= 1);

        % Análisis de GPR rules
        ruleAnalysis = analyzeGPR(coreReactions, model_p.grRules, model_p);

        % Guarda los resultados para este percentil
        results.(sprintf('P%d', p)) = struct(...
            'threshold', percentile_value, ...
            'numCoreGenes', length(coreGenes), ...
            'numCoreReactions', length(coreReactions), ...
            'numOrRules', ruleAnalysis.or, ...
            'numAndRules', ruleAnalysis.and, ...
            'numSingleRules', ruleAnalysis.single, ...
            'numCombinedRules', ruleAnalysis.combined);

    end
end

function ruleAnalysis = analyzeGPR(coreReactions, grRules, model_p)
    ruleAnalysis = struct('or', 0, 'and', 0, 'single', 0, 'combined', 0);
    for i = 1:length(coreReactions)
        reaction = coreReactions(i);
        rule = grRules{strcmp(model_p.rxns, reaction)};
        if contains(rule, 'and') && contains(rule, 'or')
            ruleAnalysis.combined = ruleAnalysis.combined + 1;
        elseif contains(rule, 'and')
            ruleAnalysis.and = ruleAnalysis.and + 1;
        elseif contains(rule, 'or')
            ruleAnalysis.or = ruleAnalysis.or + 1;
        else
            ruleAnalysis.single = ruleAnalysis.single + 1;
        end
    end
end


