function [adjusted_matrix, expression_scoreMatrix] = localGini_function(mappedDat, unpruneddata, lowerThres, upperThres)
    % Esta función ajusta los datos de expresión génica en función de los umbrales de percentil proporcionados.
    % Convierte los datos de entrada en un formato binario donde los valores iguales o mayores a 1 se establecen en 1,
    % y los valores menores a 1 se establecen en 0.

    if exist('lowerThres', 'var') && exist('upperThres', 'var') && lowerThres < 1 && upperThres < 1
        warning('The thresholds should be percentiles in the 0 to 100 range, you may be using the 0 to 1 range here.')
    end
    if ~exist('lowerThres', 'var')
        lowerThres = 25;
    end
    if ~exist('upperThres', 'var')
        upperThres = 75;
    end

    mappedDat = table2array(mappedDat);
    unpruneddata = table2array(unpruneddata);

    linData = reshape(unpruneddata(:, 1:end-1), [], 1);
    lowerThres = quantile(linData, lowerThres / 100);
    upperThres = quantile(linData, upperThres / 100);

    coreMat = false(size(mappedDat));
    expression_ths_local_Gini_25_75 = zeros(size(mappedDat, 1), 1);

    for i = 1:size(mappedDat, 1)
        expressionValue = mappedDat(i, :);

        % Calcular el coeficiente de Gini usando la fórmula del paper
        giniCoefficient = gini(expressionValue);

        % Calcular el umbral Localgini
        localGiniThreshold = prctile(expressionValue, giniCoefficient * 100);

        % Asignar el umbral adecuado basado en lowerThres y upperThres
        if localGiniThreshold >= upperThres
            expression_ths_local_Gini_25_75(i) = upperThres;
        elseif localGiniThreshold <= lowerThres
            expression_ths_local_Gini_25_75(i) = lowerThres;
        else
            expression_ths_local_Gini_25_75(i) = localGiniThreshold;
        end
    end

    expression_scoreLocal_Gini_25_75 = zeros(size(mappedDat));
    for i = 1:size(mappedDat, 2)
        expression_scoreLocal_Gini_25_75(:, i) = mappedDat(:, i) ./ expression_ths_local_Gini_25_75;
    end

    expression_scoreMatrix = expression_scoreLocal_Gini_25_75;
    adjusted_matrix = expression_scoreLocal_Gini_25_75 >= 1;
end

% Definir la función gini según la fórmula del paper
function G = gini(x)
    % Calcular el coeficiente de Gini de un vector x
    n = numel(x);
    mean_x = mean(x);
    G = sum(sum(abs(x(:) - x(:)'))) / (2 * n^2 * mean_x);
end
