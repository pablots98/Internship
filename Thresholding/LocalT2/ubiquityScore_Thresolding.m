function [coreRxns_binary, reactionFractions] = ubiquityScore_Thresholding(Rxns_local25_75)
    % Aplicar el umbral a la matriz Rxns_local25_75
    % Los valores >= 1 se convierten en 1, y los valores < 1 en 0
    coreRxns_binary = Rxns_local25_75 >= 1;

    % Calcular la fracción de 1s para cada reacción a lo largo de todas las muestras
    reactionFractions = sum(coreRxns_binary, 2) / size(coreRxns_binary, 2);

    % Devolver las matrices procesadas
    return
end

