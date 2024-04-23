function gini = giniCoefficient(x)
    x = x(x>0); % Considerando solo valores positivos para el cálculo
    n = length(x);
    if n < 2
        gini = 0; % No se puede calcular Gini con un solo dato o sin datos
        return;
    end
    x = sort(x);
    sumx = sum(x);
    sumx2 = sum((2 * (1:n) - n - 1) .* x);
    gini = 1 - 2 * sumx2 / (n * sumx);

    % Asegurarse que el resultado esté en el rango [0, 1]
    gini = max(0, min(1, gini));
end
