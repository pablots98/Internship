% Definir la función gini
function G = gini(x)
    % Calcular el coeficiente de Gini de un vector x
    assert(numel(x) > 1, 'gini expects at least two elements in the vector.');

    % Ordenar el vector en orden ascendente
    x = sort(x);

    % Número de elementos en el vector
    n = numel(x);

    % Calcular el coeficiente de Gini
    G = (2 * (1:n) * x(:)) / (n * sum(x)) - (n + 1) / n;
end