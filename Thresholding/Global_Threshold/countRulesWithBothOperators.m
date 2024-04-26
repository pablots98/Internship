function count = countRulesWithBothOperators(grRules)
    count = 0;  % Inicializa el contador

    for i = 1:length(grRules)
        rule = grRules{i};  % Extrae cada regla
        if isempty(rule)
            continue;  % Si la regla está vacía, continua con la siguiente iteración
        end
        
        hasAND = contains(rule, 'and', 'IgnoreCase', true);  % Verifica si contiene 'AND'
        hasOR = contains(rule, 'or', 'IgnoreCase', true);  % Verifica si contiene 'OR'

        if hasAND && hasOR
            count = count + 1;  % Incrementa el contador si la regla contiene ambos
        end
    end
end
