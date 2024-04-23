function [g, l, a] = giniUniform(val, makeplot)
    % Assign a uniform population vector
    pop = ones(size(val));
    [g, l, a] = gini(pop, val, makeplot);
end
