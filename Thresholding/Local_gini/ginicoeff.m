function G = ginicoeff(x)
% Calculate the Gini coefficient of a vector x
n = length(x);
x = sort(x);
G = 2 * sum((1:n)' .* x) / (n * sum(x)) - (n + 1) / n;
end