% Syntax:
%   v = gaussian_moments(monomial_exps)
%
% Description:
%   Compute int x^{moment_exps} N(x | 0,I) dx

function v = gaussian_moments(monomial_exps)
    ind = find(rem(monomial_exps,2) == 0);
    if length(ind) ~= length(monomial_exps)
        v = 0; % At least one was odd
    else
        v = prod(double_factorial(monomial_exps-1));
    end
end

