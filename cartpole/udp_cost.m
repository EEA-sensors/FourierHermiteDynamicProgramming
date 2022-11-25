function [J] = udp_cost(x, u, WN, Wk, Rk, wN, wk, rk)
% Syntax:
%   [J] = udp_cost(x, u, WN, Wk, Rk, wN, wk, rk)
%
% In:
%   x - state trajectory
%   u - input/control trajectory
%   WN, Wk, Rk - cost weighting matrices
%   wN, wk, rk - cost weighting vectors
%
% Out:
%   J - non-quadratic cost function
%
% Description:
%   cost function
% by SSH'21
    N = size(u, 2);
    J = 0.5*x(:, end)'*WN*x(:, end) + wN' * x(:, end);
    for k = 1:N-1
        J = J + 0.5*x(:, k)'*Wk*x(:, k) + wk'*x(:, k) ...
              + 0.5*u(:, k)'*Rk*u(:, k) + rk'*u(:, k);
    end
end
