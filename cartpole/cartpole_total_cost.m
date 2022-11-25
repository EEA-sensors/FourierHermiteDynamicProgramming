function [J] = cartpole_total_cost(xs, us, xg, Q, QN, R)
%% Syntax:
%   [J] = cartpole_total_cost(xs, us, Q, QN, R)
%
% In:
%   xs - states from timestep 1, ..., T+1 xdim*T+1
%   xg - target state
%   us - controls from timestep 1, ..., T xdim*T
%   Q - state transition matrix xdim*xdim
%   R - control matrix udim*udim
%   QN - terminal matrix xdim*xdim
%   xg - target state
%
% Out:
%   J - total cost
%
% Description:
%   total costs  of the cartpole dynamics.
%   J = 0.5*x_T'*Q*x_T + \sum_i (0.5*x_t'*Q*x_t + 0.5*u_t'*R*u_t)
% by SSH'21

    T = size(us, 2);
    J = 0;
    for t = 1:T
        J = J+0.5*(xs(:,t)-xg)'*Q*(xs(:,t)-xg)+0.5*us(t)'*R*us(t);
    end
    J = J+0.5*(xs(:, end)-xg)'*QN*(xs(:, end)-xg);
end