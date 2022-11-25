function [lT] = pendulum_lT(x, xg, QN)
%% Syntax:
%   [lT] = pendulum_lT(x, u, xg QN)
%
% In:
%   x - state at time t xdim*1
%   QN - state transition matrix xdim*xdim
%   xg - target state
%
% Out:
%   lT - terminal cost
%
% Description:
%   Terminal cost function of the pendulum dynamics.
% by SSH'21
    lT = 0.5*(x-xg)'*QN*(x-xg);    
end