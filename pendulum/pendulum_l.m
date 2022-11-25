function [l] = pendulum_l(x, u, xg, Q, R)
%% Syntax:
%   [l] = pendulum_l(x, u, xg, Q, R)
%
% In:
%   x - state at time t xdim*1
%   u - control at time t udim*1
%   Q - state transition matrix xdim*xdim
%   R - control matrix udim*udim
%   xg - target state
%
% Out:
%   l - intermediate cost
%
% Description:
%   Intermediate cost function of the pendulum dynamics.
% by SSH'21
    l = 0.5*(x-xg)'*Q*(x-xg)+0.5* u'*R*u;    
end