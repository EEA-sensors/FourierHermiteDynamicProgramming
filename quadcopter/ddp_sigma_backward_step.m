% Syntax:
%   [K, d, S, v] = ddp_sigma_backward_step(S,v,xn,x,u,LP,l,f,W,XI, lambda)
%
% In:
%   S,v,xn - Value function parameters and its linearization point
%   x,u    - Nominal state and control
%   LP     - Matrix such that P = LP * LP' for joint covariance P of x,u
%   l,f    - Functions l(x,u) and f(x,u)
%   W,XI   - Sigma-point weights and unit points jointly for x,u
%   lambda - Optional parameter
%
% Description:
%   Sigma-point based differential dynamic programming backward step.

function [K, d, S, v] = ddp_sigma_backward_step(S,v,xn,x,u,LP,l,f,W,XI, lambda)
    if nargin < 11
        lambda = 0;
    end
    V0 = 0; % Does not matter
    [Q0,Qx,Qu,Qxx,Qxu,Quu] = ddp_sigma_qs(x,u,LP,l,f,W,XI,S,v,V0,xn);
    [K,d,S,v,~] = ddp_control_and_value(Q0,Qx,Qu,Qxx,Qxu,Quu,lambda);
end
