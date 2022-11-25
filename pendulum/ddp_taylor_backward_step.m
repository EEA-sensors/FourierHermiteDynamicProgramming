% Syntax:
%   [K, d, S, v] = ddp_taylor_backward_step(S,v,xn,x,u,l,Lx,Lu,Lxx,Lxu,Luu,f,Fx,Fu,Fxx,Fxu,Fuu, lambda)
%
% In:
%   S,v,xn - Value function parameters and its linearization point
%   x,u    - Nominal state and control
%   l,Lx,Lu,Lxx,Lxu,Luu - l(x,u) and its derivatives (functions or values).
%                         Note that Lx and Lu are gradients i.e. column
%                         vectors.
%   f,Fx,Fu,Fxx,Fxu,Fuu - f(x,u) and its derivatives (functions or values).
%                         Note that Fx and Fu are Jacobians and the second
%                         derivatives are given as cell arrays over the
%                         output dimensions of f (they can also be empty).
%   lambda - optional parameter
%
% Description:
%   Taylor-series based differential dynamic programming backward step.

% by SS'21

function [K, d, S, v, V0] = ddp_taylor_backward_step(S,v,V0,xn,x,u,l,Lx,Lu,Lxx,Lxu,Luu,f,Fx,Fu,Fxx,Fxu,Fuu,lambda)
%     V0 = 0; % Does not matter: Commented by Sakira: probably matters,
%     however does not affect the control. 
    if nargin < 19
        lambda = 0;
    end
    
    [Q0,Qx,Qu,Qxx,Qxu,Quu] = ddp_taylor_qs(x,u,l,Lx,Lu,Lxx,Lxu,Luu,f,Fx,Fu,Fxx,Fxu,Fuu,S,v,V0,xn);
    [K,d,S,v, V0] = ddp_control_and_value(Q0,Qx,Qu,Qxx,Qxu,Quu,lambda);
end

