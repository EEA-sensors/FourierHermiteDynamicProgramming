% Syntax: 
%    [u,x,P] = ddp_taylor_forward_step(x, u, K,d,dx, f,Fx,Fu,Fxx,Fxu,Fuu,G,Q,P)
%
% In:
%   x - State xdim * 1
%   u - Nominal control udim * 1
%   K - Control gain udim * xdim
%   d - Control bias udim * 1
%   dx - difference to nominal trajectory
%   f,Fx,Fu,Fxx,Fxu,Fuu - f(x,u) and its derivatives (functions or values).
%                         Note that Fx and Fu are Jacobians and the second
%                         derivatives are given as cell arrays over the
%                         output dimensions of f (they can also be empty).
%   G - Noise gain xdim * wdim       (default: identity)
%   Q - Noise covariance wdim * wdim (default: zeros)
%   P - State covariance xdim * xdim (default: zeros)
%
% Out:
%   u - Optimal control
%   x - Optimal states
%   P - State covariances
%
% Description:
%   Compute Taylor-series based DDP forward step. If parameters from Fx
%   are not given, performs a simple deterministic step. 

function [u,x,P] = ddp_taylor_forward_step(x,u, K,d, dx, f,Fx,Fu,Fxx,Fxu,Fuu,G,Q,P)
    if nargin < 7
        % Deterministic
        P = zeros(size(x,1));
%         du = - K * dx + d;
%         u = u + du;
%         x = f(x,u);  
        u = -K * x + d;
        x = f(x,u); 
        
    else
        error('Not implemented');
        % Stochastic case is implemented as simple linearization
        % f(tx,tu) = f0 + Fx (tx - x) + Fu (tu - u)
        %          = f0 + Fx (tx - x) - Fu K (tx - x)
        %          = f0 + (Fx - Fu K) (tx - x)
        % TODO: could replace with second order version of this
        x = f(x,u);
        F = Fx(x,u) - Fu(x,u) * K;
        P = F * P * F' + G * Q * G';
    end    
end
