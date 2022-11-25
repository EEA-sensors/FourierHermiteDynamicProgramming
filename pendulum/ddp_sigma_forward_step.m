% Syntax: 
%    [u,x,P] = ddp_sigma_forward_step(x,K,d,f,Fx,Fu,Fxx,Fxu,Fuu,G,Q,P)
%
% In:
%   x - State xdim * 1
%   K - Control gain udim * xdim
%   d - Control bias udim * 1
%   f - Function f(x,u) .
%   W,XI - Sigma-point weights and unit points jointly for x
%   G - Noise gain xdim * wdim       (default: identity)
%   Q - Noise covariance wdim * wdim (default: zeros)
%   P - State covariance xdim * xdim (default: zeros)
%
% Out:
%   u - Optimal control
%   x - Optimal state
%   P - State covariance
%
% Description:
%   Compute sigma-point based DDP forward step. If parameters from Fx
%   are not given, performs a simple deterministic step. 

function [u,x,P] = ddp_sigma_forward_step(x,K,d,f,W,XI,G,Q,P)
    if nargin < 5
        % Deterministic
        P = zeros(size(x,1));
        u = -K * x + d;
        x = f(x,u);        
    else
        % Stochastic
        if nargin < 9 || all(P(:) == 0)
            u = -K * x + d;
            x = f(x,u);        
            P = G * Q * G';
        else
            % Compute prediction is mean and covariance of x |-> f(x,-K x + d)
            % by using numerical integration. Returned u is computed using
            % the nominal x though.
            u = -K * x + d;
            LP = chol(P,'lower');
            new_x = zeros(size(x));
            tf = zeros(size(x,1),length(W));
            for i=1:length(W)
                tx = LP * XI(:,i) + x;
                tf(:,i) = f(tx, -K * tx + d);
                new_x = new_x + W(i) * tf(:,i);
            end
            new_P = G * Q * G';
            for i=1:length(W)
                new_P = new_P + W(i) * (tf(:,i) - new_x) * (tf(:,i) - new_x)';
            end
            x = new_x;
            P = new_P;
        end
    end    
end

