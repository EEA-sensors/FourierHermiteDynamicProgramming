% Syntax:
%   [K,d,S,v,V0] = ddp_control_and_value(Q0,Qx,Qu,Qxx,Qxu,Quu,lambda)
%
% In:
%   Q0,Qx,Qu,Qxx,Qxu,Quu - parameters of Q(x,u)
%   lambda - optional parameter
%
% Out:
%   K,d    - parameters of control law
%   S,v,V0 - parameters of value function
%
% Description:
%   Given Q(x,u) = Q0 + Qx^T x + Qu^T u + 1/2 [x;u]^T [Qxx Qxu; Qux Quu] [x;u]
%   compute u^*(x) = arg min_u Q(x,u) and then return control law u(x) as well as
%   value function V(x) = min_u Q(x,u) which have the following forms:
%
%     u^*(x) = -K x + d
%       V(x) = 1/2 x^T S x - v^T x + V0

% by SS'21

function [K,d,S,v,V0] = ddp_control_and_value(Q0,Qx,Qu,Qxx,Qxu,Quu,lambda)
    
    % regularized term added
    if nargin < 7
        lambda = 0;
    end
    
    QuuF = Quu + lambda*eye(size(Quu,1));
    
    d  = -QuuF \ Qu;
    K  = QuuF \ Qxu';
    
    S  = Qxx + K' * Quu * K - K'*Qxu' - Qxu*K;
    v  = -Qx + K' * Quu * d + K'*Qu - Qxu*d;
    V0 = Q0 + 0.5 * d' * Qu; % Need to fix this
    
%     S  = Qxx - K' * Quu * K;
%     v  = -Qx - K' * Quu * d;
%     V0 = Q0 + 0.5 * d' * Qu;
        
end
% 
% + K_i'*Qu  + Qux'*k_i;
%         Vxx(:,:,i)  = Qxx + K_i'*Quu*K_i + K_i'*Qux + Qux'*K_i;

