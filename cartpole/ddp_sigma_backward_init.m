% Syntax:
%   [S, v, V0] = ddp_sigma_backward_init(x,lT,W,XI)
%
% In:
%   x    - Nominal state
%   LP   - Matrix such that P = LP * LP' for joint covariance P of x
%   lT   - Terminal cost function lT(x) and its derivatives
%   W,XI - Sigma-point weights and unit points for x
%
% Description:
%   Sigma-point based differential dynamic programming initialization step.

% by SS'21

function [S, v, V0] = ddp_sigma_backward_init(x,LP,lT,W,XI)
    xdim = size(x,1);
    
    al = 0;
    bl = zeros(xdim,1);
    Cl = zeros(xdim,xdim);
    
    for i=1:length(W)
        tx = LP * XI(:,i) + x;
        tl = lT(tx);
        al = al + W(i) * tl;
        bl = bl + W(i) * tl * XI(:,i);
        Cl = Cl + W(i) * tl * (XI(:,i) * XI(:,i)' - eye(xdim));
    end
    
    V0 = al - 0.5 * trace(Cl);
    v  = -LP' \ bl;
    S  = LP' \ Cl / LP;
end

