% Syntax:
%   [Q0,Qx,Qu,Qxx,Qxu,Quu] = ddp_sigma_qs(x,u,LP,l,f,W,XI,S,v,V0,xn)
%   
% In:
%   x,u - Nominal state and control
%   LP  - Matrix such that P = LP * LP' for joint covariance P of x,u
%   l,f - Functions l(x,u) and f(x,u)
%   W,XI - Weight and unit sigma points jointly for x,u
%   S,v,V0,xn - Value function parameters and its linearization point
%   
% Out:
%   Q0,Qx,Qu,Qxx,Qxu,Quu - Parameters of function Q(x,u)
%
% Description:
%   Use sigma-point differential dynamic programming to approximate
%   Q(x,u) = l(x,u) + V(f(x,u)) in form
%
%    Q(x + dx,u + du) = Q0 + Qx^T dx + Qu^T du + 1/2 [dx;du]^T [Qxx Qxu; Qux Quu] [dx;du]
%
%   The approximation is based on second order Fourier-Hermite series expansion
%   centered at x,u with covariance P = LP LP^T. The value function at next step
%   is given as
%
%    V(xn + dx, un + du) = 1/2 dx^T S dx - v^T dx + V0

% by SS'21

function [Q0,Qx,Qu,Qxx,Qxu,Quu] = ddp_sigma_qs(x,u,LP,l,f,W,XI,S,v,V0,xn)

    xdim = size(x,1);
    udim = size(u,1);

    aQ = 0;
    bQ = zeros(xdim+udim,1);
    CQ = zeros(xdim+udim,xdim+udim);
    
    for i=1:length(W)
        xu = LP * XI(:,i) + [x;u];
        tx = xu(1:xdim,1);
        tu = xu(xdim+1:end,1);
        tf = f(tx,tu);
        tQ = l(tx,tu) + 0.5 * (tf-xn)' * S * (tf-xn) - v' * (tf-xn) + V0;
        aQ = aQ + W(i) * tQ;
        bQ = bQ + W(i) * tQ * XI(:,i);
        CQ = CQ + W(i) * tQ * (XI(:,i) * XI(:,i)' - eye(xdim+udim));
    end
    
    Q0 = aQ - 0.5 * trace(CQ);
    
    tmp = LP' \ bQ;
    Qx = tmp(1:xdim,1);
    Qu = tmp(xdim+1:end,1);
    
    tmp = LP' \ CQ / LP;
    Qxx = tmp(1:xdim,1:xdim);
    Qxu = tmp(1:xdim,xdim+1:end);
    Quu = tmp(xdim+1:end,xdim+1:end);  
    
end
