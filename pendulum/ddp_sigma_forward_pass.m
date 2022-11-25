function [us,xs] = ddp_sigma_forward_pass(x0,f,Ks,ds, xnom, unom)
% Syntax: 
%    [us,xs] = ddp_sigma_forward_pass(x0,f,Ks,ds, xnom, unom)
%
% In:
%   x0 - Initial state xdim * 1
%   Ks - Control gains udim * xdim * T for k=0,...,T-1
%   ds - Control biases udim * T for k=0,...,T-1
%   f  - function signature f(x, u)
%   xnom - nominal state trajectory
%   unom - nominal control trajectory
% Out:
%   us - Optimal controls for k=0,...,T-1 as vectors udim * T
%   xs - Optimal states for k=0,...,T as vectors xdim * (T+1)
%
% Description:
%   Compute DDP forward pass corresponding to ddp_sigma_backward_pass.

    T    = size(Ks,3);
    xdim = size(Ks,2);
    udim = size(Ks,1);
        
    us = zeros(udim,T);
    xs = zeros(xdim,T+1);
    xs(:,1) = x0;
    
    x = x0;
    for k=1:T
        dx = x - xnom(:, k);
        du = -Ks(:,:,k) * dx + ds(:,k);
        u = du + unom(:, k);
        if iscell(f)
            x = f{k}(x, u);
        else
            x = f(x, u);
        end
        us(:,k) = u;
        xs(:,k+1) = x;        
    end
end

