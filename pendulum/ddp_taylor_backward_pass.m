function [Ks, ds, Ss, vs] = ddp_taylor_backward_pass(f,l,lT,xnom,unom, lambda)
% Syntax:  
%    [Ks, ds, Ss, vs] = ddp_taylor_backward_pass(f,l,lT,xnom,unom)
% 
% In:
%   f  - dynamic model function signature f(x, u) or a cell of functions
%   l  - intermediate quadratic cost function, signature l(x,u) or a cell
%   of functions
%   lT - terminal quadratic cost function, signature lT(x)
%   xnom - nominal state trajectory xdim * T
%   unom - nominal control trajectory udim * T-1
%   lambda - optional parameter
%
% Out:
%   Ks - Control gains udim * xdim * T for k=0,...,T-1
%   ds - Control biases udim * T for k=0,...,T-1
%   Ss - Value functionmatrices xdim * xdim * (T+1) for k=0,...,T
%   vs - Value function vectors xdim * (T+1) for k=0,...,T
%
% Description:
%    Perform backward pass of DDP for generic model. For example,
%
%      x[k+1] = F[k] x[k] + c[k] + L[k] u[k] + G[k] w[k]
%        J(u) = E{ 1/2 (H[T] x[T] - r[T)].T X[T] (H[T] x[T] - r[T])
%         + sum_{k=0}^{T-1} 1/2 (H[k] x[k] - r[k]).T X[k] (H[k] x[k] - r[k])
%                         + 1/2 (Z[k] u[k] - s[k]).T U[k] (Z[k] u[k] - s[k])
%                             + (H[k] x[k] - r[k]).T M[k] (Z[k] u[k] - s[k]) }
%
%    The value function is then given as (for k=0,...,T)
%
%      V_k(x) = 1/2 x^T S[k] x - v[k]^T x + constant
%
%    and the control law has the form (for k=0,...,T-1)
%
%      u_k(x) = -K[k] x + d[k]
%
%    The noise w[k] ~ N(0,Q[k]) which has no effect to control (but it has
%    to forward computation).
% by SSH'21

    if nargin < 6
        lambda = 0;
    end
    
    T = size(xnom, 2);
    xdim = size(xnom,1);
    udim = size(unom,1);
           
    Ks = zeros(udim,xdim,T-1);
    ds = zeros(udim,T-1);
    Ss = zeros(xdim,xdim,T);
    vs = zeros(xdim,T);
           
    xT = xnom(:, end);
    dl_xT = dlarray(xT);
    [lT_val, LTx, LTxx] = dlfeval(@(xT) cost_lT_diff(lT, xT), dl_xT);
    [S, v, V0] = ddp_taylor_backward_init(xT, lT_val, LTx, LTxx);
    
    
    Ss(:,:,end) = S;
    vs(:,end) = v;
           
    for k = T-1:-1:1
        % convert to dlarray before autodiff
        u = unom(:, k);
        x = xnom(:, k);
        xn = xnom(:, k+1);
        dl_x = dlarray(x);
        dl_u = dlarray(u);
        
        % autodiff
        if iscell(f) && iscell(l)
            [l_val, Lx, Lu, Lxx, Lxu, Luu] = dlfeval(@(x, u) cost_l_diff(l{k}, x, u), dl_x, dl_u);
            [f_val, Fx, Fu, Fxx, Fxu, Fuu] = dlfeval(@(x, u) dynamic_model_diff(f{k}, x, u), dl_x, dl_u);
        elseif iscell(f) && ~iscell(l)
            [l_val, Lx, Lu, Lxx, Lxu, Luu] = dlfeval(@(x, u) cost_l_diff(l, x, u), dl_x, dl_u);
            [f_val, Fx, Fu, Fxx, Fxu, Fuu] = dlfeval(@(x, u) dynamic_model_diff(f{k}, x, u), dl_x, dl_u);
        elseif ~iscell(f) && iscell(l)
            [l_val, Lx, Lu, Lxx, Lxu, Luu] = dlfeval(@(x, u) cost_l_diff(l{k}, x, u), dl_x, dl_u);
            [f_val, Fx, Fu, Fxx, Fxu, Fuu] = dlfeval(@(x, u) dynamic_model_diff(f, x, u), dl_x, dl_u);
             
        else
            [l_val, Lx, Lu, Lxx, Lxu, Luu] = dlfeval(@(x, u) cost_l_diff(l, x, u), dl_x, dl_u);
            [f_val, Fx, Fu, Fxx, Fxu, Fuu] = dlfeval(@(x, u) dynamic_model_diff(f, x, u), dl_x, dl_u);
            
        end
                
        [K, d, S, v, V0] = ddp_taylor_backward_step(S,v,V0,xn,x,u,l_val,Lx,Lu,Lxx,Lxu,Luu,f_val,Fx,Fu,Fxx,Fxu,Fuu,lambda);
        Ks(:,:,k) = K;
        ds(:,k)   = d;
        Ss(:,:,k) = S;
        vs(:,k)   = v;
    end
end