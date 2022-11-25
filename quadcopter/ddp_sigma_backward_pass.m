function [Ks, ds, Ss, vs] = ddp_sigma_backward_pass(f,l,lT,LP,LPT,p,xnom,unom, lambda, ws)
% Syntax:  
%    [Ks, ds, Ss, vs] = ddp_sigma_backward_pass(f,l,lT,LP,LPT,p,xnom,unom, lambda)
% 
% In:
%   f      - dynamic model function signature f(x, u) or a cell of functions
%   l      - intermediate quadratic cost function, signature l(x,u) or a cell
%   of functions
%   lT     - terminal quadratic cost function, signature lT(x)
%   LP     - joint covariance matrix or a cell of matrices of nominal trajectory and control
%   LPT    - convariance matrix of nominal trajectory
%   p      - order of Hermite polynomials if Gauss--Hermite gh_ws() is used.
%   xnom   - nominal state trajectory xdim * T
%   unom   - nominal control trajectory udim * T-1
%   lambda - optional parameter
%   ws     - function to generate weights W and sigma-points SX
%            signature: 
%               gh_ws(xdim, p) - Generate Gauss-Hermite cubature order p
%               based W and SX
%               utX_ws(xdim) - Generate W and SX for X-order UT
%
% Out:
%   Ks - Control gains udim * xdim * T for k=0,...,T-1
%   ds - Control biases udim * T for k=0,...,T-1
%   Ss - Value function matrices xdim * xdim * (T+1) for k=0,...,T
%   vs - Value function vectors xdim * (T+1) for k=0,...,T
%
% Description:
%    Perform backward pass of sigma-point DDP for generic model. For example,
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
    if nargin < 9
        lambda = 0;
    end
    if nargin < 10
        ws = @(xdim) gh_ws(xdim, p);
    end
    
    T = size(xnom, 2);
    xdim = size(xnom,1);
    udim = size(unom,1);
           
    Ks = zeros(udim,xdim,T-1);
    ds = zeros(udim,T-1);
    Ss = zeros(xdim,xdim,T);
    vs = zeros(xdim,T);
           
    xT = xnom(:, end);
    [WT,XIT] = ws(xdim);
    [S, v, V0] = ddp_sigma_backward_init(xT, LPT, lT, WT, XIT); 
    Ss(:,:,end) = S;
    vs(:,end) = v;
    
    [Wk,XIk] = ws(xdim+udim);
           
    for k = T-1:-1:1
        
        u = unom(:, k);
        x = xnom(:, k);
        xn = xnom(:, k+1);
        
        if iscell(f) && iscell(LP)
%             fk = @(x, u) f{k}(x, u, k);
%             lk = @(x, u) l{k}(x, u, k);
            [K, d, S, v] = ddp_sigma_backward_step(S,v,xn,x,u,LP{k},l{k},f{k},Wk,XIk,lambda);
        elseif iscell(f) && ~iscell(LP)
            [K, d, S, v] = ddp_sigma_backward_step(S,v,xn,x,u,LP,l{k},f{k},Wk,XIk, lambda);
        elseif ~iscell(f) && iscell(LP)
            [K, d, S, v] = ddp_sigma_backward_step(S,v,xn,x,u,LP{k},l,f,Wk,XIk, lambda);
        else
            [K, d, S, v] = ddp_sigma_backward_step(S,v,xn,x,u,LP,l,f,Wk,XIk, lambda);
        end
        Ks(:,:,k) = K;
        ds(:,k)   = d;
        Ss(:,:,k) = S;
        vs(:,k)   = v;
    end
end