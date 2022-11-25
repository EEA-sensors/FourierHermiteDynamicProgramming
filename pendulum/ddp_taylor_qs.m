% Syntax:
%   [Q0,Qx,Qu,Qxx,Qxu,Quu] = ddp_taylor_qs(x,u,l,Lx,Lu,Lxx,Lxu,Luu,f,Fx,Fu,Fxx,Fxu,Fuu,S,v,V0,xn)
%
% In:
%   x,u - Nominal state and control
%   l,Lx,Lu,Lxx,Lxu,Luu - l(x,u) and its derivatives (functions or values).
%                         Note that Lx and Lu are gradients i.e. column
%                         vectors.
%   f,Fx,Fu,Fxx,Fxu,Fuu - f(x,u) and its derivatives (functions or values).
%                         Note that Fx and Fu are Jacobians and the second
%                         derivatives are given as cell arrays over the
%                         output dimensions of f (they can also be empty).
%   S,v,V0,xn - Value function parameters and its linearization point.
%
% Out:
%   Q0,Qx,Qu,Qxx,Qxu,Quu - Parameters of function Q(x,u)
%
% Description:
%   Use Taylor-series-based differential dynamic programming to approximate
%   Q(x,u) = l(x,u) + V(f(x,u)) in form
%
%    Q(x + dx,u + du) = Q0 + Qx^T dx + Qu^T du + 1/2 [dx;du]^T [Qxx Qxu; Qux Quu] [dx;du]
%
%   The approximation is based on second order Taylor series expansion
%   centered at x,u. The value function at next step is given as
%
%    V(xn + dx, un + du) = 1/2 dx^T S dx - v^T dx + V0

% by SS'21

function [Q0,Qx,Qu,Qxx,Qxu,Quu] = ddp_taylor_qs(x,u,l,Lx,Lu,Lxx,Lxu,Luu,f,Fx,Fu,Fxx,Fxu,Fuu,S,v,V0,xn)
    if isa(l,'function_handle')
        l_val   = l(x,u);
        Lx_val  = Lx(x,u);
        Lu_val  = Lu(x,u);
        Lxx_val = Lxx(x,u);
        Lxu_val = Lxu(x,u);
        Luu_val = Luu(x,u);
    else
        l_val   = l;
        Lx_val  = Lx;
        Lu_val  = Lu;
        Lxx_val = Lxx;
        Lxu_val = Lxu;
        Luu_val = Luu;
    end
    
    if isa(f,'function_handle')
        f_val   = f(x,u);
        Fx_val  = Fx(x,u);
        Fu_val  = Fu(x,u);

        Fxx_val = cell(size(Fxx));
        Fxu_val = cell(size(Fxu));
        Fuu_val = cell(size(Fuu));

        for m=1:length(Fxx)
            Fxx_val{m} = Fxx{m}(x,u);
            Fxu_val{m} = Fxu{m}(x,u);
            Fuu_val{m} = Fuu{m}(x,u);
        end
    else
        f_val   = f;
        Fx_val  = Fx;
        Fu_val  = Fu;

        Fxx_val = cell(size(Fxx));
        Fxu_val = cell(size(Fxu));
        Fuu_val = cell(size(Fuu));

        for m=1:length(Fxx)
            Fxx_val{m} = Fxx{m};
            Fxu_val{m} = Fxu{m};
            Fuu_val{m} = Fuu{m};
        end
    end
    
    df  = f_val - xn;
    Vx  = S * df - v;
    Vxx = S;
    
    
    Q0  = l_val + 0.5 * df' * S * df - v' * df + V0;
    Qx  = Lx_val + Fx_val' * Vx;
    Qu  = Lu_val + Fu_val' * Vx;
    Qxx = Lxx_val + Fx_val' * Vxx * Fx_val;
    for m=1:length(Fxx_val)
        Qxx = Qxx + Fxx_val{m} * Vx(m);
    end
    Qxu = Lxu_val + Fx_val' * Vxx * Fu_val;
    for m=1:length(Fxu_val)
        Qxu = Qxu + Fxu_val{m} * Vx(m);
    end
    Quu = Luu_val + Fu_val' * Vxx * Fu_val;
    for m=1:length(Fuu_val)
        Quu = Quu + Fuu_val{m} * Vx(m);
    end
    
    
    
end
