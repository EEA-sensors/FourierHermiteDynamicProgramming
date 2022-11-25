% Syntax:
%   [S, v, V0] = ddp_taylor_backward_init(x,lT,LTx,LTxx)
%
% In:
%   x - Nominal state
%   lT,LTx,LTxx - Terminal cost function lT(x) and its derivatives
%
% Description:
%   Taylor-series based differential dynamic programming initialization step.

% by SS'21

function [S, v, V0] = ddp_taylor_backward_init(x,lT,LTx,LTxx)
    if isa(lT,'function_handle')
        lT_val   = lT(x);
        LTx_val  = LTx(x);
        LTxx_val = LTxx(x);
    else
        lT_val   = lT;
        LTx_val  = LTx;
        LTxx_val = LTxx;
    end
    
    S  = LTxx_val;
    v  = -LTx_val; 
    V0 = lT_val;
end
