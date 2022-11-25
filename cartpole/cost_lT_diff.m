function [lT_val, LTx, LTxx] = cost_lT_diff(lT, x)
%% Syntax:
%   [lT_val, LTx] = cost_lT_diff(lT, x)
%   [lT_val, LTx, LTxx] = cost_lT_diff(lT, x)
%
% In:
%   lT - quadratic cost function lT(xT)
%   x - Terminal state
%
% Out:
%   lT_val - Terminal cost
%   LTx, LTxx - derivatives
%
% Description:
%   [lT_val, LTx] = cost_lT_diff(lT, x) returns
%       Terminal cost function and its first order derivatives.
%   [lT_val, LTx, LTxx] = cost_lT_diff(lT, x) returns
%       Terminal cost function and its first and second order derivatives.
% by SSH'21
    lT_val = lT(x);
    % derivatives
    LTx = [];
    LTxx = [];
    xdim = size(x, 1);
           
    if strcmp(class(x), 'dlarray')
        LTx = dlgradient(lT_val, x, 'EnableHigherDerivatives', true);
        
        if nargout > 2
            % Lxx
            for i = 1:xdim
                [Lxx] = dlgradient(LTx(i), x); 
                LTxx = [LTxx; Lxx'];
            end
            LTxx = extractdata(LTxx);
        end
        
        
        % conversion
        lT_val = extractdata(lT_val);
        LTx = extractdata(LTx);
        
          
    end
end