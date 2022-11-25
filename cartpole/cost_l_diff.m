function [y, Lx, Lu, Lxx, Lxu, Luu] = cost_l_diff(l, x, u)
%% Syntax:
%   [y, Lx, Lu] = cost_l_diff(l, x, u)
%   [y, Lx, Lu, Lxx, Lxu, Luu] = cost_l_diff(l, x, u)
%
% In:
%   l - quadratic cost function signature l(x, u)
%   x - state
%   u - control
%
% Out:
%   y - intermediate cost
%   Lx, Lu - first order derivatives w.r.t x, u
%   Lxx, Luu, Lxu - second order derivatives w.r.t x, u
%
% Description:
%   [y, Lx, Lu] = cost_l_diff(l, x, u) returns
%       Intermediate cost function and its first order derivatives.
%
%   [y, Lx, Lu, Lxx, Lxu, Luu] = cost_l_diff(l, x, u) returns
%       Intermediate cost function and its first and second order derivatives.
%
% by SSH'21
    y = l(x, u);
    
    Lx = [];
    Lu = [];
    Lxx = [];
    Lxu = [];
    Luu = [];
    xdim = size(x, 1);
    udim = size(u, 1);
   
    if strcmp(class(x), 'dlarray')
        [Lx, Lu] = dlgradient(y, x, u, 'EnableHigherDerivatives', true);
        
        if nargout > 3
            Lxx = dlarray(zeros(xdim, xdim));
            for i = 1:xdim
                [Lxx(:, i), ~] = dlgradient(Lx(i), x, u);
            end        

            Lxu = dlarray(zeros(xdim, udim));
            Luu = dlarray(zeros(udim, udim));
            for i = 1:udim
                [Lxu(:, i), Luu(:, i)] = dlgradient(Lu(i), x, u);
            end
            % conversion of 2nd order derivatives
            Lxx = extractdata(Lxx);
            Luu = extractdata(Luu);
            Lxu = extractdata(Lxu);
        end
                
        % conversion of y and 1st order derivatives
        y = extractdata(y);
        Lx = extractdata(Lx);
        Lu = extractdata(Lu);
        
        
    end
end