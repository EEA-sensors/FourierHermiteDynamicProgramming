function [ W, XI, u ] = ut3_ws(n,kappa)
% [ W, XI, u ] = ut3_ws(n,[kappa])
%
% Return weights and sigma-points for 3rd order
% UT for dimension n with parameter kappa (default 1-n).

    if nargin < 2
        kappa = 1-n;
    end
        
    % Weights
    W = zeros(1,2*n+1);
    for j=1:2*n+1
        if j==1
            wm = kappa / (n + kappa);
        else
            wm = 1 / (2 * (n + kappa));
        end
        W(j) = wm;
    end

    % Sigma points
    XI = [zeros(n,1) eye(n) -eye(n)];
    XI = sqrt(n + kappa)*XI;
    u  = sqrt(n + kappa);
end

