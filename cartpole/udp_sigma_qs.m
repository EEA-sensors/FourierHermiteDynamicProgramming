function [Q, Qx, Qu, Qxx, Qxu, Quu] = udp_sigma_qs(S, v, rk_b, Lx, Lu, Lxx, Lxu, Luu, x, u)
% Syntax:
%   [Q, Qx, Qu, Qxx, Qxu, Quu] = udp_sigma_qs(S, v, rk_b, Lx, Lu, Lxx, Lxu, Luu, x, u)
%
% In:
%   S - Matrix Value function parameter
%   v - vector Value function parameter
%   x,u - Nominal state and control
%   rk_b - 4th order runge-kutta backward dynamic function of signature: rk_b(x,u) 
%   Lx, Lu, Lxx, Lxu, Luu - derivatives of the cost function
%
% Out:
%  Q, Qx, Qu, Qxx, Qxu, Quu - Qs coefficient
% 
% Description:
%   Unscented differential dynamic programming backward step.
%   Adapted Reference implementation from the paper: 
%   Z. Manchester and S. Kuindersma, 
%   "Derivative-free trajectory optimization with unscented dynamic programming," 
%   2016 IEEE 55th Conference on Decision and Control (CDC), 2016, 
%   pp. 3642-3647, doi: 10.1109/CDC.2016.7798817.
%   https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7798817

% by SSH'21

    b = 0.01; %2.8; % beta: The value was set in the original implementation beta = 2.8 for pendulum and beta = 0.01 for cartpole
    H = S;
    R = Luu;
    g = v;
    P = blkdiag(H^-1, R^-1);
    mu = [x; u];
    
    %Generate sigma points from Vxx(i+1) and cuu(i+1)
    L = chol(P, 'lower'); % Equation (24)
    
    n = size(x, 1);
    m = size(u, 1);
    Sig = zeros(n+m, 2*(n+m));
    for i = 1:(n+m)
        Sig(:, i) = mu + b * L(:, i);
        Sig(:, i + n + m) = mu - b * L(:, i);
    end
    
    % Equation (30)
    gs = zeros(2*(n+m), 1);

    for i = 1:(n+m)
        gs(i) = g'*b*L(1:n, i);
        gs(i+n+m) = g'*-b*L(1:n, i);
    end

    % backpropagate sample points
    Sb = Sig;
    for i = 1:2*(n+m)
        y = rk_b(Sig(1:n, i), Sig(n+1:end, i));
        Sb(1:n, i) = y;
    end
    
    
    % Equation (32): compute D matrix
    D = zeros(n+m, n+m);
    for i=1:n+m
        D(1:n, i) = Sb(1:n, i) - Sb(1:n, i+n+m);
        D(n+1:end, i) = Sb(n+1:end, i) - Sb(n+1:end, i+n+m);
    end
    D = D';

    % Equation (33): compute d vector
    d = zeros(n+m, 1);
    for i=1:n+m
        d(i) = gs(i) - gs(i+n+m);

    end
    
    % Solve the linear system in Equation (31) to retrive a-, b-
    QxQu = -D\d;
                
    % Compute mu and M for Hessian
    mu_prime = zeros(n+m, 1);
    for i = 1:2*(n+m)
        mu_prime = mu_prime + (1/(2*(m+n)))*Sb(:, i);
    end

    M = zeros(n+m);
    for i = 1:2*(n+m)
        M = M + (1/(2*b*b))*(Sb(:, i)-mu_prime)*(Sb(:, i)-mu_prime)';
    end

    %Qs = M^-1 + blkdiag(S, 0); % S = W Eq. (28)
    Qs = M^-1; 
    Qx = QxQu(1:n) + Lx;
    Qu = QxQu(n+1:end) + Lu;
             
    Qxx = Qs(1:n, 1:n) + Lxx;
    Qxu = Qs(n+1:end, 1:n)';
    Quu = Qs(n+1:end, n+1:end);
    Q = 0;
    
    
%     Qx = QxQu(1:n) + S*x + v; % Eq. (36) equivalent to Eq. (7)
%     Qu = QxQu(n+1:end) + R*u + r; % Eq. (37) equivalent to Eq. (8)
    

end