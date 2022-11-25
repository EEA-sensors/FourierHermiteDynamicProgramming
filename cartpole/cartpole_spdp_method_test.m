function cartpole_spdp_method_test(sf)
%% Syntax:
%   cartpole_spdp_method_test(sf)
% In:
%   sf - different sigma-point methods of signature sf(xdim, options)
%   %
% Out:
%
% Description:
%   spdp_method_test(sf) tests different spdp methods
%
% by SSH'21
    
    % parameters of the model
    g = 9.81;
    pl = 0.5;  
    mp = 1; 
    mc = 10;
    dt = 0.1;
    % rk method
    rk_f = @(x, u) cartpole_rk(x, u, dt, mp, mc, pl, g);
    
    % dimensions of state and control
    xdim = 4;
    udim = 1;
    
    % initial and target states and controls
    xg = [0; pi; 0; 0]; % target
    x0 = [0; 0; 0;  0]; % initial state
    u0 = 0;
    
    % define terminal and intermediate cost parameters and functions    
    Q = 0.1*eye(xdim); 
    R = 0.01*eye(udim); 
    QN = 1000*eye(xdim); 
    l = @(x, u) cartpole_l(x, u, xg, Q, R);
    lT = @(x) cartpole_lT(x, xg, QN);
    cost_f = @(x, u) cartpole_total_cost(x, u, xg, Q, QN, R);
    
    % other related parameters
    iters = 150;
    T = ceil(5/dt);
    lambda = 1;
    dlambda = 1;
    minLambda = 1e-6;
    maxLambda = 1e10;
    lambdaFactor = 1.6;
    alphas = 10.^linspace(0, -3, 11);
    
    % backward init for sigma
    p = 3; % order of Hermite polynomials
    c = 1e-6;
    P = c*eye(xdim + udim);
    LP = chol(P,'lower');
    PT = c*eye(xdim);
    LPT = chol(PT,'lower');
        
    % norminal trajectory and control
    xnom_sigma = zeros(xdim, T+1);
    unom_sigma = u0*ones(udim, T);
    % initialize nominal trajectory
    
    x = x0;
    xnom_sigma(:, 1) = x;
    for i=1:T
        u = unom_sigma(:, i);
        x = rk_f(x, u);
        xnom_sigma(:, i+1) = x;
    end
    
    % time
    times_sigma = zeros(iters, 1);
        
    JJ_sigma = zeros(iters, 1);
    XX_sigma = zeros(xdim, T+1, iters);
    UU_sigma = zeros(udim, T, iters);
    

    for i=1:iters
        c1 = cost_f(xnom_sigma, unom_sigma);
        tic;        
            [Ks2, ds2, Ss2, vs2] = ddp_sigma_backward_pass(rk_f, l, lT, LP, LPT, p, xnom_sigma, unom_sigma, lambda, sf);
            % line search
            fwdpassDone = 0;
            alphaOpt = 1;
            for alpha = alphas
                [us2, xs2] = ddp_sigma_forward_pass(x0, rk_f, Ks2, ds2*alpha, xnom_sigma, unom_sigma);
                c2 = cost_f(xs2, us2);
                if c2 < c1
                    fwdpassDone = 1;
                    alphaOpt = alpha;
                    break
                end
            end

        times_sigma(i) = toc;
        c2 = cost_f(xs2, us2);
        fprintf('At iter %d: Nom Cost: %f New cost: %f AlphaOpt: %f, lambda: %f\n', i, c1, c2, alphaOpt, lambda);
        
        if fwdpassDone
            xnom_sigma = xs2;
            unom_sigma = us2;
            % decrease lambda
            %dlambda = min(dlambda/lambdaFactor, 1/lambdaFactor);
            %lambda = lambda*dlambda *(lambda > minLambda);
            lambda = lambda/2;
        else
            % increase lambda
%             dlambda = max(dlambda/lambdaFactor, 1/lambdaFactor);
%             lambda = max(lambda*dlambda, minLambda);
              lambda = lambda*2;
        end
        
        
        JJ_sigma(i) = cost_f(xnom_sigma, unom_sigma);
        XX_sigma(:, :, i) = xnom_sigma;
        UU_sigma(:, :, i) = unom_sigma;
        

    end
    fprintf('End of Sigma at iter: %d\n', i);
    xs = xnom_sigma;
    us = unom_sigma;
    filename = ['results_fh_paper/cartpole_spdp_', func2str(sf), '.mat'];
    save(filename, 'JJ_sigma', 'XX_sigma', 'UU_sigma', 'times_sigma', 'xs', 'us', 'i');
end
    