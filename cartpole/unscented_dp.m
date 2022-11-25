function [x, u] = unscented_dp(x0, unom, rk_f, rk_b, l, lT, T, epsilon)
%unscented_dp(x0, u, rk_f, epsilon)

    iters = 20; % will be replaced by abs(delta*J) < epsilon
    xdim = size(x0, 1);
    udim = size(unom, 1);
    % initial trajectory
    xnom = zeros(xdim, T+1);
    unom = zeros(udim, T);
    x = x0;
    xnom(:, 1) = x;
    for k = 1:T
        u = unom(k);
        x = rk_f(x, u);
        xnom(:, k+1) = x;
    end
    
       
    
    for i = 1:iters
        % autodiff first
        xT = xnom(:, end);
        dl_xT = dlarray(xT);
        [lT_val, LTx, LTxx] = dlfeval(@(x) cost_lT_diff(lT, x), dl_xT);
        
        % cost-to-go approximation
        Vx = zeros(xdim, T+1);
        Vxx = zeros(xdim, xdim, T+1);
        Vx(:, end) = LTx;
        Vxx(:, :, end) = LTxx;
        
        ds     = zeros(udim, T);
        Ks     = zeros(udim, xdim, T);
        
        for k=T:-1:1
            
            % autodiff first
            x = xnom(:, k+1);
            u = unom(:, k);
            dl_x = dlarray(x);
            dl_u = dlarray(u);
            [~, Lx, Lu, Lxx, Lxu, Luu] = dlfeval(@(x, u) cost_l_diff(l, x, u), dl_x, dl_u);

            S = Vxx(:, :, k+1);
            v = Vx(:, k+1);
            [Q,Qx,Qu,Qxx,Qxu,Quu] = udp_sigma_qs(S,v,rk_b,Lx,Lu,Lxx,Lxu,Luu,x,u);
            [K,d,S,v,V0] = ddp_control_and_value(Q,Qx,Qu,Qxx,Qxu,Quu);
            Vxx(:, :, k) = S;
            Vx(:, k) = v;
            % save controls/gains
            ds(:,k)      = d;
            Ks(:,:,k)    = K;
 
            
        end
        [unom,xnom] = udp_forward_pass(x0,rk_f,Ks,ds, xnom, unom)
        
%         [K, l, dV] = udp_backward_pass();
%         [x, u, dJ] = udp_forward_pass(x, u, K, l, dV);
        
    end
end

