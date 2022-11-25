function [x] = cartpole_rk(x, u, dt, mp, mc, pl, g)
% Syntax:
%   [x] = cartpole_rk(x, u, dt,mp, mc, pl, g)
%
%   x - state xdim * 1
%   u - control udim * 1
%   dt - time step
%   mp - mass of pole in kg
%   mc - mass of cart in kg
%   g - gravitation force ~9.81
%   pl - pole length
% Out:
%   x - next state 
%
% Description:
%   Solve the ODE with Runge--Kutta method given x
% by SSH'21
    
    steps = 1;
    
    for i = 1:steps
        dx1 = cartpole_f(x, u, mp, mc, pl, g)*dt;
        dx2 = cartpole_f(x + dx1*0.5, u, mp, mc, pl, g)*dt;
        dx3 = cartpole_f(x + dx2*0.5, u, mp, mc, pl, g)*dt;
        dx4 = cartpole_f(x + dx3, u, mp, mc, pl, g)*dt;
        x = x + (1/6)*(dx1 + 2*dx2 + 2*dx3 + dx4);
    end
end