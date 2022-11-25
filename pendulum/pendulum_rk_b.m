function [x] = pendulum_rk_b(x, u, dt, g, pl)
% Syntax:
%   [x] = pendulum_rk_b(x, u, dt)
%
% In:
%   x - current state xdim * 1
%   u - current control udim * 1
%   dt - time step
%   g - gravitational force 9.81
%   pl - pole length
% Out:
%   x - next state
%
% Description:
%   Solve the ODE with Backward Runge--Kutta method given x
% by SSH'21
    
    steps = 1;
 
    for i = 1:steps
        dx1 = pendulum_f2(x, u, g, pl)*dt; 
        dx2 = pendulum_f2(x - dx1*0.5, u, g, pl)*dt;
        dx3 = pendulum_f2(x - dx2*0.5, u, g, pl)*dt;
        dx4 = pendulum_f2(x - dx3, u, g, pl)*dt;
        x = x - (1/6)*(dx1 + 2*dx2 + 2*dx3 + dx4);
    end
 
end