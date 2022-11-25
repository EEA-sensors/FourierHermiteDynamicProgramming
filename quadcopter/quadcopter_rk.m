function [x] = quadcopter_rk(x, u, Ixx, Iyy, Izz, omega_r, Jr, m, l, g, dt)
% Syntax:
%   [x] = quadcopter_rk(x, u, Ixx, Iyy, Izz, omega_r, Jr, m, l, g, dt)
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
%   Solve the ODE with Runge--Kutta method given x
% by SSH'22
    
    steps = 1;
    
    for i = 1:steps
        dx1 = quadcopter_f(x, u, Ixx, Iyy, Izz, omega_r, Jr, m, l, g)*dt; 
        dx2 = quadcopter_f(x + dx1*0.5, u, Ixx, Iyy, Izz, omega_r, Jr, m, l, g)*dt;
        dx3 = quadcopter_f(x + dx2*0.5, u, Ixx, Iyy, Izz, omega_r, Jr, m, l, g)*dt;
        dx4 = quadcopter_f(x + dx3, u, Ixx, Iyy, Izz, omega_r, Jr, m, l, g)*dt;
        x = x + (1/6)*(dx1 + 2*dx2 + 2*dx3 + dx4);
    end
end