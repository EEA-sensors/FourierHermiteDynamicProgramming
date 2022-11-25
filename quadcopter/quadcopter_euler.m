function [x] = quadcopter_euler(x, u, Ixx, Iyy, Izz, omega_r, Jr, m, l, g, dt)
% Syntax: 
%    [x] = quadcopter_euler(x, u, Ixx, Iyy, Izz, omega_r, Jr, m, l, g, dt)
%
% In:
%   x - state xdim * 1
%   u - control udim * 1
%   dt - time step
%   g - gravitation force ~9.81
%   pl - pole length
% Out:
%   y - next state 
%
% Description:
%   Solve the ODE with Euler method given x
    steps = 1;
    for i  = 1:steps
         x = x + quadcopter_f(x, u, Ixx, Iyy, Izz, omega_r, Jr, m, l, g)*dt;
    end
end