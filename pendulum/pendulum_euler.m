function [x] = pendulum_euler(x, u, dt, g, pl)
% Syntax: 
%    [x] = pendulum_euler(x, u, dt, g, pl)
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
         x = x + pendulum_f2(x, u, g, pl)*dt;
    end
end