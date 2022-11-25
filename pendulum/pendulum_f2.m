function [y] = pendulum_f2(x, u, g, l)
% Syntax: 
%    [y] = pendulum_f2(x, u)
%
% In:
%   x - state xdim * 1
%   u - control udim * 1
%   g - gravitation force
%   l - length of the pole
% Out:
%   y - next state
%
% Description:
%   The pendulum dynamic model
    b = 0.1;
    I = 0.25;
    theta = x(1);
    dtheta = x(2);
    y = [dtheta; (-g*l*sin(theta)-b*dtheta)/I] + [0; u/I];
    
end