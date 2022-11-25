function [y] = cartpole_f(x, u, mp, mc, pl, g)
% Syntax: 
%    [y] = cartpole_f(x, u,  mp, mc, pl, g)
%
% In:
%   x - state xdim * 1
%   u - control udim * 1
%   mp - mass of pole in kg
%   mc - mass of cart in kg
%   pl - length of pole in m
%   g - gravitational force m/s^2
% Out:
%   y - next state
%
% Description:
%   The cartpole dynamic model
  
    x1 = x(1);
    theta = x(2);
    dx = x(3);
    dtheta = x(4);
    s = sin(theta);
    c = cos(theta);
      
    C = [0, -mp*pl*dtheta*s; 0, 0];
    G = [0;  mp*g*pl*s];
    B = [1; 0];
    H = [mc+mp, mp*pl*c; mp*pl*c, mp*pl^2];
    
    % MATLAB autodiff does not support inv()
    % Hinv generated from symbolic expression
    Hinv = [1/(- mp*c^2 + mc + mp),-c/(- mp*pl*c^2 + mc*pl + mp*pl);
        -c/(- mp*pl*c^2 + mc*pl + mp*pl), (mc + mp)/(- c^2*mp^2*pl^2 + mp^2*pl^2 + mc*mp*pl^2)];
      
    qdd = Hinv*(B*u - C*[dx; dtheta]- G);
    y = [dx; dtheta; qdd(1); qdd(2)];
    
end