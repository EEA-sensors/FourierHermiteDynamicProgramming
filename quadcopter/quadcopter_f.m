function [y] = quadcopter_f(x, u, Ixx, Iyy, Izz, omega_r, Jr, m, l, g)
% Syntax: 
%    [y] = quadcopter_f(x, u, Ixx, Iyy, Izz, omega_r, Jr, m, l, g)
%
% In:
%   x - state xdim * 1, states comprises: [phi, phidot, theta, thetadot, psi, psidot, z, zdot, x, xdot, y, ydot]
%   u - control udim * 1, inputs comprises: [u1, u2, u3, u4]
%   phi - roll, theta - pitch, psi - yaw
%   m - mass in kg
%   l - distance
%   Jr - rotor inertia
%   Ixx, Iyy, Izz - airframe inertia of roll, pitch and yaw
%   omega_r - angular rotation
%   g - gravitational force m/s^2
% Out:
%   y - next state
%
% Description:
%   The quadcopter dynamic model (C) 2022 SSH

    [phi, phidot, theta, thetadot, psi, psidot] = deal(x(1), x(2), x(3), x(4), x(5), x(6));

    a1 = thetadot*psidot*((Iyy-Izz)/Ixx);
    a2 = thetadot*(Jr/Ixx)*omega_r;
    a3 = (l/Ixx)*u(2);

    a4 = phidot*psidot*((Ixx-Izz)/Iyy);
    a5 = phidot*(Jr/Iyy)*omega_r;
    a6 = (l/Iyy)*u(3);
    
    a7 = thetadot*phidot*((Ixx-Iyy)/Izz);
    a8 = (1/Izz)*u(4);

    c_phi = cos(phi);
    c_psi = cos(psi);
    c_theta = cos(theta);
    s_theta = sin(theta);
    s_phi = sin(phi);
    s_psi = sin(psi);
       
    
    y1 = x(2);
    y2 = a1 + a2 + a3;
    y3 = x(4);
    y4 = a4 - a5 + a6;
    y5 = x(6);
    y6 = a7 + a8;
    y7 = x(8);
    y8 = g - (c_phi*c_theta)*(1/m)*u(1);
    y9 = x(10);
    y10 = (c_phi*s_theta*c_psi + s_phi*s_psi)*(1/m)*u(1);
    y11 = x(12);
    y12 = (c_phi*s_theta*s_psi - s_phi*c_theta)*(1/m)*u(1);


    y = [y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12]';
    
end