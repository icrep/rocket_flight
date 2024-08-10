% Main script
clear; clc;

% Gather variable constants
C_v = 5500 ;    % Rocket's exhaust velocity (m/s)
C_d = 0.5 ;     % Aerodynamic drag coefficient
A_x = 1 ;       % Cross-sectional area of the rocket (m^2)
Q = 0.05;       % Rate of fuel cnsumption (s^-1)
m_0 = 2500;    % Mass of rocket (kg)
T_f = input('Enter total time of simulation (s): \n');    % Total duration of the simulation.
n = input('Enter number of time steps: \n');  % Time step size for the simulation.

[t, Vs, Hs] = performSim(C_d, C_v, A_x, Q, m_0, T_f, n);

figure;
% Plot for Altitude over time
subplot(2, 1, 1); 
plot(t, Hs, '--r');
title('Rocket Altitude Over Time');
xlabel('Time (s)');
ylabel('Altitude (m)');
grid on; 

% Plot for velocity
subplot(2, 1, 2); 
plot(t, Vs, '--b');
title('Rocket Velocity Over Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
grid on; 