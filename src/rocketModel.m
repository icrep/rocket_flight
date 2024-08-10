% Main script
clear; clc;

% Gather variable constants
C_v = 5500;    % Rocket's exhaust velocity (m/s)
C_d = 0.5;     % Aerodynamic drag coefficient
A_x = 1;       % Cross-sectional area of the rocket (m^2)
Q = 0.05;      % Rate of fuel consumption (s^-1)
m_0 = 200;     % Mass of rocket (kg)
T_f = input('Enter total time of simulation (s): \n');    % Total duration of the simulation.
n = input('Enter number of time steps: \n');  % Time step size for the simulation.

[t, Vs, Hs] = performSim(C_d, C_v, A_x, Q, m_0, T_f, n);

figure;
% Plot for Altitude over time on logarithmic scale
subplot(2, 1, 1);
loglog(t, Hs, '--r');
title('Rocket Altitude Over Time (Logarithmic Scale)');
xlabel('Time (s)');
ylabel('Altitude (m)');
grid on;

% Plot for velocity on logarithmic scale
subplot(2, 1, 2);
loglog(t, Vs, '--b');
title('Rocket Velocity Over Time (Logarithmic Scale)');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
grid on;

% Function that performs simulation
function [t, Vs, Hs] = performSim(C_d, C_v, A_x, Q, m_0, T_f, n)

    % Time set-up
    inter = [0, T_f];
    dt = (inter(2) - inter(1)) / n;
    t(1) = inter(1);

    % Initial Conditions
    Vs(1) = 0;
    Hs(1) = 0;

    % Main Euler Method Loop
    for i = 1:n
        t(i + 1) = t(i) + dt;
        [dv, dh, F_t, F_g, F_d, current_m] = rocket_step(t(i), Vs(i), Hs(i), m_0, C_d, C_v, A_x, Q, dt);
       
        Vs(i + 1) = Vs(i) + dv;
        Hs(i + 1) = Hs(i) + dh;
    end
end

function [dv, dh, F_t, F_g, F_d, current_m] = rocket_step(t, v, h, m_0, C_d, C_v, A_x, Q, dt)
    % Universal constants
    G = 6.67430e-11;    % Gravitational Constant
    m_e = 5.972e24;     % Earth's mass (kg)
    R = 6.371e6;        % Earth's radius (m)
    rho_0 = 1.225;      % Air density at sea level (kg/m^3)
    H = 8500;           % Scale height of Earth's atmosphere (m)

    % Mass of Rocket
    current_m = 50 + m_0 * exp(-Q * t);
    dm_dt = -Q * m_0 *exp(-Q*t);

    % Thrust Force
    F_t = -C_v * dm_dt;

    % Gravitational force
    F_g = G * (m_e * current_m) / (R + h)^2;

    % Drag Force
    F_d = 0.5 * rho_0 * exp(-h / H) * C_d * A_x * v^2;

    % Net Force
    net_F = F_t - F_g - F_d;

    % Differential Velocity
    dv = (net_F / current_m) * dt;

    % Differential Height
    dh = v * dt;
end
