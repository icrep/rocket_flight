% Main script
clear; clc;

T_f = input('Enter total time of simulation (s): \n');    % Total duration of the simulation.
n = input('Enter number of time steps: \n');  % Number of time steps

% Time Interval Set-Up
dt = T_f/n;
t(1) = 0;

% Initialize Position, Velocity, and Acceleration vectors (3D)
r(:,1) = [0; 0; 6.371e6];  % Initial position (at Earth's surface along z-axis)
v(:,1) = [100; 100; 500];  % Initial velocity (3D, with some components in x, y, z directions)
a(:,1) = [0; 0; 0];        % Initial acceleration

for i = 1:n
    t(i+1) = t(i) + dt;

    % Use Runge-Kutta to calculate the next step
    [v(:,i+1), r(:,i+1), a(:,i)] = RK4tStep(t(i), v(:,i), r(:,i), dt);

    % Stop condition if altitude returns to or goes below 0
    if r(3,i+1) <= 6.371e6
        break;
    end
end

% Trim arrays to the actual length of the simulation
t = t(1:i+1);
r = r(:,1:i+1);
v = v(:,1:i+1);
a = a(:,1:i);

% Plotting the 3D trajectory
figure
plot3(r(1,:), r(2,:), r(3,:) - 6.371e6, '-k');  % Plot altitude above Earth's surface
grid on;
title('Rocket 3D Trajectory');
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Altitude (m)');
legend({'Trajectory'}, 'Location', 'best');

% Plotting the altitude (z-component of position vector)
figure
subplot(3,1,1);
plot(t, r(3,:) - 6.371e6, '-b');  % Subtract Earth's radius to plot altitude
title('Rocket Altitude');
xlabel('Time (s)');
ylabel('Altitude (m)');
legend({'Altitude'}, 'Location', 'Southeast');

% Plotting the velocity (magnitude)
subplot(3,1,2);
plot(t, vecnorm(v, 2, 1), '-r');
title('Rocket Velocity');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend({'Velocity Magnitude'}, 'Location', 'Southeast');

% Plotting the acceleration (magnitude)
subplot(3,1,3);
plot(t(1:end-1), vecnorm(a, 2, 1), '-g');  % Magnitude of acceleration
title('Rocket Acceleration');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
legend({'Acceleration Magnitude'}, 'Location', 'Southeast');

%% Function Definitions

function [dv, dr, acc] = vdot(t, v, r)
    % Universal constants
    G = 6.67430e-11;    % Gravitational Constant
    m_e = 5.972e24;     % Earth's mass (kg)
    R = 6.371e6;        % Earth's radius (m)
    rho_0 = 1.225;      % Air density at sea level (kg/m^3)
    H = 8500;           % Scale height of Earth's atmosphere (m)
    
    % Rocket Constants
    C_v = -4000;    % Rocket's exhaust velocity (m/s)
    C_d = 0.5;      % Aerodynamic drag coefficient
    A_x = 1;        % Cross-sectional area of the rocket (m^2)
    Q = 0.2;        % Rate of fuel consumption (s^-1)
    m_p = 50;       % Mass of fuel (kg)
    m_d = 50;       % Mass of rocket (kg)
    m_0 = 10 * m_d;

    % Current mass calculation
    current_m = m_d + (m_0 - m_d)*exp(-Q * t);
    dm = -Q * (m_0 - m_d) * exp(-Q * t);

    % Thrust force (aligned with the velocity initially)
    if norm(v) > 0
        F_t = C_v * dm * (v / norm(v));  % Thrust direction matches velocity direction
    else
        F_t = [0; 0; abs(C_v * dm)];  % Initial thrust directed upward if velocity is zero
    end

    % Gravity force (points toward Earth's center)
    r_mag = norm(r);  % Distance from Earth's center
    F_g = -((G * m_e * current_m) / (r_mag^2)) * (r / r_mag);  % Gravity vector

    % Drag force (opposes velocity direction)
    if norm(v) > 0
        F_d = -0.5 * rho_0 * exp(-(norm(r) - R)/H) * C_d * A_x * norm(v)^2 * (v / norm(v));  % Drag vector
    else
        F_d = [0; 0; 0];  % No drag if velocity is zero
    end

    % Net force calculation
    netForce = F_t + F_g + F_d;
    dv = netForce / current_m;  % Acceleration (dv/dt)
    dr = v;                      % Change in position (dr/dt)
    acc = dv;                    % Return acceleration as dv for storing and plotting
end

function [v_next, r_next, acc] = RK4tStep(t, v, r, dt)
    % Runge-Kutta method for 3D vector-based simulation
    [k1_v, k1_r, a1] = vdot(t, v, r);
    [k2_v, k2_r, a2] = vdot(t + dt/2, v + (dt*k1_v)/2, r + (dt*k1_r)/2);
    [k3_v, k3_r, a3] = vdot(t + dt/2, v + (dt*k2_v)/2, r + (dt*k2_r)/2);
    [k4_v, k4_r, a4] = vdot(t + dt, v + dt*k3_v, r + dt*k3_r);

    dv = (dt/6) * (k1_v + 2*k2_v + 2*k3_v + k4_v);
    dr = (dt/6) * (k1_r + 2*k2_r + 2*k3_r + k4_r);

    v_next = v + dv;
    r_next = r + dr;
    acc = (a1 + 2*a2 + 2*a3 + a4) / 6;  % Average acceleration during the step
end
