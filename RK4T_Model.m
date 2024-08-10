 % Main script
clear; clc;


T_f = input('Enter total time of simulation (s): \n');    % Total duration of the simulation.
n = input('Enter number of time steps: \n');  % Time step size for the simulation.

%Time Interval Set-Up
dt = T_f/n;
t(1) = 0;

% Initialize Velocity and Height vectors
Hs(1) = 0;
Vs(1) = 0;

Hse(1) = 0;
Vse(1) = 0;

for i = 1:n
    t(i+1) = t(i) + dt;

    % Call the rkstep function to call the vdot function
    [Vs(i+1), Hs(i+1)] = RK4tStep(t(i), Vs(i), Hs(i), dt);

    % Compare to Euler's Method
    [Vse(i+1), Hse(i+1)] = eulerRocket(t(i), Vse(i), Hse(i), dt);
end


figure
subplot(2, 1, 1)
plot(t, Vs, '-r', t, Vse, '--r');
title('Rocket Velocity');
xlabel('Time(s)');
ylabel('Velocity (m/s)');
legend({'RK4 Method', 'Euler Method'}, 'Location', 'Southeast');

subplot(2,1,2)
plot(t, Hs, '-b', t, Hse, '--b');
title('Rocket Altitutude');
xlabel('Time(s)');
ylabel('Altitude (m)');
legend({'RK4 Method', 'Euler Method'}, 'Location', 'Southeast');
 
function [dv, dh] = vdot(t, v, h)
    % Universal constants
    G = 6.67430e-11;    % Gravitational Constant
    m_e = 5.972e24;     % Earth's mass (kg)
    R = 6.371e6;        % Earth's radius (m)
    rho_0 = 1.225;      % Air density at sea level (kg/m^3)
    H = 8500;           % Scale height of Earth's atmosphere (m)
    
    % Rocket Constants
    C_v = -4000;    % Rocket's exhaust velocity (m/s)
    C_d = 0.5;     % Aerodynamic drag coefficient
    A_x = 1;       % Cross-sectional area of the rocket (m^2)
    Q = 0.2;      % Rate of fuel consumption (s^-1)
    m_p = 50;     % Mass of fuel (kg)
    m_d = 50;  % Mass of rocket (kg) (to ensure proper ratio)
    m_0 = 10*m_d;

    current_m = m_d + (m_0 - m_d)*exp(-Q * t);
    dm = -Q*(m_0 - m_d )*exp(-Q*t);

    % Thrust, only positive force
    F_t = C_v*dm;
    
    % Gravity
    F_g = (G*m_e*current_m) / (R + h)^2 ;

    % Air Resistance, dependant on altitude and speed
    F_d = (.5*rho_0*exp(-h/H)*C_d*A_x*(v^2));

    % Sum of all forces...etc
    netForce = (F_t - F_g - F_d);
    dv = (netForce/current_m);

    % derivative of height is velocity
    dh = v;
end

function [v_next, h_next] = RK4tStep(t, v, h, dt)
    k1_v = vdot(t, v, h);
    k2_v = vdot(t + dt/2, v + (dt*k1_v)/2, h);
    k3_v = vdot(t + dt/2, v + (dt*k2_v)/2, h);
    k4_v = vdot(t + dt, v + dt*k3_v, h);

    dv = (dt/6) * (k1_v + 2*k2_v + 2*k3_v + k4_v);

    % calculates height with current velocity
    % readjust for efficiency
    dh = v * dt;

    v_next = v + dv;
    h_next = h + dh; 
end 

function [v_next, h_next] = eulerRocket(t, v, h, dt)
    [dv, dh] = vdot(t, v, h);
    v_next = v + dt*dv;
    h_next = h + dt*dh;
end 