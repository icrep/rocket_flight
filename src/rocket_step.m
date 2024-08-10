function [dv, dh] = rocket_step(t, v, h, m_0, C_d, C_v, A_x, Q, dt) 
    % Universal constants
    G = 6.67430e-11;    % Gravitational Constant
    m_e = 5.972e24;     % Earth's mass (kg)
    R = 6.371e6;        % Earth's radius (m)
    rho_0 = 1.225;      % air density at sea level ()
    H_0 = 8500;         % scale height of Earth's atmosphere (m)

    % Mass of Rocket
    current_m = m_0*exp(-Q*t);

    % Thrust Force
    F_t = (-C_v*-Q)*current_m;

    % Gravitational force
    F_g = G*((m_e*current_m))/((R+h)^2);

    % Drag Force
    F_d = 0.5*(rho_0)*exp(-(h)/H_0)*C_d*A_x*(v^2);

    % Differential Velocity
    dv = dt*((F_t - F_g - F_d)/(current_m));

    % Differential Height
    dh = v * dt;


end