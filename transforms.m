% Rocket Launch Initial Conditions in ECI Frame

% Constants
R = 6378137;                  % Earth's equatorial radius in meters
e = 0.08181919;               % Earth's eccentricity
omega_earth = 7.2921159e-5;   % Earth's rotation rate in rad/s

% Input: Launch Site Parameters
latitude = input('Enter the launch site latitude (degrees): ');
longitude = input('Enter the launch site longitude (degrees): ');
altitude = input('Enter the altitude above sea level (meters): ');
v_launch = input('Enter the initial launch velocity [vx, vy, vz] in m/s: ');

% Convert latitude and longitude to radians
lat_rad = deg2rad(latitude);
lon_rad = deg2rad(longitude);

% 1. Calculate ECF Position
R_lat = R / sqrt(1 - e^2 * sin(lat_rad)^2);  % Radius at latitude
x_ECF = (R_lat + altitude) * cos(lat_rad) * cos(lon_rad);
y_ECF = (R_lat + altitude) * cos(lat_rad) * sin(lon_rad);
z_ECF = ((1 - e^2) * R_lat + altitude) * sin(lat_rad);
r_ECF = [x_ECF; y_ECF; z_ECF];

% Display ECF position
disp('Initial Position in ECF coordinates (m):');
disp(r_ECF);

% 2. Time-based Rotation for ECF to ECI Transformation
t = 0;  % Start at t = 0 for initial position
theta_0 = 0;  % Assuming Greenwich Sidereal Time (GST) is 0 at t = 0
theta_t = theta_0 + omega_earth * t;

% Rotation Matrix for ECF to ECI
R_ECF_to_ECI = [cos(theta_t), -sin(theta_t), 0;
                sin(theta_t),  cos(theta_t), 0;
                0,            0,           1];

% 3. Calculate ECI Position
r_ECI = R_ECF_to_ECI * r_ECF;

% Display ECI position
disp('Initial Position in ECI coordinates (m):');
disp(r_ECI);

% 4. Calculate Initial Velocity in ECI due to Earth's Rotation
% The cross product omega_earth x r_ECI gives the velocity component due to Earth's rotation
v_rotation = cross([0; 0; omega_earth], r_ECI);

% Add initial launch velocity in the rocket's FOR (assuming it's in ECF)
v_ECI = v_rotation + R_ECF_to_ECI * v_launch(:);  % Convert launch velocity to ECI frame

% Display ECI velocity
disp('Initial Velocity in ECI coordinates (m/s):');
disp(v_ECI);

% Optional: Display both Position and Velocity in ECI
fprintf('\nInitial Position (ECI):\n');
fprintf('X: %.2f m, Y: %.2f m, Z: %.2f m\n', r_ECI(1), r_ECI(2), r_ECI(3));
fprintf('\nInitial Velocity (ECI):\n');
fprintf('VX: %.2f m/s, VY: %.2f m/s, VZ: %.2f m/s\n', v_ECI(1), v_ECI(2), v_ECI(3));
