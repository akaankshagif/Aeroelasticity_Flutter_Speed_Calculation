% Clear workspace
clc;
clear;

% Constants and initialization
% %notes: Define material properties and geometric parameters for the aeroelastic system
material_density = 800; % Density of the material (kg/m^3)
chord_length = 2; % Chord length of the aerofoil (m)
x_f = 0.45 * chord_length; % Position of aerodynamic center related to chord length
a = -0.05 * chord_length; % Offset from elastic axis
thickness = 0.01; % Thickness of the aerofoil (m)
length_fuselage = 3; % Length of the fuselage (m)
omega_h = 2 * 2 * pi; % Plunge frequency in rad/s
omega_a = 6 * 2 * pi; % Pitch frequency in rad/s
mass = material_density * chord_length * thickness; % Mass of the aerofoil per unit span (kg/m)

% Moment of inertia and stiffness values
% %notes: Calculate critical structural parameters for dynamic stability analysis
inertia_a = (mass * chord_length^2 / 12) + mass * a^2; % Moment of inertia about elastic axis
s = -mass * a; % Coupling term for plunge-pitch dynamics
K_h = mass * omega_h^2; % Plunge stiffness
K_a = inertia_a * omega_a^2; % Pitch stiffness

% Aerodynamic parameters
% %notes: Include air properties and aerodynamic span parameters
air_density = 1.225; % Density of air (kg/m^3)
b = chord_length / 2; % Half chord length (m)
identity_matrix = eye(2); % Identity matrix for state-space construction
zero_matrix = zeros(2, 2); % Zero matrix for initialization
air_speed = 1:0.1:100; % Range of airspeeds (m/s) to evaluate stability
initial_guess = [omega_h, omega_a]; % Initial frequency guesses for iterative solution
frequency_results = []; % Array to store natural frequencies
damping_ratios = []; % Array to store damping ratios

% Main loop for calculating frequencies and damping ratios
% %notes: Loop through airspeeds and compute dynamic parameters for each speed
for speed_index = 1:length(air_speed)
    speed_value = air_speed(speed_index); % Current airspeed
    for dof = 1:2 % Loop through plunge and pitch degrees of freedom
        omega = initial_guess(dof); % Initial frequency guess for current DOF
        error_val = 100; % Initialize error to start iterative solver
        
        % %notes: Refine frequency using iterative solver
        while abs(error_val) > 0.001 % Convergence criterion
            k = omega * b / speed_value; % Reduced frequency
            Ck = 1 - (0.165 / (1 - (0.0455 / k) * 1i)) - (0.335 / (1 - (0.30 / k) * 1i)); % Unsteady aerodynamic function
            
            % Construct system matrices
            stiffness_matrix = [K_h, 0; 0, K_a]; % Stiffness matrix
            mass_matrix = [mass, s; s, inertia_a]; % Mass matrix
            D_matrix = (air_density * pi * b^2) * [-1, a; a, -(a^2 + (b^2) / 8)]; % Damping matrix
            E_matrix = (air_density * pi * speed_value * b) * [-2 * Ck, -b + (2 * Ck * (a - b / 2)); 
                b - (b - (2 * a + b) * Ck), -b * (chord_length / 4) + (b - (2 * a + b) * Ck) * (a - (b / 2))];
            F_matrix = (-air_density * pi * speed_value^2 * b) * [0, 2 * Ck; 0, -((2 * a + b) * Ck)];
            A_matrix = [zero_matrix, identity_matrix; (mass_matrix - D_matrix) \ (F_matrix - stiffness_matrix), (mass_matrix - D_matrix) \ E_matrix];

            % Eigenvalue calculation and sorting
            eigenvalues = eig(A_matrix); % Calculate eigenvalues
            [~, idx] = sort(imag(eigenvalues), 'ascend'); % Sort eigenvalues by imaginary part
            eigenvalues = eigenvalues(idx); % Reorder eigenvalues
            new_omega = abs(eigenvalues(dof + 2)); % Extract updated frequency
            error_val = abs(new_omega - omega); % Compute error for convergence
            omega = new_omega; % Update frequency for next iteration
        end
        
        % %notes: Store results for natural frequencies and damping ratios
        frequency_results(dof, speed_index) = abs(eigenvalues(dof + 2)); % Natural frequency
        damping_ratios(dof, speed_index) = -real(eigenvalues(dof + 2)) / abs(eigenvalues(dof + 2)); % Damping ratio
    end
end

% Detect flutter
% %notes: Identify the airspeed where flutter occurs by observing sign changes in damping ratio
flutter_speed = NaN; % Initialize flutter speed
flutter_frequency = NaN; % Initialize flutter frequency
for i = 1:length(damping_ratios)-1
    if (sign(damping_ratios(2, i)) * sign(damping_ratios(2, i+1))) == -1 % Check for sign change
        flutter_speed = air_speed(i);
        flutter_frequency = frequency_results(2, i);
        break; % Exit loop after detecting flutter
    end
end

% Plotting frequency results
% %notes: Visualize natural frequencies as a function of airspeed
figure;
hold on;
plot(air_speed, frequency_results(1, :), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Plunge Frequency');
plot(air_speed, frequency_results(2, :), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Pitch Frequency');
plot(flutter_speed, flutter_frequency, 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Flutter Frequency');
plot([flutter_speed, flutter_speed], [-5, 30], 'k--', 'LineWidth', 1);
title('Natural Frequency vs. Airspeed');
xlabel('Airspeed (m/s)');
ylabel('Natural Frequency (rad/s)');
legend('Location', 'northeast');
grid on;

% Plotting damping ratios
% %notes: Visualize stability trends by plotting damping ratios
figure;
hold on;
plot(air_speed, damping_ratios(1, :), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Plunge Damping Ratio');
plot(air_speed, damping_ratios(2, :), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Pitch Damping Ratio');
plot(flutter_speed, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Flutter Speed');
plot([0, max(air_speed)], [0, 0], 'k--', 'LineWidth', 1);
ylim([-0.2, 1]);
title('Damping Ratio vs. Airspeed');
xlabel('Airspeed (m/s)');
ylabel('Damping Ratio, \xi');
legend('Location', 'northeast');
grid on;
