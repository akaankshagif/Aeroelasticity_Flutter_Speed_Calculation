
clear
clc

% Define symbolic variables
syms i L y C C1;

% Boundry condition  (y = 0) to resolve C1
phi_i_PP = (1/L^3) * (L - y) * (y/L)^(i-1); %defining the equation for phi i double prime
phi_i_P = int(phi_i_PP, y) + C; 
findC = solve(subs(phi_i_P, y, 0) == 0, C);
phi_i_P = subs(phi_i_P, C, findC);

% Applying boundary condition at y = 0 to find C2 from phi_i equation
phi_i = int(phi_i_P, y) + C1;
findC1 = solve(subs(phi_i, y, 0) == 0, C1);
phi_i = subs(phi_i, C1, findC1);

% defining the constants
L = 30; % Length of the wing in meters (given)
E = 2*10^11; % Young's modulus in N/m^2(given)
rho = 1600; % Density in kg/m^3        (given)
m = 6000; % approx mass of the trent1000 engine in kg (assumed)

% Create syms var for integration called y 
syms y;
Y = linspace(0, L, 100); % so as to divide the wing into seperate points for plotting
% Define the cross-sectional area (A) and moment of inertia (I)
A = -60.78 * (y/L)^3 + 121.16 * (y/L)^2 - 79.624 * (y/L) + 19.344; 
I = 2.7462 * exp(-7.958 * (y/L)); 

% INPUT THE NUMBER OF MODES HERE (Task 2)
n_modes = 8; % defining the number of modes for calculation

% Allocate M and K as matrices
M = zeros(n_modes, n_modes);
K = zeros(n_modes, n_modes);
m_noengine= zeros(n_modes,n_modes);
% Loop through and calculate mass and stiffness matrices
for i = 1:n_modes
    for j = 1:n_modes
        % Define shape functions
        phi_i_PP = (1/L^3) * (L - y) * (y/L)^(i-1);
        phi_j_PP = (1/L^3) * (L - y) * (y/L)^(j-1);

        phi_i = (y/L)^(i+1) * (2 + i - i * (y/L)) / (i * (i + 1) * (i + 2));
        phi_j = (y/L)^(j+1) * (2 + j - j * (y/L)) / (j * (j + 1) * (j + 2));

        % Mass matrix for the wing (using symbolic integration)
        m_wing = double(int(rho * A * phi_i * phi_j, y, 0, L));

        % Evaluate shape functions at the engine location (L/3)
        phi_i_engine = double(subs(phi_i, y, L/3));
        phi_j_engine = double(subs(phi_j, y, L/3));
        m_engine =  m * phi_i_engine * phi_j_engine;
        % Total mass matrix including engine
        M(i,j) = m_wing + m_engine;
        
        m_noengine(i,j)= m_wing;
        % Stiffness matrix (using symbolic integration)
        K(i,j) = double(int(E * I * phi_i_PP * phi_j_PP, y, 0, L));
    end
end

% Solve the eigenvalue problem for natural frequencies
[V, D] = eig(K, M);

% Sort eigenvalues and eigenvectors
[D, index] = sort(diag(D));
V_sor = V(:, index);

% Natural frequencies (angular)
W = sqrt(D);

% Convert to frequency in Hz
f = W / (2 * pi);
f4 = f(4); % Base case: 4th natural frequency

%% without engine
[Vhash, Dhash] = eig(K, m_noengine);
% Sort eigenvalues and eigenvectors for no engine scenario
[Dhash, index] = sort(diag(Dhash));
Vh_sorted = Vhash(:, index);

% Natural frequencies with no engine
Wknot = sqrt(Dhash);

% Convert to frequency in Hz
fknot = Wknot / (2 * pi);
fknot4 = fknot(4); % Base case: 4th natural frequency


%% Calculation loop to ensure no error (Task 2)
error = 1;
err_no_engine =1;

while abs(error) > 0.01 * f4
while abs(err_no_engine)>0.01 * fknot4
    % Re-initialize matrices
    M = zeros(n_modes, n_modes);
    K = zeros(n_modes, n_modes);
    m_noengine = zeros(n_modes,n_modes);
    for i = 1:n_modes
        for j = 1:n_modes
            % Recalculate shape functions for each iteration
            phi_i_PP = (1/L^3) * (L - y) * (y/L)^(i-1);
            phi_j_PP = (1/L^3) * (L - y) * (y/L)^(j-1);

            phi_i = (y/L)^(i+1) * (2 + i - i * (y/L)) / (i * (i + 1) * (i + 2));
            phi_j = (y/L)^(j+1) * (2 + j - j * (y/L)) / (j * (j + 1) * (j + 2));

            % Mass matrix integration
            m_wing = double(int(rho * A * phi_i * phi_j, y, 0, L));

            % Evaluate shape functions at engine location (L/3)
            phi_i_engine = double(subs(phi_i, y, L/3));
            phi_j_engine = double(subs(phi_j, y, L/3));
            m_eng = m * phi_i_engine * phi_j_engine;
            % Total mass matrix including engine
            M(i, j) = m_eng + m_wing ;
            m_noengine(i,j)=m_wing;
            % Stiffness matrix integration
            K(i, j) = double(int(E * I * phi_i_PP * phi_j_PP, y, 0, L));
        end
    end

    % Solve the eigenvalue problem for with engine
    [V, D] = eig(K, M);
    [D, index] = sort(diag(D));
    V_sorted = V(:, index);
     

    % Calculate natural frequencies
    W = sqrt(D);
    f = W / (2 * pi);

    % Update error and fourth natural frequency
    newf4 = f(4);
    error = newf4 - f4;
    f4 = newf4;

    %% without engine
    [Vhash, Dhash] = eig(K, m_noengine);

    % Sort eigenvalues and eigenvectors
    [Dhash, index] = sort(diag(Dhash));
    Vh_sorted = Vhash(:, index);

    % Natural frequencies (angular)
    Wknot = sqrt(Dhash);

    % Convert to frequency in Hz
    fknot = Wknot / (2 * pi);
    Fknew4 = fknot(4);
    err_no_engine = Fknew4 - fknot4;
    fknot4 = Fknew4;
end
end

% Display the first four natural frequencies
for mode_num = 1:4
    fprintf('Natural frequency for mode %d is %.3f Hz\n', mode_num, f(mode_num));
end

fprintf('Converged fourth natural frequency is %.3f Hz\n', f4);

%display the first four natural frq without engine
for mode_num = 1:4
    fprintf('Natural frequency WITHOUT ENGINE for mode %d is %.3f Hz\n', mode_num, fknot(mode_num));
end

fprintf('Converged fourth natural frequency without engine is %.3f Hz\n', fknot4);

%% Plotting mode shapes with engine
figure(1);
hold on;

for j = 1:4
    mode1 = zeros(1, 100); % Preallocate for mode shape calculation
    for i = 1:n_modes
        % Define the shape function symbolically
        phi_i_plot = (y/L)^(i+1) * (2 + i - i * (y/L)) / (i * (i + 1) * (i + 2));
        phi_i_eval = double(subs(phi_i_plot, y, Y)); % Evaluate for numeric points

        % Add contribution to mode shape
        mode1 = mode1 + V_sorted(i, j) * phi_i_eval;
    end
    % Normalize and plot the mode shape
    plot(Y, mode1 / mode1(end), 'DisplayName', ['Mode ', num2str(j)]);
end

% Label the plot
title('Mode Shapes with Engine at L/3');
xlabel('Span (m)');
ylabel('Normalized Displacement');
yline(0, '--');
xlim([0, L]);
legend('Location', 'northwest');

hold off;

%% Plotting Mode Shapes (Without Engine)
figure(2); % Ensure this figure is created
hold on;

for j = 1:4
    mode2 = zeros(1, 100); % Preallocate for mode shape calculation
    for i = 1:n_modes
        % Define the shape function symbolically
        phi_i_plot = (y/L)^(i+1) * (2 + i - i * (y/L)) / (i * (i + 1) * (i + 2));
        phi_i_eval = double(subs(phi_i_plot, y, Y)); % Evaluate for numeric points

        % Add contribution to mode shape
        mode2 = mode2 + Vh_sorted(i, j) * phi_i_eval;
    end
    % Normalize and plot the mode shape
    plot(Y, mode2 / mode2(end), 'DisplayName', ['Mode ', num2str(j)]);
end

% Label the plot
title('Mode Shapes WITHOUT Engine');
xlabel('Span (m)');
ylabel('Normalized Displacement');
yline(0, '--');
xlim([0, L]);
legend('Location', 'northwest');

hold off;

%% Plotting Overlapping Mode Shapes
figure(3); % Ensure this figure is created
hold on;

for j = 1:4
    mode1 = zeros(1, 100); % Preallocate for mode shape calculation with engine
    mode2 = zeros(1, 100); % Preallocate for mode shape calculation without engine
    for i = 1:n_modes
        % Define the shape function symbolically
        phi_i_plot = (y/L)^(i+1) * (2 + i - i * (y/L)) / (i * (i + 1) * (i + 2));
        phi_i_eval = double(subs(phi_i_plot, y, Y)); % Evaluate for numeric points

        % Add contribution to mode shape with engine
        mode1 = mode1 + V_sorted(i, j) * phi_i_eval;

        % Add contribution to mode shape without engine
        mode2 = mode2 + Vh_sorted(i, j) * phi_i_eval;
    end
    % Normalize and plot both mode shapes
    plot(Y, mode1 / mode1(end), 'DisplayName', ['Mode ', num2str(j), ' With Engine']);
    plot(Y, mode2 / mode2(end), 'DisplayName', ['Mode ', num2str(j), ' Without Engine'], 'LineStyle', '--');
end

% Label the plot
title('Overlapping Mode Shapes');
xlabel('Span (m)');
ylabel('Normalized Displacement');
yline(0, '--');
xlim([0, L]);
legend('Location', 'northwest');

hold off;
