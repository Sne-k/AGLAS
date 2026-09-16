% =========================================================================
% SCRIPT 3: SIMULATE GUST RESPONSE
% =========================================================================
% Description:
% Simulates the dynamic response of the wing to a "1-cos" gust load using
% modal superposition. The equations of motion are solved in modal
% coordinates using ODE45.
% =========================================================================

% --- Load Model
load('structural_model.mat', 'M', 'K', 'V_eigenvectors', 'V_deflection', 'freqs_hz', 'x_span');

% --- Setup Modal System
n_modes = 3; % Number of modes to include in the simulation

% Use the FULL eigenvectors for modal transformation
V_modes_full = V_eigenvectors(:, 1:n_modes);

% Calculate generalized modal properties
M_modal = diag(diag(V_modes_full' * M * V_modes_full)); % Should be ~Identity
K_modal = diag(diag(V_modes_full' * K * V_modes_full)); % Should be omega^2

% Define modal damping (2% of critical damping for each mode)
zeta = 0.02;
omega = 2 * pi * freqs_hz(1:n_modes);
C_modal = diag(2 * zeta * omega);

% --- Define Gust Load
V_gust = 10;  % Peak gust velocity (m/s)
t_gust = 2;   % Gust duration (s)

% Create a spatial force distribution (uniform pressure)
% Force is applied only to translational DOFs, not rotational
force_distribution = zeros(size(K, 1), 1);
force_distribution(1:2:end) = 1; % Apply unit force to deflection DOFs

% Transform spatial force into generalized modal forces
F_generalized = V_modes_full' * force_distribution;
gust_force_function = @(t) (t <= t_gust) * (0.5 * V_gust * (1 - cos(pi * t / t_gust))) * F_generalized;

% --- Solve Equations of Motion in Modal Coordinates
odefun = @(t, x) [x(n_modes+1:end); M_modal \ (gust_force_function(t) - C_modal*x(n_modes+1:end) - K_modal*x(1:n_modes))];
x0 = zeros(2 * n_modes, 1); % Initial conditions: zero displacement and velocity
tspan = [0 8];              % Simulation time
[t_sol, x_sol] = ode45(odefun, tspan, x0);

% --- Reconstruct Physical Displacement
q_sol = x_sol(:, 1:n_modes); % Extract modal displacements (generalized coordinates)
V_modes_deflection = V_deflection(:, 1:n_modes); % Get deflection shapes for the modes used
wing_tip_disp = q_sol * V_modes_deflection(end, :)'; % Reconstruct tip displacement

% --- Plot Wing Tip Response
figure('Name', 'Wing Tip Gust Response', 'NumberTitle', 'off');
plot(t_sol, wing_tip_disp, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Wing Tip Displacement (m)');
title('Wing Tip Dynamic Response to Gust Load');
grid on;

% --- Save Response Data
save('modal_response.mat', 'q_sol', 't_sol', 'n_modes');