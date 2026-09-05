clc;
clear;
close all;

%% 1. LOAD WING STRUCTURAL MODEL
fprintf('Loading structural model from c2m.m...\n');
if ~exist('structural_model.mat', 'file')
    error('structural_model.mat not found. Please run the corrected c2m.m script first.');
end
% **MODIFIED**: Loads the new V_eigenvectors and V_deflection variables
load('structural_model.mat', 'M', 'K', 'V_eigenvectors', 'V_deflection', 'freqs', 'L', 'N');

%% 2. SETUP MODAL PARAMETERS
fprintf('Setting up modal system for dynamic analysis...\n');
n_modes = 3; % Number of modes to include in the simulation

% **MODIFIED**: Use the FULL eigenvectors for modal transformation
V_modes_full = V_eigenvectors(:, 1:n_modes);

% Calculate generalized modal properties using the full eigenvectors
M_modal = diag(diag(V_modes_full' * M * V_modes_full));
K_modal = diag(diag(V_modes_full' * K * V_modes_full));

% Define modal damping
zeta = 0.02; % 2% damping ratio
omega = 2 * pi * freqs(1:n_modes);
C_modal = diag(2 * zeta * omega);

%% 3. DEFINE GUST LOAD
V_gust = 10; t_gust = 2;

% **MODIFIED**: The force distribution must match the size of the full system (M and K)
% We create a force vector that applies force only to the deflection degrees of freedom
force_distribution = zeros(size(K, 1), 1);
force_distribution(1:2:end) = 1; % Apply force to deflection DOFs, not rotation

% Transform the spatial force into a generalized modal force using full eigenvectors
F_generalized = V_modes_full' * force_distribution;

% Gust load function
gust_force_function = @(t) (t <= t_gust) * (0.5 * V_gust * (1 - cos(pi * t / t_gust))) * F_generalized;

%% 4. SOLVE THE EQUATIONS OF MOTION
% (No changes in this section)
odefun = @(t, x) [x(n_modes+1:end); M_modal \ (gust_force_function(t) - C_modal*x(n_modes+1:end) - K_modal*x(1:n_modes))];
x0 = zeros(2*n_modes, 1);
tspan = [0 5];
[t_sol, x_sol] = ode45(odefun, tspan, x0);

%% 5. RECONSTRUCT PHYSICAL DISPLACEMENT
fprintf('Reconstructing physical displacement...\n');
q_sol = x_sol(:, 1:n_modes);

% **MODIFIED**: Use the deflection-only shapes to reconstruct the final displacement
V_modes_deflection = V_deflection(:, 1:n_modes);
V_tip = V_modes_deflection(end, :); % Deflection shape at the wing tip
wing_tip_disp = q_sol * V_tip';

%% 6. PLOT AND SAVE RESULTS
% (No changes in this section)
fprintf('Plotting results and saving data...\n');
figure('Name', 'Wing Tip Gust Response', 'NumberTitle', 'off');
plot(t_sol, wing_tip_disp, 'r', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Wing Tip Displacement (m)');
title('Realistic Wing Tip Response to Gust Load'); grid on;

% Save results for c5m.m
save('modal_response.mat', 'q_sol', 't_sol', 'n_modes');
% Use the existing 'modal_shapes.mat' name but with the corrected variables
V = V_deflection;
save('modal_shapes.mat', 'V', 'freqs', 'N', 'L');

fprintf('Analysis complete. Response data saved.\n');