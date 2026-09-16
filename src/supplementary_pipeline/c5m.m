%% STRAIN ANALYSIS AND SAFETY VALIDATION MODULE (c5m.m)
% Computes displacement, curvature, stress, and strain from modal data.
% Performs a comprehensive Factor of Safety (FoS) analysis.
%
% This version is robust to mismatches in mode counts between input files.
%
% Dependencies:
%   - modal_response.mat (from c3m.m)
%   - modal_shapes.mat (from c2m.m)
%
% Generates:
%   - strain_analysis_results.mat: A struct with all calculated data.
%   - A series of plots for visualization.
%
% Last Revised: July 17, 2024
clc; clear; close all;

%% 1. Configuration and Constants
fprintf('=== Wing Strain Analysis & Safety Validation ===\n');

% Centralized configuration for analysis
config.material.E = 70e9;              % Young's modulus [Pa] (Aluminum)
config.material.nu = 0.33;             % Poisson's ratio
config.material.sigma_yield = 345e6;   % Yield strength [Pa] (Al 2024-T3)

% **FIX**: Define geometry directly, as it's not in the .mat file
config.geometry.b = 0.15;              % Wing section width [m]
config.geometry.h = 0.02;              % Wing section height [m]

config.analysis.safety_factor_req = 2.5; % Minimum required safety factor

%% 2. Load Input Data from Previous Modules
fprintf('Loading modal data...\n');
if ~exist('modal_response.mat', 'file') || ~exist('modal_shapes.mat', 'file')
    error('Required data files not found. Please run c2m.m and c3m.m first.');
end

try
    % Load modal response data (from c3m.m)
    modal_data = load('modal_response.mat', 'q_sol', 't_sol');
    q_raw = modal_data.q_sol;
    t = modal_data.t_sol;
    
    % Load modal shapes and geometry (from c2m.m)
    shape_data = load('modal_shapes.mat', 'V', 'freqs', 'N', 'L');
    Phi_raw = shape_data.V;         % Mode shapes (N_nodes x n_modes)
    x = linspace(0, shape_data.L, shape_data.N);
    
catch ME
    error('Failed to load or process data: %s', ME.message);
end

%% 3. Data Preprocessing and Synchronization
% **FIX**: Resolve mode count mismatch between input files.
n_modes_in_response = size(q_raw, 2);
n_modes_in_shapes = size(Phi_raw, 2);

% Use the minimum number of modes available from both sources.
n_modes_to_use = min(n_modes_in_response, n_modes_in_shapes);

fprintf('Data loaded. Found %d modes in response data and %d in shape data.\n', ...
        n_modes_in_response, n_modes_in_shapes);
fprintf('--> Synchronizing to use %d modes for analysis.\n', n_modes_to_use);

% Trim data to the synchronized mode count
q = q_raw(:, 1:n_modes_to_use);       % (n_time x n_modes_to_use)
Phi = Phi_raw(:, 1:n_modes_to_use);   % (N_nodes x n_modes_to_use)

%% 4. Calculate Dynamic and Structural Response
% Geometric property for stress calculation
I = (config.geometry.b * config.geometry.h^3) / 12; % Moment of inertia
z_outer = config.geometry.h / 2;               % Max distance from neutral axis

% Reconstruct full displacement field w(x,t)
w = Phi * q'; % Displacement field (N_nodes x n_time)
w(~isfinite(w)) = 0; % Sanity check

% Calculate curvature (d^2w/dx^2) using a 2nd-order finite difference
dx = x(2) - x(1);
curvature = zeros(size(w));
% Central difference for interior points
curvature(2:end-1, :) = (w(3:end, :) - 2*w(2:end-1, :) + w(1:end-2, :)) / dx^2;
% Forward difference for the fixed root (x=0)
curvature(1, :) = (w(3, :) - 2*w(2, :) + w(1, :)) / dx^2;
% Curvature is zero at the free tip (x=L)
curvature(end, :) = 0;

% Calculate Stress and Strain
stress_matrix = -config.material.E * z_outer * curvature;
stress_matrix(~isfinite(stress_matrix)) = 0;
strain_matrix = stress_matrix / config.material.E;

%% 5. Safety and Performance Analysis
safety_factor_stress = config.material.sigma_yield ./ (abs(stress_matrix) + eps);
[max_stress, max_stress_idx] = max(abs(stress_matrix(:)));
[min_sf, min_sf_idx] = min(safety_factor_stress(:));
[max_stress_x_idx, max_stress_t_idx] = ind2sub(size(stress_matrix), max_stress_idx);
[min_sf_x_idx, min_sf_t_idx] = ind2sub(size(safety_factor_stress), min_sf_idx);

%% 6. Visualization
fprintf('Generating analysis plots...\n');
[T_mesh, X_mesh] = meshgrid(t, x);

% Create a new figure window for all subplots
figure('Name', 'Comprehensive Strain & Safety Analysis', 'NumberTitle', 'off', 'WindowState', 'maximized');

% Plot 1: Spanwise Stress Distribution at Peak Time
subplot(2,2,1);
[~, peak_time_idx] = max(max(abs(stress_matrix), [], 1));
peak_time = t(peak_time_idx);
plot(x, stress_matrix(:, peak_time_idx)/1e6, 'b', 'LineWidth', 2);
hold on;
yline(config.material.sigma_yield/1e6, 'k--', 'Yield Strength');
yline(-config.material.sigma_yield/1e6, 'k--');
xlabel('Spanwise Location [m]'); ylabel('Stress [MPa]');
title(sprintf('Stress Distribution at Peak Response (t=%.2f s)', peak_time));
legend('Actual Stress', 'Yield Strength', 'Location', 'southeast');
grid on; axis tight;

% Plot 2: 3D Stress Surface
subplot(2,2,2);
surf(T_mesh, X_mesh, stress_matrix/1e6, 'EdgeColor', 'none');
xlabel('Time [s]'); ylabel('Span Position [m]'); zlabel('Stress [MPa]');
title('3D Stress Distribution Over Time');
colorbar; colormap('jet'); view(45, 30); shading interp;
ylabel(colorbar, 'Stress [MPa]');

% Plot 3: Factor of Safety Heatmap
subplot(2,2,3);
imagesc(t, x, safety_factor_stress);
colorbar;
title('Factor of Safety (FoS) Distribution');
xlabel('Time [s]'); ylabel('Span Position [m]');
colormap('hot');
caxis([0, min(15, max(safety_factor_stress(:)))]);
ylabel(colorbar, 'Factor of Safety');
hold on;
plot(t(min_sf_t_idx), x(min_sf_x_idx), 'g+', 'MarkerSize', 15, 'LineWidth', 3);
legend('Min FoS Location', 'Location', 'northwest');

% Plot 4: Minimum Safety Factor Along Span
subplot(2,2,4);
plot(x, min(safety_factor_stress, [], 2), 'r-', 'LineWidth', 2);
hold on;
yline(config.analysis.safety_factor_req, 'k--', 'Required FoS');
xlabel('Span Position [m]'); ylabel('Minimum Safety Factor');
title('Minimum FoS Along Span');
legend('Actual Min FoS', 'Required FoS', 'Location', 'best');
grid on; axis tight;

%% 7. Summary Report
fprintf('\n--- ANALYSIS SUMMARY ---\n');
fprintf('  Max stress      : %.1f MPa at x=%.2f m, t=%.2f s\n', ...
    max_stress/1e6, x(max_stress_x_idx), t(max_stress_t_idx));
fprintf('  Min safety factor : %.2f at x=%.2f m, t=%.2f s\n', ...
    min_sf, x(min_sf_x_idx), t(min_sf_t_idx));

if min_sf >= config.analysis.safety_factor_req
    fprintf('  STATUS: ✓ DESIGN PASSES safety requirement (FoS >= %.1f)\n', config.analysis.safety_factor_req);
else
    fprintf('  STATUS: ✗ DESIGN FAILS safety requirement (FoS < %.1f)\n', config.analysis.safety_factor_req);
end
fprintf('-------------------------\n');