clc; clear; close all;

%% 1. WING STRUCTURAL MODEL - FINITE ELEMENT METHOD (FEM)
fprintf('Building structural model...\n');

% Parameters
L = 10;      % Wing span (m)
E = 70e9;    % Young's modulus (Pa)
I = 0.02;    % Moment of inertia (m^4)
rho = 1600;  % Material density (kg/m^3)
A = 0.05;    % Cross-sectional area (m^2)
N = 10;      % Number of finite elements

dx = L / N;  % Length of each element

% Define local element stiffness and mass matrices for a beam element
% (4x4 matrix for 2 nodes, each with 2 degrees of freedom: deflection & rotation)
k_local = (E * I / dx^3) * [12,   6*dx,    -12,   6*dx;
                            6*dx, 4*dx^2,  -6*dx, 2*dx^2;
                           -12,  -6*dx,   12,    -6*dx;
                            6*dx, 2*dx^2,  -6*dx, 4*dx^2];

m_local = (rho * A * dx / 420) * [156,   22*dx,   54,    -13*dx;
                                  22*dx, 4*dx^2,  13*dx,  -3*dx^2;
                                  54,    13*dx,   156,    -22*dx;
                                 -13*dx, -3*dx^2, -22*dx,  4*dx^2];

% Initialize global matrices (2 degrees of freedom per node)
n_nodes = N + 1;
K_global = zeros(2 * n_nodes);
M_global = zeros(2 * n_nodes);

% Assemble global matrices from local element matrices
for i = 1:N
    dof_indices = (2*i-1):(2*i+2); % Get the 4 relevant DOFs for the element
    K_global(dof_indices, dof_indices) = K_global(dof_indices, dof_indices) + k_local;
    M_global(dof_indices, dof_indices) = M_global(dof_indices, dof_indices) + m_local;
end

%% 2. APPLY BOUNDARY CONDITIONS
% For a cantilever beam, the wing root is fixed.
% This means deflection and rotation at the first node are zero.
% We remove the first two rows and columns, corresponding to these DOFs.
dofs_to_keep = 3:(2 * n_nodes);
K = K_global(dofs_to_keep, dofs_to_keep);
M = M_global(dofs_to_keep, dofs_to_keep);

%% 3. SOLVE EIGENVALUE PROBLEM
% Solve (K - omega^2 * M)v = 0 for natural frequencies and mode shapes
[V_raw, D] = eig(K, M);
omegas = sqrt(diag(D)); % Angular frequencies (rad/s)
freqs = omegas / (2 * pi); % Convert to Hz

% Normalize mode shapes for consistency
for i = 1:size(V_raw, 2)
    V_raw(:,i) = V_raw(:,i) / max(abs(V_raw(:,i)));
end

% The eigenvector `V_raw` contains both deflection and rotation.
% Extract only the deflection mode shapes for analysis.
V = V_raw(1:2:end, :);

%% 4. DISPLAY AND PLOT RESULTS
fprintf('Natural Frequencies of Wing (Hz):\n');
disp(freqs(1:5)); % Display first 5 frequencies

% Plot the first 3 mode shapes
figure('Name', 'Wing Mode Shapes', 'NumberTitle', 'off');
x_nodes = linspace(0, L, n_nodes);
for i = 1:3
    subplot(3, 1, i);
    plot(x_nodes, [0; V(:,i)], 'o-', 'LineWidth', 1.5);
    xlabel('Spanwise Location (m)');
    ylabel('Normalized Deflection');
    title(sprintf('Mode %d Shape (%.2f Hz)', i, freqs(i)));
    grid on;
end
%% 5. SAVE MODEL FOR SUBSEQUENT SCRIPTS
% **MODIFIED SECTION**
% Save the key variables, including the FULL eigenvectors (V_raw)
b = 0.15; h = 0.02;
V_eigenvectors = V_raw; % Assign to a clear variable name
V_deflection = V;       % Assign to a clear variable name
save('structural_model.mat', 'M', 'K', 'V_eigenvectors', 'V_deflection', 'L', 'N', 'b', 'h', 'freqs');
fprintf('\nStructural model data saved successfully to structural_model.mat\n');