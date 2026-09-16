% =========================================================================
% SCRIPT 2: BUILD STRUCTURAL MODEL
% =========================================================================
% Description:
% Constructs a 2D Euler-Bernoulli beam finite element model of the wing.
% It calculates the global mass (M) and stiffness (K) matrices, applies
% fixed-root boundary conditions, and solves the eigenvalue problem to
% find the natural frequencies and mode shapes.
% =========================================================================

% --- Load Parameters
load('wing_geom.mat');
dx = L / N; % Length of one element

% --- Define Local Element Matrices
% Local stiffness matrix for an Euler-Bernoulli beam element
k_local = (E*I/dx^3) * [12, 6*dx, -12, 6*dx;
                       6*dx, 4*dx^2, -6*dx, 2*dx^2;
                       -12, -6*dx, 12, -6*dx;
                       6*dx, 2*dx^2, -6*dx, 4*dx^2];

% Consistent local mass matrix for an Euler-Bernoulli beam element
m_local = (rho*A*dx/420) * [156, 22*dx, 54, -13*dx;
                           22*dx, 4*dx^2, 13*dx, -3*dx^2;
                           54, 13*dx, 156, -22*dx;
                           -13*dx, -3*dx^2, -22*dx, 4*dx^2];

% --- Assemble Global Matrices
n_nodes = N + 1;
K_global = zeros(2 * n_nodes);
M_global = zeros(2 * n_nodes);

for i = 1:N
    dof_indices = (2*i-1):(2*i+2); % DOFs for element i
    K_global(dof_indices, dof_indices) = K_global(dof_indices, dof_indices) + k_local;
    M_global(dof_indices, dof_indices) = M_global(dof_indices, dof_indices) + m_local;
end

% --- Apply Boundary Conditions (Fixed Root at x=0)
% The first two DOFs (deflection and rotation at node 1) are fixed.
dofs_to_keep = 3:(2 * n_nodes);
K = K_global(dofs_to_keep, dofs_to_keep);
M = M_global(dofs_to_keep, dofs_to_keep);

% --- Solve Eigenvalue Problem: K*v = omega^2*M*v
[V_eigenvectors, D] = eig(K, M);
omega = sqrt(diag(D));
freqs_hz = omega / (2 * pi);

% --- Normalize and Extract Mode Shapes
% Normalize eigenvectors for consistent scaling
for i = 1:size(V_eigenvectors, 2)
    V_eigenvectors(:,i) = V_eigenvectors(:,i) / max(abs(V_eigenvectors(:,i)));
end
% Extract translational components for plotting and reconstruction
V_deflection = V_eigenvectors(1:2:end, :);

% --- Plot First 3 Mode Shapes
figure('Name', 'Wing Mode Shapes', 'NumberTitle', 'off');
x_span = linspace(0, L, n_nodes);
for i = 1:3
    subplot(3, 1, i);
    plot(x_span, [0; V_deflection(:,i)], 'o-', 'LineWidth', 1.5, 'MarkerSize', 5);
    title(sprintf('Mode %d Shape (%.2f Hz)', i, freqs_hz(i)));
    ylabel('Norm. Deflection');
    grid on;
    if i == 3, xlabel('Spanwise Location (m)'); end
end

% --- Save Model Data
save('structural_model.mat', 'M', 'K', 'V_eigenvectors', 'V_deflection', 'freqs_hz', 'x_span');