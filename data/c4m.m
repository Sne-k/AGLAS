clc;
clear;
close all;

%% Wing and Beam Properties
E = 70e9;            % Young's modulus [Pa]
I = 5e-6;            % Moment of inertia [m^4]
z = 0.05;            % Distance from neutral axis [m]
L = 5;               % Wing span [m]
N = 50;              % Number of spatial nodes
x = linspace(0, L, N);
dx = x(2) - x(1);

%% Time and Mode Setup
T = 5;                       % Total time [s]
dt = 0.01;
t = 0:dt:T;
n_t = length(t);
n_modes = 3;

%% Generate Synthetic Mode Shapes (Fixed-Free Beam)
Phi = zeros(N, n_modes);
for i = 1:n_modes
    Phi(:, i) = sin((i * pi * x) / (2 * L)); % Simple sinusoidal shapes
end

%% Generate Synthetic Modal Response (e.g., damped response to gust)
q = zeros(n_modes, n_t);
for i = 1:n_modes
    q(i, :) = 0.02 * sin(2 * pi * i * t) .* exp(-0.3 * t);  % Damped sine
end

%% Compute Displacement Field w(x,t)
w = Phi * q;   % Reconstruct full displacement

%% Compute Stress Over Time
stress = zeros(N, n_t);
for i = 1:n_t
    w_i = w(:, i);
    curvature = gradient(gradient(w_i, dx), dx);  % ∂²w/∂x²
    stress(:, i) = -E * z * curvature;            % σ = -Ez ∂²w/∂x²
end

%% Plot Wing Tip Stress Over Time
figure;
plot(t, stress(end, :), 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Stress (Pa)');
title('Wing Tip Bending Stress Over Time');
grid on;

%% Optional: Plot Stress Distribution Along Span at Final Time
figure;
plot(x, stress(:, end), 'b', 'LineWidth', 2);
xlabel('Spanwise Location (m)');
ylabel('Stress (Pa)');
title('Stress Distribution at Final Time');
grid on;

%% Optional: 3D Surface Plot
figure;
surf(t, x, stress, 'EdgeColor', 'none');
xlabel('Time (s)');
ylabel('Span Position (m)');
zlabel('Stress (Pa)');
title('Bending Stress Distribution Over Wing Span and Time');
colorbar;
view(45, 30);
%% 
% --- INPUTS (reuse from your previous script) ---
L = 5;              % Wing span (m)
N = 50;             % Number of spanwise nodes
x = linspace(0, L, N);     % Spanwise positions
b = 0.15;           % Width of wing section (m)
h = 0.02;           % Height of wing section (m)
I = (b * h^3) / 12; % Moment of inertia
E = 70e9;           % Young’s modulus (Pa)
c = h / 2;          % Distance to outer fiber

% --- Load or define q and Phi ---
% q (n_modes x time), Phi (N x n_modes), t (1 x time)

% --- Compute stress ---
M = Phi * q;                         % Bending moment at each span point over time (N x time)
stress_matrix = -c * M / I;         % Stress (Pa) (N x time)

% --- Compute strain ---
strain_matrix = stress_matrix / E;  % Strain (unitless)

% --- Plot Strain over Time at Wing Tip ---
figure;
plot(t, strain_matrix(end, :), 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Strain');
title('Wing Tip Strain Over Time');
max(abs(strain_matrix(:, end)))

% --- Plot Spanwise Strain at Final Time ---
figure;
plot(x, strain_matrix(:, end), 'b', 'LineWidth', 1.5);
xlabel('Spanwise Location (m)');
ylabel('Strain');
title('Strain Distribution at Final Time');

% --- 3D Surface Plot of Strain Over Time and Span ---
figure;
surf(t, x, strain_matrix, 'EdgeColor', 'none');
xlabel('Time (s)');
ylabel('Span Position (m)');
zlabel('Strain');
title('Strain Distribution Over Wing Span and Time');
colorbar;
view(135, 30);

% === Load or use your stress matrix ===
% Assume `stress_matrix` is already in workspace

% === Material Property ===
sigma_yield = 345e6; % Pa, for Aluminum 2024-T3

% === Calculate FoS ===
FoS_matrix = sigma_yield ./ abs(stress_matrix);  % Element-wise

% === Clip extreme values for plotting ===
FoS_matrix(FoS_matrix > 10) = 10;

% === Plotting ===
figure;
imagesc(FoS_matrix);
colorbar;
title('Factor of Safety Distribution');
xlabel('Spanwise Elements');
ylabel('Chordwise Elements');
colormap jet;
% Example FoS matrix (replace with your real FoS data)
FoS_matrix = rand(50, 500)*3 + 8;  % Dummy safe data: FoS ~ [8–11]

% Define custom colormap
custom_map = [1 0 0;    % Red for FoS < 1.5
              1 1 0;    % Yellow for 1.5 <= FoS < 2.5
              0 1 0];   % Green for FoS >= 2.5

% Classify FoS regions
FoS_classified = zeros(size(FoS_matrix));
FoS_classified(FoS_matrix < 1.5) = 1;       % Red
FoS_classified(FoS_matrix < 2.5 & FoS_matrix >= 1.5) = 2;  % Yellow
FoS_classified(FoS_matrix >= 2.5) = 3;      % Green

% Plot with warning zones
figure;
imagesc(FoS_classified);
colormap(custom_map);
colorbar('Ticks', [1 2 3], 'TickLabels', {'<1.5 (Critical)', '1.5-2.5 (Caution)', '≥2.5 (Safe)'});
xlabel('Spanwise Elements');
ylabel('Chordwise Elements');
title('Factor of Safety Zones');
