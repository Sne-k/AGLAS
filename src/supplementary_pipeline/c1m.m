clc; clear; close all;

% Parameters
V_max = 10; % Maximum gust velocity (m/s)
t_g = 2;   % Gust penetration time (seconds)
T = 5;     % Total simulation time (seconds)
ts = 0.01; % Time step

t = 0:ts:T;

% 1-Cosine Gust Model
V_g_cosine = (V_max / 2) * (1 - cos(pi * t / t_g));
V_g_cosine(t > t_g) = 0; % Gust ends after t_g

% Von Kármán Turbulence Model (Random gusts)
N = length(t);
white_noise = randn(N,1); % Generate white noise
fc = 1 / t_g; % Cutoff frequency based on gust penetration time
[b, a] = butter(2, fc * 2 * ts, 'low'); % Low-pass filter
V_g_turbulence = filter(b, a, white_noise) * V_max; % Scale by max gust velocity

% Plot results
figure;
subplot(2,1,1);
plot(t, V_g_cosine, 'b', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Gust Velocity (m/s)');
title('1-Cosine Gust Model'); grid on;

subplot(2,1,2);
plot(t, V_g_turbulence, 'r', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Gust Velocity (m/s)');
title('Von Kármán Turbulence Model'); grid on;
