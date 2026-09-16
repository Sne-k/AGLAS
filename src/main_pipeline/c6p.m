% =========================================================================
% SCRIPT 6: STRESS ANALYSIS
% =========================================================================
% Description:
% Calculates the bending stress in the wing resulting from the gust load.
% It uses the modal stress recovery method, where stress is related to the
% second spatial derivative of the mode shapes (curvature). It finds the
% maximum stress at the wing root and calculates the Factor of Safety (FOS).
% =========================================================================

% --- Load Data
load('wing_geom.mat', 'E', 'I', 'h', 'yield_strength');
load('structural_model.mat', 'V_deflection', 'x_span');
load('modal_response.mat', 'q_sol', 't_sol', 'n_modes');

% --- Modal Stress Recovery Method
% Bending Stress (sigma) = M*y/I, where M = E*I*(d^2(deflection)/dx^2)
% We first find the "modal curvature" (d^2(phi)/dx^2) for each mode shape phi.
d2V_dx2 = zeros(size(V_deflection, 1), n_modes);
V_modes_deflection_full = [zeros(1, n_modes); V_deflection(:, 1:n_modes)]; % Add fixed root

for mode = 1:n_modes
    % Use central finite differences to find the second derivative
    V = V_modes_deflection_full(:, mode);
    for node = 2:length(x_span)-1
        d2V_dx2(node-1, mode) = (V(node+1) - 2*V(node) + V(node-1)) / (x_span(2)-x_span(1))^2;
    end
end

% Modal Bending Moment: M_modal = E*I * (modal curvature)
M_modal = E * I * d2V_dx2;

% Reconstruct time history of bending moment at each node
% M(x,t) = sum over modes [ q(t) * M_modal(x) ]
M_t = q_sol * M_modal'; % Result is a [time x nodes] matrix

% --- Calculate Bending Stress
% sigma = M*c/I, where c is the distance from the neutral axis (h/2)
c = h / 2;
stress_t = (M_t * c) / I;

% --- Analyze Stress at the Wing Root (Node 2, index 1 in our matrices)
stress_root_t = stress_t(:, 1);
[max_stress_root, time_idx] = max(abs(stress_root_t));
time_of_max_stress = t_sol(time_idx);
stress_dist_at_max = stress_t(time_idx, :);

% --- Calculate Factor of Safety (FOS)
FOS = yield_strength / max_stress_root;

% --- Plot Results
figure('Name', 'Stress Analysis Results', 'NumberTitle', 'off');

% Plot 1: Stress at Wing Root vs. Time
subplot(2,1,1);
plot(t_sol, stress_root_t / 1e6, 'b', 'LineWidth', 1.5);
hold on;
plot(time_of_max_stress, max_stress_root / 1e6, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('Bending Stress at Wing Root');
xlabel('Time (s)');
ylabel('Stress (MPa)');
legend(sprintf('Max Stress: %.2f MPa', max_stress_root / 1e6), 'Location', 'southeast');
grid on;

% Plot 2: Stress Distribution along Span at Time of Max Stress
subplot(2,1,2);
plot(x_span(2:end), stress_dist_at_max / 1e6, 'r', 'LineWidth', 1.5);
title(sprintf('Stress Distribution at t = %.2f s (Time of Peak Stress)', time_of_max_stress));
xlabel('Spanwise Location (m)');
ylabel('Stress (MPa)');
grid on;

% --- Display Final Results in Command Window
fprintf('--- Stress Analysis Summary ---\n');
fprintf('Max Bending Stress at Root: %.2f MPa\n', max_stress_root / 1e6);
fprintf('Yield Strength of Material: %.0f MPa\n', yield_strength / 1e6);
fprintf('Factor of Safety (FOS): %.2f\n', FOS);
if FOS < 1.5
    fprintf('WARNING: Factor of Safety is below the typical minimum of 1.5.\n');
end

% --- Save Results
save('fos_data.mat', 'max_stress_root', 'FOS', 'yield_strength');