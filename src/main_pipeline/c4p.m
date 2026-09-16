% =========================================================================
% SCRIPT 4: FLUTTER ANALYSIS
% =========================================================================
% Description:
% Performs a V-g-f (velocity-damping-frequency) flutter analysis using a
% 2D unsteady aerodynamic model (Theodorsen's theory). It calculates the
% eigenvalues of the aeroelastic system for a range of airspeeds and plots
% the results to find the flutter speed and frequency.
% =========================================================================

% --- Load Models
load('wing_geom.mat', 'L', 'b');
load('structural_model.mat', 'M', 'K', 'freqs_hz');

% --- Aeroelastic Parameters
rho_air = 1.225;                % Air density (kg/m^3)
U = 10:1:250;                   % Range of airspeeds to test (m/s)
n_modes = 2;                    % Number of modes for flutter analysis (typically 2-3 is enough)
omega = 2 * pi * freqs_hz(1:n_modes); % Frequencies in rad/s

% --- Build Generalized Aerodynamic Force Matrix (Q)
% Using 2D strip theory and Theodorsen's function.
% This is a frequency-domain model, so we approximate it for the state-space analysis.
% We calculate Q at a representative reduced frequency k = omega*b/U.
% For simplicity, we use k=0.2 for this example.
k_ref = 0.2;
C_k = 0.5; % Theodorsen's function C(k) for k=0.2 (real part)
Q_aero = zeros(n_modes);
for i = 1:n_modes
    for j = 1:n_modes
        % This is a simplified integral of phi_i * L_j over the span
        % Assumes simple mode shapes for this example calculation
        Q_aero(i, j) = (0.5 * rho_air * b^2 * (2*pi)) * (L / (i+j-1));
    end
end

% Extract generalized mass and stiffness for the selected modes
M_gen = eye(n_modes); % Assuming mass-normalized eigenvectors
K_gen = diag(omega.^2);

% --- Perform V-g-f Analysis
eigen_values = zeros(length(U), 2*n_modes);

% *** THIS IS THE CORRECTED LINE ***
for i = 1:length(U)
    % Current airspeed
    U_current = U(i);
    
    % Aerodynamic Damping and Stiffness matrices
    C_aero = (U_current / b) * Q_aero * C_k;
    K_aero = (U_current / b)^2 * Q_aero;
    
    % Assemble the State-Space Matrix A = [0, I; -M^-1*(K+Ka), -M^-1*Ca]
    A = [zeros(n_modes), eye(n_modes);
         -M_gen\(K_gen + K_aero), -M_gen\C_aero];
     
    eigen_values(i, :) = eig(A);
end

% --- Process and Plot Results
damping = real(eigen_values);
frequency = imag(eigen_values) / (2*pi);

% Find flutter point (where damping of any mode becomes positive)
flutter_speed = nan;
flutter_freq = nan;
flutter_mode_idx = -1;

for mode = 1:2*n_modes
    % Find where damping crosses from negative to positive
    cross_idx = find(diff(sign(damping(:, mode))) > 0, 1);
    if ~isempty(cross_idx)
        % Interpolate for a more accurate crossing point
        v1 = U(cross_idx); v2 = U(cross_idx+1);
        d1 = damping(cross_idx, mode); d2 = damping(cross_idx+1, mode);
        flutter_speed_candidate = interp1([d1, d2], [v1, v2], 0);
        
        if isnan(flutter_speed) || flutter_speed_candidate < flutter_speed
            flutter_speed = flutter_speed_candidate;
            flutter_freq = interp1(U, frequency(:, mode), flutter_speed);
            flutter_mode_idx = mode;
        end
    end
end

% --- Create V-g-f Plot
figure('Name', 'Flutter Analysis (V-g-f Plot)', 'NumberTitle', 'off');

% Plot Frequency vs. Airspeed
subplot(2,1,1);
plot(U, frequency, 'LineWidth', 1.5);
hold on;
if ~isnan(flutter_speed)
    plot([flutter_speed flutter_speed], get(gca,'Ylim'), 'r--', 'LineWidth', 1.5);
    text(flutter_speed*1.05, flutter_freq, sprintf('Flutter Speed: %.1f m/s', flutter_speed), 'Color', 'r');
end
title('V-g-f Flutter Analysis');
ylabel('Frequency (Hz)');
grid on; legend('Mode 1', 'Mode 2', 'Location', 'northwest');

% Plot Damping vs. Airspeed
subplot(2,1,2);
plot(U, damping, 'LineWidth', 1.5);
hold on;
plot(get(gca,'Xlim'), [0 0], 'k-'); % Zero damping line
if ~isnan(flutter_speed)
    plot([flutter_speed flutter_speed], get(gca,'Ylim'), 'r--', 'LineWidth', 1.5);
end
xlabel('Airspeed (m/s)');
ylabel('Damping (Real part of Eigenvalue)');
ylim([-max(abs(damping(:)))*0.5, max(abs(damping(:)))*0.1]); % Zoom into the critical region near zero
grid on;

% --- Save Flutter Data
save('aero_data.mat', 'U', 'damping', 'frequency', 'flutter_speed', 'flutter_freq');