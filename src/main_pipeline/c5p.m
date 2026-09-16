% =========================================================================
% SCRIPT 5: ANIMATE GUST RESPONSE
% =========================================================================
% Description:
% Creates an animation of the wing's deflection over time based on the
% results of the gust response simulation. The animation can optionally
% be saved as a GIF file.
% =========================================================================

% --- Configuration
SAVE_GIF = true; % Set to true to save animation, false to just display
filename = 'wing_gust_animation.gif'; % Name for the GIF file
frame_delay = 0.05; % Time between animation frames

% --- Load Data
load('modal_response.mat'); % Loads q_sol, t_sol, n_modes
load('structural_model.mat', 'V_deflection', 'x_span');

% --- Setup Figure
fig = figure('Name', 'Wing Gust Response Animation', 'NumberTitle', 'off');
ax = gca;
ax.NextPlot = 'replaceChildren';

% Calculate dynamic plot limits
V_modes_deflection = V_deflection(:, 1:n_modes);
full_deflection = q_sol * V_modes_deflection';
max_def = max(abs(full_deflection(:))) * 1.2; % 20% margin
if max_def == 0, max_def = 1; end % Avoid zero limits if no deflection

% --- Animation Loop
for i = 1:5:length(t_sol)
    % Reconstruct the wing's deflected shape at the current time step
    q_instantaneous = q_sol(i, :);
    deflection = q_instantaneous * V_modes_deflection';
    
    % Plot the wing shape
    plot(x_span, [0, deflection], 'b-', 'LineWidth', 2);
    hold on;
    plot(x_span, zeros(size(x_span)), 'k--'); % Show undeflected position
    hold off;
    
    % Format plot
    grid on;
    ylim([-max_def, max_def]);
    xlim([x_span(1), x_span(end)]);
    xlabel('Spanwise Location (m)');
    ylabel('Deflection (m)');
    title(sprintf('Wing Deflection at Time: %.2f s', t_sol(i)));
    
    drawnow;
    
    % Capture frame for GIF
    if SAVE_GIF
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', frame_delay);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', frame_delay);
        end
    else
        pause(frame_delay);
    end
end