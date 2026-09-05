% Load the FoS data
load('fos_data.mat');  % Make sure the variable is named FoS_matrix

% Create binary mask: 1 where FoS < 2.5, else 0
FoS_mask = FoS_matrix < 2.5;

% Plot the binary mask
binary_mask = FoS_matrix >= 2.5;  % 1 = Safe, 0 = Unsafe

figure;
imagesc(binary_mask);
colormap([1 0 0; 0 1 0]);  % 0 = Red (Unsafe), 1 = Green (Safe)
colorbar('Ticks', [0,1], 'TickLabels', {'Unsafe', 'Safe'});
title('Binary Mask for Unsafe FoS Regions (FoS < 2.5)');
xlabel('Spanwise Elements');
ylabel('Chordwise Elements');
%% 

load('fos_data.mat');  % Assuming FoS_matrix is loaded

% Target FoS
target_max_fos = 3.0;

% Initialize
thickness_factor = 1.0;
step = 0.01;
min_fos_margin = 0.05;  % margin from the target to avoid dipping too low

% Optimization loop
while true
    scaled_FoS = FoS_matrix * thickness_factor;
    current_max_fos = max(scaled_FoS(:));
    
    if current_max_fos <= target_max_fos + min_fos_margin
        break;
    end
    
    thickness_factor = thickness_factor - step;
    if thickness_factor <= 0.1
        warning('Stopped early to avoid unrealistic reduction');
        break;
    end
end

% Display final thickness factor
fprintf('Optimal thickness scaling factor: %.3f\n', thickness_factor);
fprintf('Max FoS after scaling: %.3f\n', max(scaled_FoS(:)));

% Plot optimized FoS
figure;
imagesc(scaled_FoS);
colorbar;
title('Optimized Factor of Safety Distribution');
xlabel('Spanwise Elements'); ylabel('Chordwise Elements');
