% =========================================================================
% SCRIPT 1: SETUP WING PARAMETERS
% =========================================================================
% Description:
% Defines the geometric, material, and discretization parameters for the
% wing model. These parameters are saved to 'wing_geom.mat' for use in
% subsequent scripts.
% =========================================================================

% --- Wing Geometry
L = 10;     % Length of the wing (m)
b = 0.15;   % Cross-section width (base) (m)
h = 0.02;   % Cross-section height (m)

% --- Material Properties (Aluminum)
E   = 70e9;   % Young's Modulus (Pa)
rho = 2700;   % Density (kg/m^3)
yield_strength = 276e6; % Yield strength for Aluminum (Pa)

% --- FEM Discretization
N = 10;     % Number of finite elements

% --- Calculated Properties
A = b * h;          % Cross-sectional area (m^2)
I = b * h^3 / 12;   % Area moment of inertia for rectangular section (m^4)

% --- Save data to file
save('wing_geom.mat', 'L', 'b', 'h', 'E', 'rho', 'yield_strength', 'N', 'A', 'I');