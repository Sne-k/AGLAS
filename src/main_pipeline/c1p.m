% =========================================================================
% AGLAS STAGE 1 (main pipeline) : WING DEFINITION
% =========================================================================
% Resolves the wing geometry, material and section properties from the
% central configuration and writes them to data/wing_geom.mat for the rest
% of the pipeline.
%
% Change parameters in src/common/aglas_config.m, not here. Nothing in this
% file hard-codes a physical number.
%
% Outputs : data/wing_geom.mat
% =========================================================================

clear; close all; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'common'));
paths = aglas_paths();

PRESET = 'baseline';          % 'baseline' or 'legacy'

cfg = aglas_config(PRESET);
g   = cfg.geom;
m   = cfg.material;
s   = cfg.section;

fprintf('=== AGLAS Stage 1: wing definition (preset "%s") ===\n\n', PRESET);

fprintf('Planform\n');
fprintf('  semi-span L          : %8.3f m\n',      g.L);
fprintf('  chord                : %8.3f m\n',      g.chord);
fprintf('  aspect ratio (full)  : %8.2f\n',        cfg.aero.AR);
fprintf('  elastic axis         : %8.1f %% chord\n', 100*g.ea_frac);
fprintf('  section CG           : %8.1f %% chord\n', 100*g.cg_frac);

fprintf('\nMaterial: %s\n', m.name);
fprintf('  E                    : %8.3f GPa\n',    m.E/1e9);
fprintf('  G                    : %8.3f GPa\n',    m.G/1e9);
fprintf('  density              : %8.1f kg/m^3\n', m.rho);
fprintf('  yield strength       : %8.1f MPa\n',    m.sigma_yield/1e6);

fprintf('\nSection (%s)\n', g.section);
fprintf('  area                 : %10.6f m^2\n',   s.A);
fprintf('  I (chordwise axis)   : %10.4e m^4\n',   s.I);
fprintf('  J (torsion constant) : %10.4e m^4\n',   s.J);
fprintf('  outer fibre distance : %10.6f m\n',     s.c_outer);
fprintf('  EI                   : %10.4e N m^2\n', s.EI);
fprintf('  GJ                   : %10.4e N m^2\n', s.GJ);
fprintf('  GJ/EI                : %10.3f\n',       s.GJ/s.EI);
fprintf('  mass per span        : %10.3f kg/m  (structural %.3f + other %.3f)\n', ...
        s.m_total, s.m_struct, g.m_nonstruct);
fprintf('  total semi-span mass : %10.1f kg\n',    s.m_total*g.L);
fprintf('  polar inertia/span   : %10.4f kg m^2/m\n', s.I_theta);

fprintf('\nDiscretisation\n');
fprintf('  elements             : %8d\n',  cfg.fem.n_elements);
fprintf('  torsion DOF included : %8d\n',  cfg.fem.include_torsion);
fprintf('  modes retained       : %8d\n',  cfg.fem.n_modes);

out_file = fullfile(paths.data, 'wing_geom.mat');
save(out_file, 'cfg');
fprintf('\nSaved %s\n', out_file);
