function cfg = aglas_config(preset)
%AGLAS_CONFIG  Single source of truth for every AGLAS model parameter.
%
%   cfg = AGLAS_CONFIG()          returns the default ('baseline') preset.
%   cfg = AGLAS_CONFIG(preset)    returns a named preset.
%
%   Presets
%   -------
%     'baseline' : Physically consistent 10 m semi-span aluminium wing with a
%                  thin-walled box spar. This is the model the analyses are
%                  validated against and the default for all scripts.
%     'legacy'   : The original solid 0.15 x 0.02 m rectangular section, kept
%                  so the pre-existing results in data/ remain reproducible.
%                  It is a flat plank, not a wing: first bending mode is
%                  0.16 Hz and it deflects ~19 % of span under load, which is
%                  far outside the small-deflection range Euler-Bernoulli beam
%                  theory is valid in. Use for regression checks only.
%
%   Every downstream script reads its numbers from here. Do not hard-code a
%   parameter anywhere else; that is what produced the L = 5 m vs L = 10 m and
%   I = 5e-6 vs 1e-7 contradictions in the original code.
%
%   See also SECTION_PROPERTIES, BEAM_FEM.

    if nargin < 1 || isempty(preset)
        preset = 'baseline';
    end
    preset = lower(strtrim(preset));

    % ---------------------------------------------------------------- common
    cfg = struct();
    cfg.preset = preset;

    % --- Discretisation
    cfg.fem.n_elements   = 20;     % beam elements along the semi-span
    cfg.fem.n_modes      = 3;      % modes retained for modal superposition
    cfg.fem.include_torsion = true;% adds a twist DOF per node (3 DOF/node)

    % --- Atmosphere / flight condition (ISA sea level)
    cfg.flight.rho_air   = 1.225;  % air density                        [kg/m^3]
    cfg.flight.U_inf     = 50;     % reference true airspeed            [m/s]
    cfg.flight.a_sound   = 340.3;  % speed of sound                     [m/s]

    % --- Discrete gust ("1-cosine"), CS-25.341(a) form
    cfg.gust.U_ds        = 10;     % peak gust velocity                 [m/s]
    cfg.gust.t_g         = 2.0;    % gust duration (full 1-cos cycle)   [s]

    % --- Continuous turbulence (von Karman)
    cfg.turb.sigma       = 3.0;    % turbulence RMS velocity            [m/s]
    cfg.turb.L_scale     = 762;    % von Karman length scale (MIL-8785) [m]
    cfg.turb.seed        = 20250907; % fixed RNG seed => reproducible

    % --- Structural damping (modal, fraction of critical)
    cfg.damping.zeta     = 0.02;   % 2 % of critical in every retained mode

    % --- Time integration (ode45)
    cfg.time.t_end       = 8.0;    % simulation duration                [s]
    cfg.time.dt_out      = 0.005;  % output sample interval             [s]
    cfg.time.rel_tol     = 1e-8;
    cfg.time.abs_tol     = 1e-10;

    % --- Flutter sweep
    cfg.flutter.U        = 10:2:400;  % airspeed sweep                  [m/s]
    cfg.flutter.n_modes  = 6;         % modes carried into the p-k solution
    cfg.flutter.k_tol    = 1e-8;      % reduced-frequency convergence
    cfg.flutter.max_iter = 80;

    % --- Safety / optimisation thresholds
    cfg.safety.fos_critical = 1.5;  % below this: red
    cfg.safety.fos_target   = 2.5;  % aerospace acceptance level
    cfg.safety.fos_max_goal = 3.0;  % thickness optimiser aims for this peak
    cfg.safety.fos_plot_cap = 10;   % clip for readable colour maps

    % --- Model-validity guards
    cfg.limits.max_tip_defl_frac = 0.10;  % small-deflection limit, fraction of span
    cfg.limits.max_delta_alpha   = 10;    % linear-aero limit, degrees

    % ------------------------------------------------------------- presets --
    switch preset
        case 'baseline'
            cfg.geom.L       = 10.0;   % semi-span                      [m]
            cfg.geom.chord   = 1.5;    % chord (rectangular planform)   [m]
            cfg.geom.section = 'box';
            cfg.geom.box_w   = 0.525;  % spar box width,  35 % of chord [m]
            cfg.geom.box_h   = 0.180;  % spar box height, 12 % of chord [m]
            cfg.geom.wall_t  = 0.009;  % uniform wall thickness         [m]
            cfg.geom.m_nonstruct = 4.8;% skin/fuel/systems mass per span[kg/m]
            cfg.geom.ea_frac = 0.35;   % elastic axis, fraction of chord from LE
            cfg.geom.cg_frac = 0.40;   % section CG,  fraction of chord from LE

            cfg.material.name  = 'Aluminium 2024-T3';
            cfg.material.E     = 73.1e9;
            cfg.material.nu    = 0.33;
            cfg.material.rho   = 2780;
            cfg.material.sigma_yield = 345e6;

        case 'legacy'
            % Original repository parameters, preserved verbatim.
            cfg.geom.L       = 10.0;
            cfg.geom.chord   = 1.0;    % never defined in the original; the old
                                       % c4p.m wrongly reused the 0.15 m section
                                       % width as an aerodynamic semi-chord.
            cfg.geom.section = 'solid';
            cfg.geom.box_w   = 0.15;   % solid section width
            cfg.geom.box_h   = 0.02;   % solid section height
            cfg.geom.wall_t  = NaN;
            cfg.geom.m_nonstruct = 0.0;
            cfg.geom.ea_frac = 0.35;
            cfg.geom.cg_frac = 0.40;

            cfg.material.name  = 'Aluminium (generic)';
            cfg.material.E     = 70e9;
            cfg.material.nu    = 0.33;
            cfg.material.rho   = 2700;
            cfg.material.sigma_yield = 276e6;

            cfg.fem.n_elements = 10;   % as originally used
            cfg.fem.include_torsion = false;

        otherwise
            error('aglas_config:unknownPreset', ...
                  'Unknown preset "%s". Use ''baseline'' or ''legacy''.', preset);
    end

    % --- Derived material constants
    cfg.material.G = cfg.material.E / (2 * (1 + cfg.material.nu));

    % --- Derived aerodynamics
    cfg.aero.semi_chord = cfg.geom.chord / 2;
    % Lift-curve slope corrected for finite aspect ratio (Helmbold, incompressible).
    AR = 2 * cfg.geom.L / cfg.geom.chord;          % full-wing aspect ratio
    cfg.aero.AR      = AR;
    cfg.aero.cl_alpha_2d = 2 * pi;                  % thin-airfoil value  [1/rad]
    cfg.aero.CL_alpha    = 2*pi*AR / (2 + sqrt(AR^2 + 4));
    % Distance from mid-chord to elastic axis, non-dimensional on semi-chord.
    % Positive aft. Theodorsen's 'a'.
    cfg.aero.a_ea  = (cfg.geom.ea_frac - 0.5) / 0.5;
    % Distance from elastic axis to CG, non-dimensional on semi-chord ('x_alpha').
    cfg.aero.x_alpha = (cfg.geom.cg_frac - cfg.geom.ea_frac) / 0.5;

    % --- Derived section properties
    cfg.section = section_properties(cfg);
end
