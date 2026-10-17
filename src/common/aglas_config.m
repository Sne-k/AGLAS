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

    % --- Aerodynamics of the control plant (plant_statespace)
    % 'quasisteady' is the k -> 0 limit of Theodorsen strip theory: no
    % circulation lag, no apparent mass, gust applied without Kussner lag.
    % 'unsteady' adds a two-state indicial lag on every circulatory force
    % (wing motion, flap and gust) and the apparent-mass terms, so the plant
    % carries the C(k) phase lag that drives flutter.
    cfg.aero.plant_model = 'quasisteady';
    % Indicial lift function for the circulatory lag, as A_i, b_i in
    % phi(s) = 1 - A_1*exp(-b_1*beta^2*s) - A_2*exp(-b_2*beta^2*s), s in
    % semi-chords. 'jones' is R.T. Jones' fit to the Wagner function and
    % reproduces Theodorsen's C(k) to a few percent. 'leishman' is the
    % compressible fit of Leishman & Nguyen (AIAA J 28(5), 1990).
    cfg.aero.indicial    = 'jones';
    % Prandtl-Glauert compressibility in the unsteady plant: lift slope
    % cl_alpha/beta and lag rates scaled by beta^2, beta = sqrt(1 - M^2).
    % Linear subsonic theory; not valid approaching M = 1.
    cfg.aero.compressible = false;

    % --- Discrete gust ("1-cosine"), CS-25.341(a) form
    cfg.gust.U_ds        = 10;     % peak gust velocity                 [m/s]
    cfg.gust.t_g         = 2.0;    % gust duration (full 1-cos cycle)   [s]

    % --- Continuous turbulence (von Karman)
    cfg.turb.sigma       = 3.0;    % turbulence RMS velocity            [m/s]
    cfg.turb.L_scale     = 762;    % von Karman length scale (MIL-8785) [m]
    cfg.turb.seed        = 20250907; % fixed RNG seed => reproducible
    cfg.turb.T_spectrum  = 200;    % record for the spectrum check: long enough
                                   % for a converged PSD estimate          [s]
    cfg.turb.T_response  = 30;     % record for the wing response: stable RMS,
                                   % short enough to integrate quickly     [s]
    cfg.turb.dt_kussner  = 0.02;   % Kussner recurrence step for turbulence;
                                   % the discrete gust uses t_g/20, which is
                                   % far too coarse for broadband input   [s]

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

    % --- Gust load alleviation control
    % Trailing-edge control surface. The span limits are snapped to element
    % boundaries by control_surface so the consistent load vector stays exact.
    cfg.control.flap_chord_frac = 0.25;  % flap chord / wing chord      [-]
    cfg.control.span_start_frac = 0.60;  % inboard edge, fraction of semi-span
    cfg.control.span_end_frac   = 0.95;  % outboard edge
    cfg.control.delta_max       = 20;    % deflection limit             [deg]
    cfg.control.rate_max        = 150;   % rate limit                   [deg/s]
    cfg.control.act_freq_hz     = 15;    % actuator natural frequency   [Hz]
    cfg.control.act_zeta        = 0.70;  % actuator damping ratio       [-]

    % Sensors used by the estimator.
    cfg.sensor.strain_noise     = 2e-6;  % root strain noise, RMS       [-]
    cfg.sensor.accel_noise      = 0.02;  % tip accelerometer noise, RMS [m/s^2]
    % Add that noise to the measurements in closed_loop_sim, white and
    % independent each integration step at the RMS levels above. Off by
    % default so every documented result stays reproducible as measured.
    cfg.sensor.simulate_noise   = false;
    cfg.sensor.probe_noise      = 0.1;   % gust probe noise, RMS        [m/s]
    cfg.sensor.noise_seed       = 7;
    % Process noise intensity on the gust state. This sets the estimator
    % bandwidth: too high and the filter tries to differentiate its sensors,
    % producing an observer far faster than the structure and impossible to
    % integrate. 0.05 places the fastest observer pole near 150 Hz, about seven
    % times the highest retained mode.
    cfg.sensor.gust_psd         = 0.05;

    % LQR cost weights. R is swept during tuning; these are the defaults the
    % design scripts start from. See docs/control-design-log.md section 4 for
    % the trade study. Every control result quoted in the README uses
    % r_command = 0.3 (set by aglas_sim_setup), not this default of 1.0,
    % which gives a far more timid controller (about 26 % load reduction).
    cfg.lqr.q_moment            = 1.0;   % weight on root bending moment
    cfg.lqr.q_tip               = 0.0;   % weight on tip deflection
    cfg.lqr.r_command           = 1.0;   % weight on commanded deflection

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
    cfg.aero.AR      = 2 * cfg.geom.L / cfg.geom.chord;   % full-wing aspect ratio
    cfg.aero.cl_alpha_2d = 2 * pi;                  % thin-airfoil value  [1/rad]
    % Distance from mid-chord to elastic axis, non-dimensional on semi-chord.
    % Positive aft. Theodorsen's 'a'.
    cfg.aero.a_ea  = (cfg.geom.ea_frac - 0.5) / 0.5;

    % --- Derived section properties
    cfg.section = section_properties(cfg);
end
