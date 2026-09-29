function state = test_control(state)
%TEST_CONTROL  Riccati/Lyapunov solvers and the closed-loop control design.
%
%   Covers the part of the project the rest of the suite does not touch:
%   care_solve, lyap_solve, ctrb_obsv_rank, and the LQG/MRAC/policy-search
%   closed loop against the peak-moment reductions quoted in the README.

fprintf('\n-- control design --\n');

% -------------------------------------------------- Riccati solver --
% Classic double integrator, Q = I, R = 1: the textbook closed-form gain is
% K = [1, sqrt(3)]. This is the check the README's verification table has
% long quoted without it existing in the suite.
A = [0 1; 0 0]; B = [0; 1]; Q = eye(2); R = 1;
[~, K, info] = care_solve(A, B, Q, R);
K_exact = [1, sqrt(3)];
state = aglas_assert(state, 'Riccati solver matches the double-integrator closed form', ...
    norm(K - K_exact) < 1e-6, sprintf('K = [%.6f, %.6f] vs [1, %.6f]', K(1), K(2), sqrt(3)));
state = aglas_assert(state, 'Riccati residual is negligible', ...
    info.residual < 1e-8, sprintf('%.2e', info.residual));

% -------------------------------------------------- Lyapunov solver --
% Scalar case with an exact answer: A = -1, Q = 2 gives -P - P = -2, P = 1.
P_lyap = lyap_solve(-1, 2);
state = aglas_assert(state, 'Lyapunov solver matches the scalar closed form', ...
    abs(P_lyap - 1) < 1e-10, sprintf('%.10f vs 1', P_lyap));

% ------------------------------------------------- plant, controller --
% r_command = 0.3 is the design point used throughout the README and the
% Simulink model (see aglas_sim_setup.m); the config default of 1.0 is a
% much more timid, undertuned controller and is not what "LQG" refers to
% anywhere else in the project.
cfg   = aglas_config('baseline');
cfg.lqr.r_command = 0.3;
fem   = beam_fem(cfg);
modes = modal_analysis(fem, cfg.fem.n_modes);
cs    = control_surface(cfg, fem, modes);
P     = plant_statespace(cfg, fem, modes, cs);
base  = lqg_design(cfg, P);
mrac  = mrac_design(cfg, P, base, struct('gamma', 0.5, 'sigma', 0.05));

[rc, ro] = ctrb_obsv_rank(P.A, P.B, P.C);
state = aglas_assert(state, 'aeroservoelastic plant is fully controllable', ...
    rc == P.n_states, sprintf('rank %d of %d', rc, P.n_states));
state = aglas_assert(state, 'aeroservoelastic plant is fully observable', ...
    ro == P.n_states, sprintf('rank %d of %d', ro, P.n_states));

% ------------------------------------------------------- closed loop --
% Same discrete 1-cosine gust and timing the Simulink comparison and the
% README table use, so these numbers are directly comparable to both.
dt = 2e-4;
t  = (0:dt:6).';
wg = gust_one_minus_cos(t.', cfg.gust.U_ds, cfg.gust.t_g).';

r_open = closed_loop_sim(cfg, P, struct('type', 'none'), t, wg);
r_lqg  = closed_loop_sim(cfg, P, base, t, wg);
r_mrac = closed_loop_sim(cfg, P, mrac, t, wg);

M0    = r_open.peak.moment;
r_lqg_pct  = 100*(1 - r_lqg.peak.moment/M0);
r_mrac_pct = 100*(1 - r_mrac.peak.moment/M0);

state = aglas_assert(state, 'open-loop response is bounded', ...
    isfinite(M0) && M0 > 0, sprintf('%.1f kN m', M0/1e3));
state = aglas_assert(state, 'open-loop peak moment matches the README baseline', ...
    abs(M0/1e3 - 145.2) < 2, sprintf('%.1f kN m vs 145.2', M0/1e3));

state = aglas_assert(state, 'LQG cuts the peak root moment by about half', ...
    abs(r_lqg_pct - 52.9) < 3, sprintf('%.1f %% vs 52.9 %%', r_lqg_pct));
state = aglas_assert(state, 'LQG respects the actuator deflection limit', ...
    r_lqg.peak.cmd_deg <= cfg.control.delta_max + 0.5, ...
    sprintf('%.1f deg vs %.1f deg limit', r_lqg.peak.cmd_deg, cfg.control.delta_max));

state = aglas_assert(state, 'MRAC at the design point is close to, and not better than, LQG', ...
    r_mrac.peak.moment >= r_lqg.peak.moment*0.98 && abs(r_mrac_pct - 50.8) < 4, ...
    sprintf('%.1f %% vs 50.8 %%, LQG %.1f %%', r_mrac_pct, r_lqg_pct));
state = aglas_assert(state, 'MRAC respects the actuator deflection limit', ...
    r_mrac.peak.cmd_deg <= cfg.control.delta_max + 0.5, ...
    sprintf('%.1f deg vs %.1f deg limit', r_mrac.peak.cmd_deg, cfg.control.delta_max));

% ------------------------------------------------------ policy search --
% Starts its search at the LQG solution, so by construction it can never do
% worse: this is the one invariant that must hold regardless of tuning.
pol   = policy_search(cfg, P, base, t, wg, struct('verbose', false, 'max_eval', 80));
r_pol = closed_loop_sim(cfg, P, pol.ctrl, t, wg);
r_pol_pct = 100*(1 - r_pol.peak.moment/M0);

state = aglas_assert(state, 'policy search never does worse than the LQG it starts from', ...
    r_pol.peak.moment <= r_lqg.peak.moment*1.001, ...
    sprintf('%.1f vs %.1f kN m', r_pol.peak.moment/1e3, r_lqg.peak.moment/1e3));
state = aglas_assert(state, 'policy search improves materially on LQG', ...
    abs(r_pol_pct - 59.6) < 5, sprintf('%.1f %% vs 59.6 %%', r_pol_pct));
state = aglas_assert(state, 'policy search respects the actuator deflection limit', ...
    r_pol.peak.cmd_deg <= cfg.control.delta_max + 0.5, ...
    sprintf('%.1f deg vs %.1f deg limit', r_pol.peak.cmd_deg, cfg.control.delta_max));
end
