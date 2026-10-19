function state = test_feedforward(state)
%TEST_FEEDFORWARD  Probe-referenced gust feedforward, fixed and adaptive.

fprintf('\n-- gust feedforward --\n');

cfg = aglas_config('baseline');
cfg.lqr.r_command = 0.3;
fem  = beam_fem(cfg);
md   = modal_analysis(fem, cfg.fem.n_modes);
P    = plant_statespace(cfg, fem, md, control_surface(cfg, fem, md));
base = lqg_design(cfg, P);

dt = 2e-4;
tg = (0:dt:6).';
wg = gust_one_minus_cos(tg.', cfg.gust.U_ds, cfg.gust.t_g).';

% A feedforward path with zero weights must not change anything.
ff0 = feedforward_design(cfg, P, base);
c = base; c.ff = ff0;
r_l = closed_loop_sim(cfg, P, base, tg, wg);
r_0 = closed_loop_sim(cfg, P, c,    tg, wg);
state = aglas_assert(state, 'zero feedforward leaves the LQG response unchanged', ...
    isequal(r_l.moment, r_0.moment), '');

% An explicit probe signal equal to what the built-in probe sees must
% reproduce it exactly.
ffw0 = ff0; ffw0.w = 1e-3 * ones(size(ff0.w));
tau  = ffw0.probe_ahead / cfg.flight.U_inf;
fr   = ffw0; fr.ref = interp1(tg, wg, tg + tau, 'linear', wg(end)).';
ca = base; ca.ff = ffw0;  cb = base; cb.ff = fr;
ra = closed_loop_sim(cfg, P, ca, tg, wg);
rb = closed_loop_sim(cfg, P, cb, tg, wg);
state = aglas_assert(state, 'explicit probe signal reproduces the built-in probe', ...
    max(abs(ra.moment - rb.moment)) < 1e-9*max(abs(ra.moment)), '');

% Probe model: with no decay it is the frozen-turbulence advance; with decay
% the correlation with the advanced gust falls as the sensor moves further
% ahead, and the measured signal keeps the spectrum's variance.
[w1, t1] = gust_von_karman(cfg, 60, 0.005);
c3 = cfg; c3.turb.seed = cfg.turb.seed + 2;
w2 = gust_von_karman(c3, 60, 0.005);
t1 = t1(:); w1 = w1(:);
r0 = probe_signal(cfg, t1, w1, 50);
adv = interp1(t1, w1, t1 + 50/cfg.flight.U_inf, 'linear', w1(end));
state = aglas_assert(state, 'probe with no decay is the frozen-turbulence advance', ...
    isequal(r0, adv), '');
rho = zeros(1, 3);  ds = [10 50 200];
for i = 1:3
    ri = probe_signal(cfg, t1, w1, ds(i), w2, 4);
    ai = interp1(t1, w1, t1 + ds(i)/cfg.flight.U_inf, 'linear', w1(end));
    cm = corrcoef(ri, ai);  rho(i) = cm(1,2);
end
state = aglas_assert(state, 'probe coherence falls with look-ahead distance', ...
    all(diff(rho) < 0) && rho(1) < 1, sprintf('correlation %.3f, %.3f, %.3f at 10, 50, 200 m', rho));

% Variance is preserved in expectation only: the turbulence generator fixes
% the spectral amplitudes and randomises the phases, so one realisation
% carries a random cross-term dominated by a few low-frequency bins.
% Average over independent fields.
vr = zeros(1, 20);
for s = 1:20
    cs = cfg; cs.turb.seed = cfg.turb.seed + 100 + s;
    wi = gust_von_karman(cs, 60, 0.005);
    ri = probe_signal(cfg, t1, w1, 200, wi, 4);
    vr(s) = var(ri) / var(interp1(t1, w1, t1 + 200/cfg.flight.U_inf, 'linear', w1(end)));
end
state = aglas_assert(state, 'probe signal keeps the turbulence variance on average', ...
    abs(mean(vr) - 1) < 0.1, sprintf('mean variance ratio %.3f over 20 fields', mean(vr)));

% A command delay must reach the plant exactly nd steps late.
cd_ = cfg; cd_.control.delay = 0.010;
rd  = closed_loop_sim(cd_, P, base, tg, wg);
nd  = round(cd_.control.delay / (tg(2) - tg(1)));
state = aglas_assert(state, 'command delay reaches the plant exactly that late', ...
    isequal(rd.u_applied(nd+1:end-1), rd.delta_cmd(1:end-nd-1)), sprintf('%d steps', nd));

% The secondary-path model must predict the strain a feedforward command
% produces. Check its static gain against the closed-loop steady state.
K = base.K; L = base.L; G = base.Be - L*base.De + L*P.D;
Acl = [P.A, -P.B*K; L*P.C, base.Ae - L*base.Ce - G*K];
Bcl = [P.B; G];
dc  = [P.C(1,:), -P.D(1)*K] * (-Acl \ Bcl);
ffl = feedforward_design(cfg, P, base, struct('S_len', 15));   % long enough to settle
state = aglas_assert(state, 'secondary-path model has the closed-loop static gain', ...
    abs(sum(ffl.S) - dc) < 1e-3*abs(dc), sprintf('%.6e vs %.6e', sum(ffl.S), dc));

% Wiener feedforward from a probe 10 m ahead, trained and tested on
% different turbulence records.
T = 40;
[wa, ta] = gust_von_karman(cfg, T, cfg.time.dt_out);
c2 = cfg; c2.turb.seed = cfg.turb.seed + 1;
[wb, tb] = gust_von_karman(c2, T, cfg.time.dt_out);
t  = (0:dt:T-1).';
wA = interp1(ta(:), wa(:), t, 'linear', 0);
wB = interp1(tb(:), wb(:), t, 'linear', 0);

ffw = feedforward_design(cfg, P, base, struct('train', {{t, wA}}));
c = base; c.ff = ffw;
rl = closed_loop_sim(cfg, P, base, t, wB);
rf = closed_loop_sim(cfg, P, c,    t, wB);
gain = 1 - rms(rf.moment)/rms(rl.moment);
state = aglas_assert(state, 'Wiener feedforward cuts turbulence load RMS beyond LQG', ...
    gain > 0.3, sprintf('%.1f %% further reduction on an unseen record', 100*gain));
state = aglas_assert(state, 'feedforward respects the actuator deflection limit', ...
    rf.peak.cmd_deg <= cfg.control.delta_max + 1e-9, ...
    sprintf('%.2f deg vs %.1f deg limit', rf.peak.cmd_deg, cfg.control.delta_max));

% Filtered-x LMS starting from zero must learn a useful filter.
fa = ffw; fa.w = zeros(size(ffw.w)); fa.adapt = true; fa.mu = 0.003;
c = base; c.ff = fa;
ra = closed_loop_sim(cfg, P, c, t, wB);
late = t >= T - 15;
ga = 1 - rms(ra.moment(late))/rms(rl.moment(late));
state = aglas_assert(state, 'adaptive feedforward learns to reduce the load', ...
    ga > 0.1, sprintf('%.1f %% over the last 15 s', 100*ga));
end
