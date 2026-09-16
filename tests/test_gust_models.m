function state = test_gust_models(state)
%TEST_GUST_MODELS  Discrete gust shape, turbulence spectrum, Kussner lag.

fprintf('\n-- gust models --\n');
cfg = aglas_config('baseline');

% ---------------------------------------------------- 1-cosine gust --
U_ds = 10; t_g = 2;
t  = 0:0.001:4;
wg = gust_one_minus_cos(t, U_ds, t_g);

state = aglas_assert(state, 'gust starts at zero', ...
    abs(wg(1)) < 1e-12, sprintf('%.2e', wg(1)));
state = aglas_assert(state, 'gust returns to zero at t_g', ...
    abs(interp1(t, wg, t_g)) < 1e-9, sprintf('%.2e', abs(interp1(t, wg, t_g))));
state = aglas_assert(state, 'gust is zero after t_g', ...
    all(abs(wg(t > t_g)) < 1e-12), '');
state = aglas_assert(state, 'peak equals the design gust velocity', ...
    abs(max(wg) - U_ds) < 1e-9, sprintf('%.6f vs %.1f', max(wg), U_ds));

[~, ip] = max(wg);
state = aglas_assert(state, 'peak occurs at the midpoint t_g/2', ...
    abs(t(ip) - t_g/2) < 2e-3, sprintf('t = %.4f s vs %.4f s', t(ip), t_g/2));

state = aglas_assert(state, 'gust is non-negative throughout', ...
    all(wg >= -1e-12), '');

% A full 1-cosine cycle has zero slope at both ends, so no velocity step is
% injected. The half-cosine form used originally ends at the peak and steps
% discontinuously to zero.
dwg = diff(wg)/ (t(2)-t(1));
slope_end = abs(interp1(t(1:end-1), dwg, t_g - 1e-3));
state = aglas_assert(state, 'gust slope vanishes at the end (no velocity step)', ...
    slope_end < 0.1, sprintf('|dw/dt| = %.4f m/s^2', slope_end));

% Impulse of the full cycle is U_ds*t_g/2, half that of the ramp form.
imp = trapz(t(t <= t_g), wg(t <= t_g));
state = aglas_assert(state, 'gust impulse equals U_ds*t_g/2', ...
    abs(imp - U_ds*t_g/2) < 1e-3, sprintf('%.6f vs %.6f m', imp, U_ds*t_g/2));

% ------------------------------------------------- von Karman spectrum --
[w1, ~, vk] = gust_von_karman(cfg, 400, 0.005);
state = aglas_assert(state, 'turbulence has essentially zero mean', ...
    abs(mean(w1)) < 1e-9, sprintf('%.2e m/s', mean(w1)));
state = aglas_assert(state, 'turbulence RMS is within 10 % of target', ...
    abs(vk.rms_actual - vk.rms_target)/vk.rms_target < 0.10, ...
    sprintf('%.4f vs %.4f m/s', vk.rms_actual, vk.rms_target));

msk = vk.omega > 5 & vk.omega < 60;
slope = polyfit(log(vk.omega(msk)), log(vk.S(msk)), 1);
state = aglas_assert(state, 'PSD follows the -5/3 inertial-subrange slope', ...
    abs(slope(1) + 5/3) < 0.01, sprintf('slope %.4f vs -1.6667', slope(1)));

w2 = gust_von_karman(cfg, 400, 0.005);
state = aglas_assert(state, 'turbulence realisation is reproducible', ...
    isequal(w1, w2), '');

c2 = cfg; c2.turb.seed = cfg.turb.seed + 1;
w3 = gust_von_karman(c2, 400, 0.005);
state = aglas_assert(state, 'a different seed gives a different realisation', ...
    ~isequal(w1, w3), '');

% -------------------------------------------------------- Kussner lag --
fem = beam_fem(cfg);
md  = modal_analysis(fem, 3);
tg  = 0:0.005:8;
wgd = gust_one_minus_cos(tg, cfg.gust.U_ds, cfg.gust.t_g);

[~, qs] = gust_modal_force(cfg, fem, md, tg, wgd, 'quasisteady');
[~, ku] = gust_modal_force(cfg, fem, md, tg, wgd, 'kussner');

state = aglas_assert(state, 'Kussner lag reduces the peak load', ...
    max(ku.q_line) < max(qs.q_line), ...
    sprintf('%.1f vs %.1f N/m', max(ku.q_line), max(qs.q_line)));

[~, i1] = max(qs.q_line); [~, i2] = max(ku.q_line);
state = aglas_assert(state, 'Kussner lag delays the peak', ...
    tg(i2) > tg(i1), sprintf('lag %.3f s', tg(i2) - tg(i1)));

tgf = 0:0.0005:8;
[~, kf] = gust_modal_force(cfg, fem, md, tgf, ...
             gust_one_minus_cos(tgf, cfg.gust.U_ds, cfg.gust.t_g), 'kussner');
d = abs(max(ku.q_line) - max(kf.q_line))/max(kf.q_line);
state = aglas_assert(state, 'Kussner convolution is grid independent', ...
    d < 1e-3, sprintf('%.3e relative difference at 10x resolution', d));

% Quasi-steady load must equal the closed-form strip-theory value.
q_hand = 0.5*cfg.flight.rho_air*cfg.flight.U_inf*cfg.geom.chord * ...
         cfg.aero.cl_alpha_2d*cfg.gust.U_ds;
state = aglas_assert(state, 'quasi-steady load equals 0.5*rho*U*c*cl_alpha*w_g', ...
    abs(max(qs.q_line) - q_hand) < 1e-6*q_hand, ...
    sprintf('%.4f vs %.4f N/m', max(qs.q_line), q_hand));
end
