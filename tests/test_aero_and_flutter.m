function state = test_aero_and_flutter(state)
%TEST_AERO_AND_FLUTTER  Theodorsen function, steady limits, divergence, flutter.

fprintf('\n-- unsteady aerodynamics and stability --\n');
cfg = aglas_config('baseline');

% ------------------------------------------------- Theodorsen's C(k) --
state = aglas_assert(state, 'C(0) = 1 (quasi-steady limit)', ...
    abs(theodorsen_Ck(0) - 1) < 1e-12, '');

Cinf = theodorsen_Ck(1e7);
state = aglas_assert(state, 'C(k) tends to 1/2 as k grows', ...
    abs(real(Cinf) - 0.5) < 1e-5 && abs(imag(Cinf)) < 1e-5, ...
    sprintf('%.6f%+.6fi', real(Cinf), imag(Cinf)));

% Published values (Theodorsen 1935; Bisplinghoff, Ashley and Halfman).
ref_k = [0.1, 0.2, 0.5, 1.0];
ref_C = [0.8320 - 0.1723i, 0.7276 - 0.1886i, 0.5979 - 0.1507i, 0.5394 - 0.1003i];
for i = 1:numel(ref_k)
    Ck = theodorsen_Ck(ref_k(i));
    err = abs(Ck - ref_C(i));
    state = aglas_assert(state, sprintf('C(%.1f) matches published tables', ref_k(i)), ...
        err < 5e-4, sprintf('|error| = %.2e', err));
end

state = aglas_assert(state, 'C(k) has negative imaginary part (lift lags motion)', ...
    all(imag(theodorsen_Ck([0.05 0.1 0.3 0.8 2])) < 0), '');

% --------------------------------------- steady limit of the aero matrix --
% At k = 0 the generalised lift must reduce to the steady strip-theory result.
fem = beam_fem(cfg);
md  = modal_analysis(fem, 6);
G   = modal_span_integrals(fem, md.Phi);
U   = 50;
Q0  = aero_generalised(cfg, G, 0, U);
state = aglas_assert(state, 'Q(0) is real (no phase lag in steady flow)', ...
    max(abs(imag(Q0(:)))) < 1e-12*max(abs(real(Q0(:)))), '');

% Isolate each steady coefficient with a synthetic span-integral set, so that
% one term is exercised at a time. Using a real mode for this would not work:
% a torsion mode still carries bending content, so its diagonal entry mixes
% the lift-from-twist and moment-from-twist contributions.
rho = cfg.flight.rho_air;
b   = cfg.aero.semi_chord;
a   = cfg.aero.a_ea;

Gt = struct('ww', 0, 'wp', 1, 'pw', 0, 'pp', 0);
lift_from_twist = aero_generalised(cfg, Gt, 0, U);
% Steady strip theory: dL/dalpha per unit span = q_inf*c*cl_alpha
expect_L = 0.5*rho*U^2*cfg.geom.chord*cfg.aero.cl_alpha_2d;
state = aglas_assert(state, 'steady lift per unit twist equals q_inf*c*cl_alpha', ...
    abs(lift_from_twist - expect_L) < 1e-9*expect_L, ...
    sprintf('%.6e vs %.6e N/m/rad', lift_from_twist, expect_L));

Gm = struct('ww', 0, 'wp', 0, 'pw', 0, 'pp', 1);
mom_from_twist = aero_generalised(cfg, Gm, 0, U);
% That lift acts at the quarter chord, so the moment about the elastic axis
% is the lift times the distance from the aerodynamic centre to that axis.
expect_M = expect_L * (cfg.geom.ea_frac - 0.25)*cfg.geom.chord;
state = aglas_assert(state, 'steady moment equals lift times the AC-to-EA arm', ...
    abs(mom_from_twist - expect_M) < 1e-9*abs(expect_M), ...
    sprintf('%.6e vs %.6e N m/m/rad', mom_from_twist, expect_M));

state = aglas_assert(state, 'steady plunge produces no lift', ...
    abs(aero_generalised(cfg, struct('ww',1,'wp',0,'pw',0,'pp',0), 0, U)) < 1e-9, ...
    'a wing at constant height in steady flow feels no incremental lift');

% ------------------------------------------------------- divergence --
res = flutter_pk(cfg, fem, md);

% Classical uniform cantilever wing: q_div = GJ*(pi/(2L))^2/(c*e*cl_alpha).
e_arm = (cfg.geom.ea_frac - 0.25)*cfg.geom.chord;
q_div = cfg.section.GJ*(pi/(2*cfg.geom.L))^2 / ...
        (cfg.geom.chord*e_arm*cfg.aero.cl_alpha_2d);
U_div_exact = sqrt(2*q_div/cfg.flight.rho_air);
err_div = abs(res.divergence_speed - U_div_exact)/U_div_exact;
state = aglas_assert(state, 'divergence speed matches the classical closed form', ...
    err_div < 2e-3, sprintf('%.3f vs %.3f m/s, error %.4f %%', ...
    res.divergence_speed, U_div_exact, 100*err_div));

% ---------------------------------------------------------- flutter --
state = aglas_assert(state, 'a flutter speed is found in the sweep', ...
    isfinite(res.flutter_speed), sprintf('%.2f m/s', res.flutter_speed));
state = aglas_assert(state, 'flutter frequency is positive and finite', ...
    isfinite(res.flutter_freq_hz) && res.flutter_freq_hz > 0, ...
    sprintf('%.3f Hz', res.flutter_freq_hz));
state = aglas_assert(state, 'flutter occurs below divergence', ...
    res.flutter_speed < res.divergence_speed, ...
    sprintf('%.1f < %.1f m/s', res.flutter_speed, res.divergence_speed));
state = aglas_assert(state, 'every mode is damped at the lowest sweep speed', ...
    all(res.damping(1,:) < 0), sprintf('max %.4f', max(res.damping(1,:))));
state = aglas_assert(state, 'damping becomes positive above the flutter speed', ...
    any(res.damping(res.U > res.flutter_speed, :) > 0, 'all') || ...
    any(any(res.damping(res.U > res.flutter_speed, :) > 0)), '');
state = aglas_assert(state, 'p-k iteration converges everywhere', ...
    all(res.converged(:)), sprintf('%.1f %% converged', 100*mean(res.converged(:))));
state = aglas_assert(state, 'flutter frequency lies between two structural modes', ...
    res.flutter_freq_hz > md.freq_hz(1) && res.flutter_freq_hz < md.freq_hz(end), ...
    sprintf('%.2f Hz within %.2f to %.2f Hz', res.flutter_freq_hz, ...
            md.freq_hz(1), md.freq_hz(end)));

% A bending-only model cannot flutter; the analysis must refuse rather than
% silently return NaN, which is what the original script did.
cnt = cfg; cnt.fem.include_torsion = false;
fnt = beam_fem(cnt);
mnt = modal_analysis(fnt, 4);
rnt = flutter_pk(cnt, fnt, mnt);
state = aglas_assert(state, 'bending-only model finds no flutter, as theory requires', ...
    isnan(rnt.flutter_speed), sprintf('%s', mat2str(rnt.flutter_speed)));
end
