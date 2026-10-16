function state = test_unsteady_plant(state)
%TEST_UNSTEADY_PLANT  Indicial-lag control plant against Theodorsen and p-k.

fprintf('\n-- unsteady control plant --\n');

cfg = aglas_config('baseline');
fem = beam_fem(cfg);
md  = modal_analysis(fem, cfg.fem.n_modes);
cs  = control_surface(cfg, fem, md);

cq = cfg; cq.aero.plant_model = 'quasisteady';
cu = cfg; cu.aero.plant_model = 'unsteady';
cu.aero.indicial = 'jones'; cu.aero.compressible = false;
Pq = plant_statespace(cq, fem, md, cs);
Pu = plant_statespace(cu, fem, md, cs);

% Every lag is 1 at zero frequency and apparent mass carries no static
% force, so the static response to a steady gust must be identical.
dc_q = Pq.Cz(1,:) * (-Pq.A \ Pq.Bw);
dc_u = Pu.Cz(1,:) * (-Pu.A \ Pu.Bw);
state = aglas_assert(state, 'unsteady plant has the quasi-steady static gust response', ...
    abs(dc_u - dc_q) < 1e-8*abs(dc_q), sprintf('%.6e vs %.6e N m per m/s', dc_u, dc_q));

% The quasi-steady stiffness must not depend on a small-k stand-in: it is
% the C = 1 strip-theory value, which the k -> 0 limit of the full
% generalised aerodynamic matrix approaches.
G  = modal_span_integrals(fem, md.Phi);
Kk = real(aero_generalised(cfg, G, 1e-7, cfg.flight.U_inf));
state = aglas_assert(state, 'quasi-steady aerodynamic stiffness is the C = 1 limit', ...
    norm(Pq.K_aero - Kk) < 1e-4*norm(Kk), sprintf('relative %.2e', norm(Pq.K_aero - Kk)/norm(Kk)));

% Jones' two-pole fit reproduces Theodorsen's C(k) to about 0.015.
A = Pu.indicial.A; b = Pu.indicial.b;
k = [0.05 0.1 0.2 0.5 1.0];
err = zeros(size(k));
for i = 1:numel(k)
    Cj = (1 - sum(A)) + sum(A .* b ./ (1i*k(i) + b));
    err(i) = abs(Cj - theodorsen_Ck(k(i)));
end
state = aglas_assert(state, 'indicial lag reproduces Theodorsen C(k)', ...
    max(err) < 0.02, sprintf('max |error| %.4f over k = 0.05 to 1', max(err)));

% Stability boundary of the incompressible unsteady plant against p-k.
md6 = modal_analysis(fem, cfg.flutter.n_modes);
pk  = flutter_pk(cfg, fem, md6);
[U_ss, f_ss] = boundary(cu, fem, md6, 280:1:340);
state = aglas_assert(state, 'unsteady plant goes unstable at the p-k flutter speed', ...
    abs(U_ss - pk.flutter_speed) < 0.02*pk.flutter_speed, ...
    sprintf('%.0f vs %.1f m/s', U_ss, pk.flutter_speed));
state = aglas_assert(state, '... and at the p-k flutter frequency', ...
    abs(f_ss - pk.flutter_freq_hz) < 0.02*pk.flutter_freq_hz, ...
    sprintf('%.2f vs %.2f Hz', f_ss, pk.flutter_freq_hz));

% Without the circulation lag the boundary is badly wrong, which is why the
% quasi-steady plant must not be used near it.
U_qs = boundary(cq, fem, md6, 150:2:420);
state = aglas_assert(state, 'quasi-steady plant misses the flutter speed by over 20 %', ...
    abs(U_qs - pk.flutter_speed) > 0.2*pk.flutter_speed, ...
    sprintf('%.0f vs %.1f m/s', U_qs, pk.flutter_speed));

% Compressibility raises the lift slope and lowers the boundary.
cc = cu; cc.aero.compressible = true;
U_c = boundary(cc, fem, md6, 220:1:300);
state = aglas_assert(state, 'compressibility lowers the flutter speed', ...
    U_c < U_ss, sprintf('%.0f vs %.0f m/s incompressible', U_c, U_ss));
end

function [Uc, fc] = boundary(cfg, fem, md, U_range)
    Uc = NaN; fc = NaN;
    for U = U_range
        c = cfg; c.flight.U_inf = U;
        P = plant_statespace(c, fem, md, control_surface(c, fem, md));
        ev = eig(P.A);
        [mr, i] = max(real(ev));
        if mr > 0
            Uc = U; fc = abs(imag(ev(i))) / (2*pi);
            return;
        end
    end
end
