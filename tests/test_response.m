function state = test_response(state)
%TEST_RESPONSE  Dynamic response and energy behaviour of the modal solver.

fprintf('\n-- dynamic response --\n');
cfg = aglas_config('baseline');
fem = beam_fem(cfg);
modes = modal_analysis(fem, cfg.fem.n_modes);

t_grid = 0:cfg.time.dt_out:cfg.time.t_end;
wg = gust_one_minus_cos(t_grid, cfg.gust.U_ds, cfg.gust.t_g);
[t_sol, q_sol, U_dof, gi] = solve_modal_response(cfg, fem, modes, t_grid, wg, ...
                                'kussner', cfg.gust.t_g/20);

tip = U_dof(:, fem.idx_w(end));

state = aglas_assert(state, 'response starts from rest', ...
    abs(tip(1)) < 1e-12, sprintf('%.2e m', tip(1)));
state = aglas_assert(state, 'response is bounded (no divergence)', ...
    all(isfinite(tip)) && max(abs(tip)) < cfg.geom.L, ...
    sprintf('peak %.4f m', max(abs(tip))));

% The gust ends at t_g; afterwards the wing must ring down, not grow.
post = t_sol > cfg.gust.t_g;
tp = tip(post); ts = t_sol(post);
first_half  = max(abs(tp(ts < (ts(1)+ts(end))/2)));
second_half = max(abs(tp(ts >= (ts(1)+ts(end))/2)));
state = aglas_assert(state, 'free response decays after the gust passes', ...
    second_half < first_half, ...
    sprintf('%.4e then %.4e m', first_half, second_half));

% Logarithmic decrement must recover the prescribed damping ratio. The decay
% is dominated by the fundamental mode once the higher modes have died away.
env = abs(tp);
pk = find(env(2:end-1) > env(1:end-2) & env(2:end-1) >= env(3:end)) + 1;
pk = pk(env(pk) > 0.05*max(env));          % ignore the noise floor
if numel(pk) >= 4
    % Fit log(peak) against time; the slope is -zeta*omega_1.
    cf = polyfit(ts(pk), log(env(pk)), 1);
    zeta_est = -cf(1) / modes.omega(1);
    state = aglas_assert(state, 'decay rate recovers the prescribed damping ratio', ...
        abs(zeta_est - cfg.damping.zeta) < 0.004, ...
        sprintf('estimated %.4f vs prescribed %.4f over %d peaks', ...
                zeta_est, cfg.damping.zeta, numel(pk)));
else
    state = aglas_assert(state, 'decay rate recovers the prescribed damping ratio', ...
        false, sprintf('only %d envelope peaks found', numel(pk)));
end

% Quasi-static check: the peak dynamic deflection must be close to, and at
% least as large as, the static deflection under the same peak load, because
% the gust is slow relative to the structure.
f_stat = consistent_line_load(fem, max(abs(gi.q_line)), max(abs(gi.m_line)));
w_stat = (fem.K \ f_stat);
w_stat = w_stat(fem.idx_w(end));
daf = max(abs(tip))/w_stat;
state = aglas_assert(state, 'dynamic amplification factor is between 1.0 and 1.3', ...
    daf > 0.99 && daf < 1.3, sprintf('%.4f', daf));

% Doubling the gust amplitude must double the response: the model is linear.
wg2 = gust_one_minus_cos(t_grid, 2*cfg.gust.U_ds, cfg.gust.t_g);
c2 = cfg; c2.gust.U_ds = 2*cfg.gust.U_ds;
[~, ~, U2] = solve_modal_response(c2, fem, modes, t_grid, wg2, ...
                                  'kussner', cfg.gust.t_g/20);
tip2 = U2(:, fem.idx_w(end));
ratio = max(abs(tip2))/max(abs(tip));
state = aglas_assert(state, 'response scales linearly with gust amplitude', ...
    abs(ratio - 2) < 1e-3, sprintf('ratio %.6f', ratio));

% Stress recovered from the response must peak at the root for a cantilever
% under a distributed upward load.
rec = recover_bending(cfg, fem, U_dof.');
[~, k] = max(max(abs(rec.stress), [], 2));
state = aglas_assert(state, 'peak stress occurs at the root', ...
    k == 1, sprintf('peak at node %d, x = %.2f m', k, rec.x(k)));

% Adding modes must not change the answer much: the solution is converged in
% the modal basis.
c3 = cfg; c3.fem.n_modes = 6;
m6 = modal_analysis(fem, 6);
[~, ~, U6] = solve_modal_response(c3, fem, m6, t_grid, wg, ...
                                  'kussner', cfg.gust.t_g/20);
tip6 = U6(:, fem.idx_w(end));
d = abs(max(abs(tip6)) - max(abs(tip)))/max(abs(tip6));
state = aglas_assert(state, 'response is converged in the number of modes', ...
    d < 0.02, sprintf('%.3f %% change from 3 to 6 modes', 100*d));
end
