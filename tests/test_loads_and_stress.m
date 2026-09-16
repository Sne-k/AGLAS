function state = test_loads_and_stress(state)
%TEST_LOADS_AND_STRESS  Consistent loads and stress recovery against statics.

fprintf('\n-- consistent loads and stress recovery --\n');
cfg = aglas_config('baseline');
cfg.fem.include_torsion = false;
fem = beam_fem(cfg);

q  = 1000;                                  % uniform line load [N/m]
L  = cfg.geom.L;
EI = cfg.section.EI;

f = consistent_line_load(fem, q, 0);
u = fem.K \ f;

% Static tip deflection of a uniformly loaded cantilever, q*L^4/(8*EI).
w_tip   = u(fem.idx_w(end));
w_exact = q*L^4/(8*EI);
e_w = abs(w_tip - w_exact)/w_exact;
state = aglas_assert(state, 'tip deflection equals q*L^4/(8*EI)', ...
    e_w < 1e-9, sprintf('relative error %.2e', e_w));

% Tip slope, q*L^3/(6*EI).
th_tip   = u(fem.idx_th(end));
th_exact = q*L^3/(6*EI);
e_t = abs(th_tip - th_exact)/th_exact;
state = aglas_assert(state, 'tip slope equals q*L^3/(6*EI)', ...
    e_t < 1e-9, sprintf('relative error %.2e', e_t));

% Recovered root moment, q*L^2/2, and its convergence order.
rec = recover_bending(cfg, fem, u);
M_exact = q*L^2/2;
e_M = abs(rec.moment(1) - M_exact)/M_exact;
state = aglas_assert(state, 'root moment within 0.1 % of q*L^2/2', ...
    e_M < 1e-3, sprintf('error %.4f %%', 100*e_M));

state = aglas_assert(state, 'root carries the largest bending moment', ...
    abs(rec.moment(1)) == max(abs(rec.moment)), '');

state = aglas_assert(state, 'tip moment is small compared with the root', ...
    abs(rec.moment(end)) < 1e-3*abs(rec.moment(1)), ...
    sprintf('%.3e vs %.3e N m', abs(rec.moment(end)), abs(rec.moment(1))));

% Second-order convergence of the recovered moment.
errs = zeros(1,3); nels = [10 20 40];
for i = 1:3
    c = cfg; c.fem.n_elements = nels(i);
    fm = beam_fem(c);
    ui = fm.K \ consistent_line_load(fm, q, 0);
    ri = recover_bending(c, fm, ui);
    errs(i) = abs(ri.moment(1) - M_exact)/M_exact;
end
order = log(errs(1)/errs(3))/log(nels(3)/nels(1));
state = aglas_assert(state, 'moment recovery converges at second order', ...
    abs(order - 2) < 0.15, sprintf('observed order %.3f', order));

% Stress is consistent with sigma = M*c/I at every station.
sig_from_M = rec.moment * cfg.section.c_outer / cfg.section.I;
state = aglas_assert(state, 'stress is consistent with M*c/I', ...
    max(abs(rec.stress - sig_from_M)) < 1e-6*max(abs(rec.stress)), '');

% Total applied force. The root node's share is reacted by the clamp and so
% does not appear among the free DOFs, hence the expected q*L - q*dx/2.
F_free = sum(f(fem.idx_w));
F_expect = q*L - q*fem.dx/2;
state = aglas_assert(state, 'consistent load vector totals q*L minus the root share', ...
    abs(F_free - F_expect) < 1e-9*abs(F_expect), ...
    sprintf('%.4f vs %.4f N', F_free, F_expect));

% Refining the mesh must not change the total load: the original code applied
% a unit force per node with no dx scaling, which failed exactly this test.
c2 = cfg; c2.fem.n_elements = 40;
fm2 = beam_fem(c2);
f2 = consistent_line_load(fm2, q, 0);
tot1 = sum(f(fem.idx_w))  + q*fem.dx/2;
tot2 = sum(f2(fm2.idx_w)) + q*fm2.dx/2;
state = aglas_assert(state, 'total load is mesh independent', ...
    abs(tot1 - tot2) < 1e-9*tot1, sprintf('%.4f vs %.4f N', tot1, tot2));

% Torsional load path.
ct = aglas_config('baseline');
ft = beam_fem(ct);
mt = 500;                                   % uniform torque [N m/m]
ftq = consistent_line_load(ft, 0, mt);
ut  = ft.K \ ftq;
phi_tip   = ut(ft.idx_ph(end));
phi_exact = mt*ct.geom.L^2/(2*ct.section.GJ);
e_p = abs(phi_tip - phi_exact)/phi_exact;
state = aglas_assert(state, 'tip twist equals m*L^2/(2*GJ)', ...
    e_p < 1e-9, sprintf('relative error %.2e', e_p));
end
