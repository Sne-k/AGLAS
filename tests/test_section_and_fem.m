function state = test_section_and_fem(state)
%TEST_SECTION_AND_FEM  Section properties and finite element assembly.

fprintf('\n-- section properties and FE assembly --\n');
cfg = aglas_config('baseline');
g = cfg.geom; s = cfg.section;

% Box area and second moment against direct hand formulas.
wi = g.box_w - 2*g.wall_t;  hi = g.box_h - 2*g.wall_t;
A_hand = g.box_w*g.box_h - wi*hi;
I_hand = (g.box_w*g.box_h^3 - wi*hi^3)/12;
state = aglas_assert(state, 'box area matches hand calculation', ...
    abs(s.A - A_hand) < 1e-12, sprintf('%.3e vs %.3e', s.A, A_hand));
state = aglas_assert(state, 'box I matches hand calculation', ...
    abs(s.I - I_hand) < 1e-15, sprintf('%.6e vs %.6e', s.I, I_hand));

% Bredt-Batho torsion constant.
Am = (g.box_w - g.wall_t)*(g.box_h - g.wall_t);
sm = 2*((g.box_w - g.wall_t) + (g.box_h - g.wall_t));
J_hand = 4*Am^2*g.wall_t/sm;
state = aglas_assert(state, 'box J matches Bredt-Batho formula', ...
    abs(s.J - J_hand) < 1e-15, sprintf('%.6e vs %.6e', s.J, J_hand));

% Solid section: closed forms.
cl = aglas_config('legacy'); gl = cl.geom; sl = cl.section;
state = aglas_assert(state, 'solid section I = b*h^3/12', ...
    abs(sl.I - gl.box_w*gl.box_h^3/12) < 1e-18, sprintf('%.6e', sl.I));
state = aglas_assert(state, 'solid J below the thin-strip bound b*h^3/3', ...
    sl.J > 0 && sl.J < gl.box_w*gl.box_h^3/3, sprintf('%.6e', sl.J));

% FE assembly invariants.
fem = beam_fem(cfg);
state = aglas_assert(state, 'K is symmetric', ...
    norm(fem.K - fem.K.', 'fro') < 1e-9*norm(fem.K, 'fro'), '');
state = aglas_assert(state, 'M is symmetric', ...
    norm(fem.M - fem.M.', 'fro') < 1e-9*norm(fem.M, 'fro'), '');
state = aglas_assert(state, 'K is positive definite', all(eig(fem.K) > 0), '');
state = aglas_assert(state, 'M is positive definite', all(eig(fem.M) > 0), '');

m_assembled = sum(sum(fem.M_full(1:fem.dpn:end, 1:fem.dpn:end)));
m_exact = s.m_total * g.L;
state = aglas_assert(state, 'translational mass is conserved exactly', ...
    abs(m_assembled - m_exact) < 1e-9*m_exact, ...
    sprintf('%.6f vs %.6f kg', m_assembled, m_exact));

% Torsion DOF count.
state = aglas_assert(state, 'torsion adds a third DOF per node', ...
    fem.dpn == 3 && numel(fem.idx_ph) == fem.n_nodes-1, ...
    sprintf('dpn = %d', fem.dpn));

cnt = cfg; cnt.fem.include_torsion = false;
fnt = beam_fem(cnt);
state = aglas_assert(state, 'torsion can be disabled, leaving 2 DOF per node', ...
    fnt.dpn == 2 && isempty(fnt.idx_ph), sprintf('dpn = %d', fnt.dpn));
end
