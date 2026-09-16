function state = test_modal(state)
%TEST_MODAL  Eigensolution against the closed-form cantilever.

fprintf('\n-- modal analysis --\n');
cfg = aglas_config('baseline');
cfg.fem.include_torsion = false;      % compare like with like
cfg.fem.n_elements = 20;

fem = beam_fem(cfg);
md  = modal_analysis(fem, 6);

beta = [1.875104068, 4.694091133, 7.854757438, 10.99554073];
f_ex = beta.^2/(2*pi) * sqrt(cfg.section.EI/(cfg.section.m_total*cfg.geom.L^4));

for i = 1:4
    e = 100*abs(md.freq_hz(i) - f_ex(i))/f_ex(i);
    state = aglas_assert(state, sprintf('bending mode %d within 0.01 %% of exact', i), ...
        e < 0.01, sprintf('error %.2e %%', e));
end

state = aglas_assert(state, 'frequencies are sorted ascending', ...
    all(diff(md.freq_hz) >= 0), '');

G = md.Phi.'*fem.M*md.Phi;
state = aglas_assert(state, 'mode shapes are mass-orthonormal', ...
    norm(G - eye(size(G)), 'fro') < 1e-9, ...
    sprintf('||Phi''*M*Phi - I|| = %.2e', norm(G - eye(size(G)), 'fro')));

Kd = md.Phi.'*fem.K*md.Phi;
state = aglas_assert(state, 'modal stiffness equals diag(omega^2)', ...
    norm(Kd - diag(md.omega.^2), 'fro') < 1e-6*max(md.omega.^2), '');

state = aglas_assert(state, 'clamped root has zero displacement in every mode', ...
    all(abs(md.w_shape(1,:)) < 1e-14), '');

% With torsion on, both families must be present and correctly labelled.
ct = aglas_config('baseline');
mdt = modal_analysis(beam_fem(ct), 8);
state = aglas_assert(state, 'both bending and torsion modes are identified', ...
    any(strcmp(mdt.type,'bending')) && any(strcmp(mdt.type,'torsion')), ...
    sprintf('%d bending, %d torsion', sum(strcmp(mdt.type,'bending')), ...
            sum(strcmp(mdt.type,'torsion'))));
end
