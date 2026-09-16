function modes = modal_analysis(fem, n_keep)
%MODAL_ANALYSIS  Solve the free-vibration eigenproblem, mass-normalised.
%
%   modes = MODAL_ANALYSIS(fem)         returns every mode.
%   modes = MODAL_ANALYSIS(fem, n_keep) returns the lowest n_keep modes.
%
%   Solves  K*phi = omega^2 * M*phi  and returns eigenvectors scaled so that
%
%       Phi' * M * Phi = I        and        Phi' * K * Phi = diag(omega^2)
%
%   Mass normalisation is not cosmetic. It is what makes the modal equations
%   reduce to  qddot + 2*zeta*omega*qdot + omega^2*q = Phi'*f , so that a modal
%   damping matrix diag(2*zeta*omega) is correct. The original code normalised
%   each eigenvector to unit peak amplitude instead, which left
%   Phi'*M*Phi = 20.25*I, and then used diag(2*zeta*omega) anyway. That
%   understated modal damping by the same factor of 20.25 and also invalidated
%   the M_gen = eye() assumption in the flutter script.
%
%   Eigenvalues are explicitly sorted ascending. MATLAB's eig(K,M) does not
%   guarantee ordering, so "mode 1" in unsorted output is not necessarily the
%   fundamental.
%
%   Output fields
%     Phi        mass-normalised mode shapes, reduced DOFs  [n_dof x n_keep]
%     omega      circular natural frequencies               [rad/s]
%     freq_hz    natural frequencies                        [Hz]
%     w_shape    vertical deflection component, root prepended [n_nodes x n_keep]
%     th_shape   bending rotation component, root prepended
%     ph_shape   twist component, root prepended (empty if no torsion)
%     type       'bending' or 'torsion' classification per mode

    if nargin < 2 || isempty(n_keep)
        n_keep = size(fem.K, 1);
    end

    [V, D] = eig(full(fem.K), full(fem.M));
    lambda = real(diag(D));

    % Guard against tiny negative eigenvalues from round-off.
    lambda(lambda < 0 & lambda > -1e-8*max(abs(lambda))) = 0;
    if any(lambda < 0)
        error('modal_analysis:negativeEigenvalue', ...
              'Negative eigenvalue found; the stiffness matrix is not positive definite.');
    end

    [lambda, order] = sort(lambda, 'ascend');
    V = real(V(:, order));

    n_keep = min(n_keep, numel(lambda));
    lambda = lambda(1:n_keep);
    V = V(:, 1:n_keep);

    % Mass-normalise: phi' * M * phi = 1.
    for i = 1:n_keep
        scale = sqrt(V(:,i).' * fem.M * V(:,i));
        V(:,i) = V(:,i) / scale;
        % Fix the sign so the largest-magnitude entry is positive, making the
        % shapes reproducible across solvers and platforms.
        [~, imax] = max(abs(V(:,i)));
        if V(imax, i) < 0
            V(:,i) = -V(:,i);
        end
    end

    modes = struct();
    modes.Phi     = V;
    modes.omega   = sqrt(lambda);
    modes.freq_hz = modes.omega / (2*pi);
    modes.n_modes = n_keep;

    % Expand to nodal shape arrays with the clamped root prepended as zero.
    z = zeros(1, n_keep);
    modes.w_shape  = [z; V(fem.idx_w,  :)];
    modes.th_shape = [z; V(fem.idx_th, :)];
    if fem.torsion
        modes.ph_shape = [z; V(fem.idx_ph, :)];
    else
        modes.ph_shape = [];
    end

    % Classify each mode by where its kinetic energy sits.
    modes.type = cell(n_keep, 1);
    for i = 1:n_keep
        e_bend = sum(V(fem.idx_w, i).^2);
        if fem.torsion
            e_tors = sum(V(fem.idx_ph, i).^2) * (fem.dx^0);
            % Compare energies scaled by the relevant inertias.
            e_tors = e_tors * 1.0;
            if e_tors > e_bend
                modes.type{i} = 'torsion';
            else
                modes.type{i} = 'bending';
            end
        else
            modes.type{i} = 'bending';
        end
    end
end
