function res = flutter_pk(cfg, fem, modes)
%FLUTTER_PK  Aeroelastic stability by the p-k method with mode tracking.
%
%   res = FLUTTER_PK(cfg, fem, modes)
%
%   Sweeps airspeed and, at each speed, solves the aeroelastic eigenproblem
%
%       [ p^2*I + p*(C_s - (b/(k*U))*Q_I) + (diag(omega^2) - Q_R) ] q = 0
%
%   where Q(k) = Q_R + i*Q_I is the generalised aerodynamic matrix from
%   Theodorsen strip theory and the mass matrix is the identity because the
%   modes are mass-normalised.
%
%   The split of Q is the heart of the p-k method. Q(k) is derived for harmonic
%   motion, so its imaginary part represents a force in quadrature with
%   displacement, i.e. proportional to velocity. Writing i*q as qdot*b/(k*U)
%   turns Im(Q) into an aerodynamic damping matrix and Re(Q) into an
%   aerodynamic stiffness matrix, which is what makes the system solvable for
%   a complex p. Because Q depends on k = Im(p)*b/U, each root is iterated to
%   self-consistency.
%
%   Flutter is the speed at which the real part of any root first crosses zero
%   from negative (damped) to positive (growing).
%
%   Divergence, a separate static instability, is found from the k = 0 limit:
%   the speed at which the aeroelastic stiffness diag(omega^2) - Q_R(0) first
%   becomes singular. That reduces to a generalised eigenproblem in U^2.
%
%   Mode tracking
%     eig returns roots in no guaranteed order, and the order changes as the
%     speed sweep progresses, so raw output produces damping curves that jump
%     between physical modes. Each root here is matched to the previous speed's
%     root by eigenvector correlation, so a branch follows one physical mode.
%
%   Output fields
%     U             airspeed sweep                                  [m/s]
%     p             tracked complex roots         [n_U x n_modes]
%     freq_hz       Im(p)/(2*pi)                  [n_U x n_modes]
%     damping       Re(p)                         [n_U x n_modes]
%     zeta          -Re(p)/|p|, positive = stable [n_U x n_modes]
%     g             2*Re(p)/|Im(p)|, the classical structural damping measure
%     k             converged reduced frequency   [n_U x n_modes]
%     flutter_speed, flutter_freq_hz, flutter_mode  (NaN if none in range)
%     divergence_speed                              (NaN if none)
%     converged     per-point iteration success flag

    U_sweep = cfg.flutter.U(:).';
    nU  = numel(U_sweep);
    nm  = modes.n_modes;
    b   = cfg.aero.semi_chord;
    zeta_s = cfg.damping.zeta;

    G  = modal_span_integrals(fem, modes.Phi);
    Ks = diag(modes.omega.^2);
    Cs = diag(2 * zeta_s * modes.omega);
    Im_ = eye(nm);
    Z  = zeros(nm);

    p_all    = nan(nU, 2*nm) + 1i*nan(nU, 2*nm);
    k_all    = nan(nU, 2*nm);
    conv_all = false(nU, 2*nm);
    vec_prev = [];

    % Initial reduced-frequency guesses from the in-vacuo frequencies.
    k_guess = [modes.omega.' * b / U_sweep(1), modes.omega.' * b / U_sweep(1)];

    for iu = 1:nU
        U = U_sweep(iu);
        p_here   = nan(1, 2*nm) + 1i*nan(1, 2*nm);
        k_here   = nan(1, 2*nm);
        vec_here = nan(2*nm, 2*nm) + 1i*nan(2*nm, 2*nm);
        conv_here = false(1, 2*nm);

        for m = 1:2*nm
            k = max(k_guess(m), 1e-6);
            ok = false;
            for it = 1:cfg.flutter.max_iter
                [p_sorted, v_sorted] = pk_roots(k, U);
                % Pick the root closest in reduced frequency to the current k.
                k_cand = abs(imag(p_sorted)) * b / U;
                [~, pick] = min(abs(k_cand - k));
                p_new = p_sorted(pick);
                k_new = max(abs(imag(p_new)) * b / U, 1e-6);
                if abs(k_new - k) < cfg.flutter.k_tol * max(1, k)
                    k = k_new;
                    p_here(m)     = p_new;
                    vec_here(:,m) = v_sorted(:,pick);
                    k_here(m)     = k;
                    ok = true;
                    break;
                end
                % Damped update keeps the iteration stable near mode crossings.
                k = k + 0.6*(k_new - k);
            end
            if ~ok
                [p_sorted, v_sorted] = pk_roots(k, U);
                k_cand = abs(imag(p_sorted)) * b / U;
                [~, pick] = min(abs(k_cand - k));
                p_here(m)     = p_sorted(pick);
                vec_here(:,m) = v_sorted(:,pick);
                k_here(m)     = k;
            end
            conv_here(m) = ok;
        end

        % --- track branches against the previous speed by eigenvector match
        if isempty(vec_prev)
            order = 1:2*nm;
        else
            order = match_modes(vec_prev, vec_here);
        end

        p_all(iu,:)    = p_here(order);
        k_all(iu,:)    = k_here(order);
        conv_all(iu,:) = conv_here(order);
        vec_prev       = vec_here(:, order);
        k_guess        = k_all(iu,:);
        k_guess(~isfinite(k_guess) | k_guess <= 0) = 1e-3;
    end

    % Keep only the roots with non-negative imaginary part; complex roots come
    % in conjugate pairs and the pair carries no extra information.
    keep = [];
    for j = 1:2*nm
        if mean(imag(p_all(:,j)), 'omitnan') >= 0
            keep(end+1) = j; %#ok<AGROW>
        end
    end
    if numel(keep) > nm
        keep = keep(1:nm);
    end

    res.U        = U_sweep(:);
    res.p        = p_all(:, keep);
    res.k        = k_all(:, keep);
    res.converged= conv_all(:, keep);
    res.freq_hz  = imag(res.p) / (2*pi);
    res.damping  = real(res.p);
    res.zeta     = -real(res.p) ./ max(abs(res.p), eps);
    res.g        = 2*real(res.p) ./ max(abs(imag(res.p)), eps);

    % ------------------------------------------------- flutter crossing
    res.flutter_speed   = NaN;
    res.flutter_freq_hz = NaN;
    res.flutter_mode    = NaN;
    for j = 1:size(res.damping, 2)
        d = res.damping(:, j);
        for i = 1:nU-1
            if isfinite(d(i)) && isfinite(d(i+1)) && d(i) < 0 && d(i+1) >= 0
                w = -d(i) / (d(i+1) - d(i));
                Uf = U_sweep(i) + w*(U_sweep(i+1) - U_sweep(i));
                if isnan(res.flutter_speed) || Uf < res.flutter_speed
                    res.flutter_speed   = Uf;
                    res.flutter_freq_hz = res.freq_hz(i,j) + ...
                        w*(res.freq_hz(i+1,j) - res.freq_hz(i,j));
                    res.flutter_mode    = j;
                end
                break;
            end
        end
    end

    % ------------------------------------------------- static divergence
    % diag(omega^2) q = U^2 * Qd q, with Qd the k = 0 aerodynamic stiffness
    % divided by U^2. Smallest positive eigenvalue gives U_div^2.
    Qd = real(aero_generalised(cfg, G, 0, 1));
    res.divergence_speed = NaN;
    if any(abs(Qd(:)) > 0)
        ev = eig(Ks, Qd);
        ev = ev(isfinite(ev) & imag(ev) == 0 & real(ev) > 0);
        if ~isempty(ev)
            res.divergence_speed = sqrt(min(real(ev)));
        end
    end

    res.n_modes = nm;
    res.modes   = modes;

    % =================================================== nested helpers
    function [ps, vs] = pk_roots(kq, Uq)
        Q  = aero_generalised(cfg, G, kq, Uq);
        QR = real(Q);
        QI = imag(Q);
        C_ae = Cs - (b/(kq*Uq)) * QI;
        K_ae = Ks - QR;
        A = [ Z,      Im_;
             -K_ae,  -C_ae ];
        [vv, dd] = eig(A);
        ps = diag(dd).';
        vs = vv;
    end
end

% ------------------------------------------------------------------------
function order = match_modes(V_prev, V_now)
%MATCH_MODES  Greedy assignment of current eigenvectors to previous branches.
    n = size(V_prev, 2);
    Cmat = zeros(n);
    for a = 1:n
        for bb = 1:n
            va = V_prev(:,a); vb = V_now(:,bb);
            na = norm(va); nb = norm(vb);
            if na == 0 || nb == 0
                Cmat(a,bb) = 0;
            else
                Cmat(a,bb) = abs(va' * vb) / (na*nb);   % modal assurance
            end
        end
    end
    order = zeros(1, n);
    used  = false(1, n);
    for step = 1:n
        [~, idx] = max(Cmat(:));
        [a, bb]  = ind2sub(size(Cmat), idx);
        order(a) = bb;
        used(bb) = true;
        Cmat(a,:)  = -Inf;
        Cmat(:,bb) = -Inf;
    end
    miss = find(order == 0);
    free = find(~used);
    order(miss) = free(1:numel(miss));
end
