function ff = feedforward_design(cfg, P, base, opts)
%FEEDFORWARD_DESIGN  Gust feedforward from a probe ahead of the wing.
%
%   ff = FEEDFORWARD_DESIGN(cfg, P, base)
%   ff = FEEDFORWARD_DESIGN(cfg, P, base, opts)
%
%   Builds the feedforward path that CLOSED_LOOP_SIM adds to the LQG command
%   when ctrl.ff is set. A probe probe_ahead metres ahead of the wing measures
%   the gust, so the wing sees it tau = probe_ahead/U later. An FIR filter of
%   n_taps samples at Ts maps the recent reference history to a surface
%   command:
%
%       u_ff(j) = w' * [r(j); r(j-1); ...; r(j-n_taps+1)]
%
%   Two ways to set w, which is the comparison this exists for:
%
%     fixed     the Wiener solution at the design condition, from
%               FEEDFORWARD_DESIGN(..., 'train', t, wg): the weights that
%               minimise root strain in a training record, held thereafter
%     adaptive  normalised filtered-x LMS (ff.adapt = true), the algorithm of
%               Wildschek et al. (AIAA 2006-6054) and Zeng et al. (JGCD 2010):
%
%                   w <- (1-leak)*w - mu * e * rf / (eps + rf'*rf)
%
%               with e the measured root strain and rf the reference filtered
%               through S, a model of the path from u_ff to root strain.
%
%   S is the impulse response of the NOMINAL closed loop (plant plus LQG
%   estimator and gain, actuator limits ignored) sampled at Ts. It is
%   deliberately not updated off-design: filtered-x LMS needs only an
%   approximate secondary path, and that robustness is part of what is being
%   tested.
%
%   opts fields (defaults)
%     probe_ahead  probe distance ahead of the wing        10 m
%     Ts           feedforward sample period               0.005 s
%     n_taps       FIR length                              100 (0.5 s)
%     S_len        secondary-path model length             1.0 s
%     mu           LMS step                                0.01
%     leak         LMS leakage per update                  0
%     adapt        adapt online                            false
%     w0           initial weights                         zeros
%     train        {t, wg}: solve the Wiener weights from this record

    if nargin < 4, opts = struct(); end
    g = @(f, d) getf(opts, f, d);

    ff.probe_ahead = g('probe_ahead', 10);
    ff.Ts          = g('Ts', 0.005);
    n_taps         = g('n_taps', 100);
    S_len          = g('S_len', 1.0);
    ff.mu          = g('mu', 0.01);
    ff.leak        = g('leak', 0);
    ff.adapt       = g('adapt', false);
    ff.err_row     = 1;                    % measurement 1, root strain
    ff.eps         = 1e-14;
    ff.S           = secondary_path(P, base, ff.Ts, round(S_len/ff.Ts), ff.err_row);
    ff.w           = g('w0', zeros(n_taps, 1));
    ff.w           = ff.w(:);

    if isfield(opts, 'train') && ~isempty(opts.train)
        ff.w = wiener(cfg, P, base, ff, opts.train{1}, opts.train{2});
    end
end

% ------------------------------------------------------------------------
function S = secondary_path(P, base, Ts, n, row)
%SECONDARY_PATH  Sampled impulse response from feedforward command to strain.
%   Closed loop with u = -K*xhat + v; the estimator sees the total command.
    K  = base.K;  L = base.L;
    G  = base.Be - L*base.De + L*P.D;
    Acl = [P.A,      -P.B*K;
           L*P.C,     base.Ae - L*base.Ce - G*K];
    Bcl = [P.B; G];
    Ccl = [P.C(row,:), -P.D(row)*K];
    Dcl = P.D(row);

    nx = size(Acl, 1);
    E  = expm([Acl, Bcl; zeros(1, nx + 1)] * Ts);
    Ad = E(1:nx, 1:nx);
    Bd = E(1:nx, nx+1);

    S = zeros(n, 1);
    S(1) = Dcl;
    v = Bd;
    for i = 2:n
        S(i) = Ccl * v;
        v = Ad * v;
    end
end

% ------------------------------------------------------------------------
function w = wiener(cfg, P, base, ff, t, wg)
%WIENER  Least-squares feedforward weights from an LQG-only training run.
%   With feedforward the sampled strain is, linearly, e = e0 + X*w, where e0
%   is the strain under feedback alone and X stacks the delayed filtered
%   references. Minimising |e|^2 gives the normal equations below.
    r0   = closed_loop_sim(cfg, P, base, t, wg);
    dt   = t(2) - t(1);
    nsub = round(ff.Ts/dt);
    ks   = 1:nsub:numel(t);
    tau  = ff.probe_ahead / cfg.flight.U_inf;
    r    = interp1(t(:), wg(:), t(ks) + tau, 'linear', wg(end));
    rf   = filter(ff.S, 1, r(:));
    e0   = r0.strain(ks).';

    N = numel(ff.w);
    X = zeros(numel(ks), N);
    for m = 1:N
        X(m:end, m) = rf(1:end-m+1);
    end
    H   = X.'*X;
    reg = 1e-8 * trace(H) / N;
    w   = -(H + reg*eye(N)) \ (X.'*e0);
end

function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
