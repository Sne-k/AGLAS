function ctrl = lqg_design(cfg, P, weights)
%LQG_DESIGN  Regulator and gust-estimating observer for load alleviation.
%
%   ctrl = LQG_DESIGN(cfg, P)
%   ctrl = LQG_DESIGN(cfg, P, weights)
%
%   weights overrides cfg.lqr:
%     .q_moment   price on root bending moment
%     .q_tip      price on tip deflection
%     .r_command  price on commanded surface deflection
%     .moment_ref, .tip_ref   Bryson normalisers (optional)
%
%   Why the gust is a state
%     A first attempt used a plain observer of the structural states alone.
%     That observer has no way to represent the gust, yet the gust appears
%     directly in the tip accelerometer through the feedthrough term Dw. The
%     filter could only explain that signal as structural motion, so it
%     produced a badly wrong estimate and drove the surface into its stops. The
%     symptom was unmistakable: output feedback appeared to beat full-state
%     feedback, which is impossible, because full-state LQR is the optimum for
%     this cost and nothing measuring less can do better.
%
%     The gust is therefore modelled explicitly as an extra state, a
%     first-order Markov process
%
%         wdot = -(1/T)*w + noise,     T = gust duration / 2
%
%     which is the standard shaping filter for a band-limited vertical gust.
%     The observer now estimates the gust alongside the structure and can
%     separate the two.
%
%   Why that also improves the regulator
%     With the gust in the state vector, the regulator Riccati solution
%     produces a gain component on it. That component is disturbance
%     feedforward: the controller acts on the estimated gust rather than
%     waiting for the wing to respond to it. For load alleviation this is the
%     part that does most of the work, because a purely reactive loop can never
%     get ahead of a disturbance arriving at the structure's own timescale.
%     The gust state is uncontrollable but stable, so the pair stays
%     stabilisable and the Riccati solution exists.
%
%   Bryson scaling
%     A bending moment of order 1e5 N m and a deflection of order 0.3 rad
%     cannot be compared with raw unit weights. Each term is divided by the
%     square of a representative maximum, so the weights are dimensionless and
%     of order one, and r_command becomes the single meaningful tuning knob.
%
%   Both Riccati equations go through CARE_SOLVE, so no toolbox is needed.

    if nargin < 3 || isempty(weights)
        weights = cfg.lqr;
    end

    ns = P.n_states;

    % ------------------------------------------------- gust shaping filter
    T_gust = max(cfg.gust.t_g/2, 1e-3);
    a_g    = 1/T_gust;

    Aa = [ P.A,          P.Bw;
           zeros(1,ns), -a_g ];
    Ba = [ P.B; 0 ];
    Ca = [ P.C, P.Dw ];
    Da = P.D;
    Cza = [ P.Cz, zeros(size(P.Cz,1), 1) ];

    na = ns + 1;

    % --------------------------------------------------- Bryson normalisers
    if isfield(weights,'moment_ref') && ~isempty(weights.moment_ref)
        scale.moment = weights.moment_ref;
    else
        scale.moment = 1.45e5;              % N m, order of the open-loop peak
    end
    if isfield(weights,'tip_ref') && ~isempty(weights.tip_ref)
        scale.tip = weights.tip_ref;
    else
        scale.tip = 0.7;                    % m
    end
    scale.delta = deg2rad(cfg.control.delta_max);

    Qz = diag([ weights.q_moment / scale.moment^2, ...
                weights.q_tip    / scale.tip^2 ]);
    R  = weights.r_command / scale.delta^2;

    Q = Cza.' * Qz * Cza;
    Q = (Q + Q.')/2;

    % ------------------------------------------------------- regulator
    [~, K, info_lqr] = care_solve(Aa, Ba, Q, R);

    % ------------------------------------------------------- estimator
    % Process noise: intensity on the gust state, plus a small floor on the
    % rest so the filter Riccati stays well conditioned.
    sigma_w = cfg.gust.U_ds;                 % scale of the gust being rejected
    W = 1e-8 * eye(na);
    W(na, na) = cfg.sensor.gust_psd * sigma_w^2 * 2 * a_g;
    W = (W + W.')/2;
    V = diag([cfg.sensor.strain_noise^2, cfg.sensor.accel_noise^2]);

    [~, Lt, info_kf] = care_solve(Aa.', Ca.', W, V);
    L = Lt.';

    ctrl.type     = 'lqg';
    ctrl.K        = K;            % gain on the augmented estimate
    ctrl.L        = L;
    ctrl.Ae       = Aa;           % estimator model
    ctrl.Be       = Ba;
    ctrl.Ce       = Ca;
    ctrl.De       = Da;
    ctrl.Cza      = Cza;
    ctrl.n_aug    = na;
    ctrl.idx_gust = na;
    ctrl.Q        = Q;
    ctrl.R        = R;
    ctrl.Qz       = Qz;
    ctrl.W        = W;
    ctrl.V        = V;
    ctrl.scale    = scale;
    ctrl.weights  = weights;
    ctrl.T_gust   = T_gust;
    ctrl.info_lqr = info_lqr;
    ctrl.info_kf  = info_kf;
    ctrl.A_cl     = Aa - Ba*K;
    ctrl.A_est    = Aa - L*Ca;

    % Bandwidths, so a caller can see whether a chosen integration step can
    % actually resolve the loops it is about to simulate.
    ctrl.bw.regulator_hz = max(abs(eig(ctrl.A_cl)))  / (2*pi);
    ctrl.bw.observer_hz  = max(abs(eig(ctrl.A_est))) / (2*pi);
    ctrl.bw.plant_hz     = max(abs(eig(P.A)))        / (2*pi);
    ctrl.P        = P;

    % Split the gain for reporting: feedback on the structure, feedforward on
    % the estimated gust.
    ctrl.K_fb = K(1:ns);
    ctrl.K_ff = K(na);
end
