function ctrl = mrac_design(cfg, P, base, opts)
%MRAC_DESIGN  Model reference adaptive augmentation of the LQG baseline.
%
%   ctrl = MRAC_DESIGN(cfg, P, base)
%   ctrl = MRAC_DESIGN(cfg, P, base, opts)
%
%   base is the LQG controller from LQG_DESIGN. opts may set:
%     .gamma   adaptation rate (scalar or ns-vector), default 5e3
%     .sigma   sigma-modification leakage, default 1e-2
%     .theta_max  bound on the adaptive parameter norm, default 5
%
%   Why adapt at all
%     Adaptive control is only worth its complexity when there is real
%     parametric variation to adapt to, and here there is. Control surface
%     effectiveness scales with dynamic pressure, so a surface authority
%     calibrated at 50 m/s is off by a factor of 2.6 at 80 m/s and by half at
%     35 m/s. Wing mass changes with fuel state, moving every modal frequency.
%     A fixed-gain LQG is optimal at exactly one flight condition and merely
%     adequate away from it. That is the honest justification for the word
%     "adaptive" in this project's name, rather than adopting the technique
%     because the acronym promises it.
%
%   Architecture: augmentation, not replacement
%     The adaptive element is added to a working LQG loop rather than used
%     alone:
%
%         u = -K*xhat_aug  +  theta_hat' * Phi(xhat)
%
%     The baseline does the nominal job and provides the gust estimate that
%     makes feedforward possible; the adaptive term absorbs whatever the
%     baseline's model got wrong. This is the standard aerospace pattern, and
%     it fails gracefully: if adaptation is switched off, or the parameters
%     drift to zero through the leakage term, what remains is the LQG design,
%     not an unstabilised aircraft.
%
%   Reference model and adaptive law
%     The reference model is the nominal closed loop driven by the estimated
%     gust,  xm_dot = (A - B*K_fb)*xm + Bw*w_hat,  so the error
%     e = xhat_struct - xm is the amount by which the real aircraft fails to
%     behave like its design model. With a matched uncertainty the Lyapunov
%     function V = e'*P_lyap*e + trace(theta_tilde'*inv(Gamma)*theta_tilde*Lambda)
%     gives the update
%
%         theta_hat_dot = -Gamma * ( Phi * (B'*P_lyap*e)/(1 + Phi'*Phi)
%                                    + sigma*theta_hat )
%
%     The 1 + Phi'*Phi denominator is normalised adaptation: it bounds the
%     drive when the states are large, so a big gust cannot produce a parameter
%     step large enough to destabilise the loop it is supposed to be helping.
%
%     with P_lyap from  Am'*P + P*Am = -Q_lyap.  Without the sigma term the
%     parameters can drift without bound whenever the error is driven by
%     something unmatched, such as sensor noise or the gust itself; sigma
%     modification trades exact asymptotic cancellation for bounded parameters,
%     which is the right trade on real hardware. A hard norm bound backs it up.
%
%   Output fields
%     Everything the LQG carries, plus Am, P_lyap, Gamma, sigma, theta_max.

    if nargin < 4, opts = struct(); end
    ns = P.n_states;

    % Adaptation rate. The drive term scales as Gamma*|Phi|*|B'*P_lyap*e|, and
    % with the Bryson-scaled weights that product is of order |Phi| ~ 5 and
    % |P_lyap*B| ~ 0.76. A rate of 5e3 therefore moves the parameters a
    % thousand times faster than the error they are meant to track: an early
    % run pinned theta at its norm bound within the first gust and made the
    % closed loop markedly worse than the LQG it was augmenting. Order 5 is the
    % scale that gives a parameter of order one over roughly a second.
    gamma     = get_opt(opts, 'gamma',     5);
    sigma     = get_opt(opts, 'sigma',     2e-1);
    theta_max = get_opt(opts, 'theta_max', 5);

    if isscalar(gamma)
        Gamma = gamma * eye(ns);
    else
        Gamma = diag(gamma(:));
    end

    % Reference model: the nominal closed loop on the structural states.
    Am = P.A - P.B * base.K_fb;
    if any(real(eig(Am)) >= 0)
        error('mrac_design:unstableReference', ...
              'The nominal closed loop is not stable; fix the LQG design first.');
    end

    % Lyapunov weight. Using the regulator's own state weight keeps the
    % adaptive error metric consistent with what the baseline is optimising.
    Q_lyap = base.Q(1:ns, 1:ns) + 1e-6*eye(ns);
    P_lyap = lyap_solve(Am, Q_lyap);

    ctrl           = base;
    ctrl.type      = 'mrac';
    ctrl.base      = base;
    ctrl.Am        = Am;
    ctrl.P_lyap    = P_lyap;
    ctrl.Gamma     = Gamma;
    ctrl.sigma     = sigma;
    ctrl.theta_max = theta_max;
    ctrl.n_struct  = ns;
    ctrl.PB        = P_lyap * P.B;        % precomputed, used every step
    ctrl.normalise = get_opt(opts, 'normalise', true);
end

% ------------------------------------------------------------------------
function v = get_opt(s, name, default)
    if isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end
