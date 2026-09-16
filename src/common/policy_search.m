function out = policy_search(cfg, P, base, t, wg, opts)
%POLICY_SEARCH  Tune the control policy directly against the simulated plant.
%
%   out = POLICY_SEARCH(cfg, P, base, t, wg)
%   out = POLICY_SEARCH(cfg, P, base, t, wg, opts)
%
%   opts fields:
%     .max_eval   evaluation budget, default 200
%     .h0         initial step in the scale parameters, default 0.4
%     .h_tol      convergence step size, default 0.01
%     .objective  'peak' (default) or 'rms' root bending moment
%     .verbose    print progress, default true
%
%   Why this is worth doing at all
%     LQG is exactly optimal, so at first sight there is nothing to search for.
%     It is optimal for a specific problem though, and that problem is not the
%     one being solved here. Two mismatches:
%
%       1. It minimises an integral quadratic cost. The quantity the structure
%          is actually sized by is the single largest root bending moment in
%          the event. Those are different objectives and they have different
%          optima.
%       2. It is derived for a linear plant. The real loop has a deflection
%          limit and a slew limit, and LQR has no representation of either. The
%          tuning sweep in the design log shows how badly that can bite: at
%          r_command = 0.001 the linear-optimal design saturates so hard that
%          the peak load becomes six times worse than no control at all.
%
%     So there is a genuine gap, and it is exactly the gap a learned policy can
%     close: optimise the real objective on the real, saturated plant.
%
%   What this is, and what it is not
%     This is direct policy search, the derivative-free end of reinforcement
%     learning: a parameterised policy, a black-box return, and a search over
%     parameters. It is not deep reinforcement learning and is not presented as
%     such. Deep RL would be the wrong tool here and would be worth saying so
%     plainly: the plant is known, low order and linear apart from two
%     saturations, so learning a value function from scratch would burn orders
%     of magnitude more evaluations to rediscover a controller that a Riccati
%     equation already gives in closed form. The useful move is to start from
%     that controller and search only over what it cannot model.
%
%   Policy class
%     Four multiplicative scales on blocks of the LQG gain, starting at the
%     LQG solution itself:
%
%       u = -(s_q*K_q*qhat + s_qd*K_qd*qdothat + s_a*K_a*actuator + s_ff*K_w*what)
%
%     Starting at [1 1 1 1] guarantees the search can never do worse than LQG,
%     because that point is in the search space and is evaluated first. A
%     larger policy class would fit the objective better and generalise worse;
%     four parameters against a single gust is already close to the limit of
%     what can be claimed without a separate validation disturbance.
%
%   Optimiser
%     Compass search: try plus and minus one step in each coordinate, move to
%     the best improvement, halve the step when none improves. Deterministic,
%     no toolbox, no gradients, and it cannot diverge. fminsearch would also
%     work but its simplex wanders through parameter combinations that
%     destabilise the loop, and every such evaluation costs a full simulation.
%
%   Output fields
%     scales, J, J0, improvement_pct, n_eval, history, ctrl (tuned controller)

    if nargin < 6, opts = struct(); end
    max_eval  = getf(opts, 'max_eval', 200);
    h         = getf(opts, 'h0',       0.4);
    h_tol     = getf(opts, 'h_tol',    0.01);
    objective = getf(opts, 'objective','peak');
    verbose   = getf(opts, 'verbose',  true);

    nm = P.n_modes;
    ns = P.n_states;
    na = base.n_aug;

    % Index blocks of the augmented gain: modal position, modal rate,
    % actuator states, gust feedforward.
    blocks = { 1:nm, nm+(1:nm), (2*nm+1):ns, na };

    K0 = base.K;

    J = @(s) evaluate(s);
    x = [1 1 1 1];
    fx = J(x);
    J0 = fx;
    n_eval = 1;
    history = [x, fx];

    if verbose
        fprintf('policy search, objective = %s\n', objective);
        fprintf('  start (LQG) : %.4f kN m\n', fx/1e3);
    end

    while h > h_tol && n_eval < max_eval
        improved = false;
        for i = 1:numel(x)
            for sgn = [1 -1]
                xt = x;
                xt(i) = max(xt(i) + sgn*h, 0);
                ft = J(xt);
                n_eval = n_eval + 1;
                history(end+1, :) = [xt, ft]; %#ok<AGROW>
                if ft < fx
                    x = xt; fx = ft; improved = true;
                    break;
                end
                if n_eval >= max_eval, break; end
            end
            if n_eval >= max_eval, break; end
        end
        if ~improved
            h = h/2;
        end
    end

    if verbose
        fprintf('  tuned       : %.4f kN m after %d evaluations\n', fx/1e3, n_eval);
        fprintf('  scales      : [%.3f %.3f %.3f %.3f]\n', x);
        fprintf('  improvement : %.2f %%\n', 100*(1 - fx/J0));
    end

    out.scales          = x;
    out.J               = fx;
    out.J0              = J0;
    out.improvement_pct = 100*(1 - fx/J0);
    out.n_eval          = n_eval;
    out.history         = history;
    out.objective       = objective;
    out.ctrl            = apply_scales(base, K0, blocks, x);

    % ------------------------------------------------- nested evaluation
    function val = evaluate(s)
        c = apply_scales(base, K0, blocks, s);
        try
            r = closed_loop_sim(cfg, P, c, t, wg);
            if strcmpi(objective, 'rms')
                val = sqrt(mean(r.moment.^2));
            else
                val = max(abs(r.moment));
            end
            if ~isfinite(val), val = Inf; end
        catch
            % A destabilising parameter set simply scores infinitely badly.
            val = Inf;
        end
    end
end

% ------------------------------------------------------------------------
function c = apply_scales(base, K0, blocks, s)
    c = base;
    K = K0;
    for i = 1:numel(blocks)
        K(blocks{i}) = K0(blocks{i}) * s(i);
    end
    c.K = K;
    c.K_fb = K(1:end-1);
    c.K_ff = K(end);
end

% ------------------------------------------------------------------------
function v = getf(s, name, default)
    if isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end
