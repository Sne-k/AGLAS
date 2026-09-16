function [F_fun, info] = gust_modal_force(cfg, fem, modes, t_grid, wg, unsteady)
%GUST_MODAL_FORCE  Generalised modal force produced by a vertical gust.
%
%   [F_fun, info] = GUST_MODAL_FORCE(cfg, fem, modes, t_grid, wg)
%   [F_fun, info] = GUST_MODAL_FORCE(..., unsteady)
%
%   Converts a vertical gust velocity history wg(t) into a function handle
%   F_fun(t) returning the generalised force vector [n_modes x 1] for the modal
%   equations of motion.
%
%   Aerodynamics (strip theory, unswept prismatic wing)
%     A vertical gust w_g changes the local incidence by dalpha = w_g/U, so the
%     incremental lift per unit span is
%
%         l'(t) = q_inf * c * cl_alpha * dalpha
%               = 0.5 * rho * U^2 * c * cl_alpha * (w_g/U)
%               = 0.5 * rho * U * c * cl_alpha * w_g(t)          [N/m]
%
%     The lift acts at the quarter-chord aerodynamic centre. With the elastic
%     axis further aft, it also applies a nose-up torque about that axis:
%
%         m'(t) = l'(t) * (x_ea - x_ac)                          [N.m/m]
%
%   Unsteady option ('kussner', the default)
%     A wing does not develop its full gust lift instantaneously. The Kussner
%     function psi(s) gives the indicial lift growth as the wing penetrates a
%     sharp-edged gust, with s = 2*U*t/c the distance travelled in semi-chords.
%     Using R.T. Jones's two-exponential fit
%
%         psi(s) = 1 - 0.5*exp(-0.13*s) - 0.5*exp(-s)
%
%     the response to an arbitrary gust follows by Duhamel superposition,
%
%         l'(t) = 0.5*rho*U*c*cl_alpha * integral_0^s psi'(s-sigma) w_g(sigma) dsigma
%
%     Pass 'quasisteady' to skip the lag and use l' proportional to w_g(t).
%
%   What the original code did
%     gust_force_function = @(t) (t<=t_gust) * (0.5*V_gust*(1-cos(pi*t/t_gust))) * F_generalized
%     fed the gust velocity, in m/s, straight in as a force in newtons. There is
%     no dynamic pressure, no chord, no lift-curve slope and no air density
%     anywhere in it, so the result is dimensionally meaningless and cannot be
%     compared against a yield stress. It also has no torsional component at all.
%
%   Outputs
%     F_fun  function handle, F_fun(t) -> [n_modes x 1] generalised force
%     info   struct: q_line (line load history), F_modal (on t_grid),
%            shape vector, peak incidence change and validity flags

    if nargin < 6 || isempty(unsteady)
        unsteady = 'kussner';
    end

    U       = cfg.flight.U_inf;
    rho     = cfg.flight.rho_air;
    chord   = cfg.geom.chord;
    cl_a    = cfg.aero.cl_alpha_2d;

    t_grid = t_grid(:).';
    wg     = wg(:).';

    % --- circulatory lift growth -----------------------------------------
    switch lower(unsteady)
        case 'quasisteady'
            wg_eff = wg;
        case 'kussner'
            wg_eff = kussner_convolve(t_grid, wg, U, chord);
        otherwise
            error('gust_modal_force:unknownModel', ...
                  'unsteady must be ''kussner'' or ''quasisteady''.');
    end

    q_line = 0.5 * rho * U * chord * cl_a * wg_eff;      % [N/m]

    % --- torque about the elastic axis -----------------------------------
    x_ac = 0.25 * chord;
    x_ea = cfg.geom.ea_frac * chord;
    arm  = x_ea - x_ac;                    % >0 => lift ahead of EA => nose-up
    m_line = q_line * arm;                 % [N.m/m]

    % --- spatial distribution (uniform gust => constant shape) -----------
    shape_lift   = consistent_line_load(fem, 1, 0);
    if fem.torsion
        shape_torque = consistent_line_load(fem, 0, 1);
    else
        shape_torque = zeros(size(shape_lift));
    end

    Phi = modes.Phi;
    g_lift   = Phi.' * shape_lift;         % [n_modes x 1] per unit N/m
    g_torque = Phi.' * shape_torque;       % [n_modes x 1] per unit N.m/m

    F_modal = g_lift * q_line + g_torque * m_line;   % [n_modes x n_t]

    % Linear interpolation inside the grid, held constant outside it.
    % The transposed copy and the uniform-grid test are hoisted out of the
    % closure: this handle is evaluated once per ODE stage, tens of thousands
    % of times per run, so allocating or re-deriving anything inside it
    % dominates the cost of the whole simulation.
    Fm_t = F_modal.';
    dt_u = NaN;
    if numel(t_grid) > 2
        d = diff(t_grid);
        if max(abs(d - d(1))) <= 1e-12 * max(1, abs(d(1)))
            dt_u = d(1);
        end
    end
    F_fun = @(t) interp_force(t, t_grid, F_modal, Fm_t, dt_u);

    % --- diagnostics ------------------------------------------------------
    info.q_line       = q_line;
    info.m_line       = m_line;
    info.F_modal      = F_modal;
    info.t_grid       = t_grid;
    info.g_lift       = g_lift;
    info.g_torque     = g_torque;
    info.arm          = arm;
    info.wg_effective = wg_eff;
    info.d_alpha_deg  = atan2(max(abs(wg)), U) * 180/pi;
    info.peak_lift_N  = max(abs(q_line)) * cfg.geom.L;
    info.linear_aero_ok = info.d_alpha_deg <= cfg.limits.max_delta_alpha;
end

% ------------------------------------------------------------------------
function F = interp_force(t, tg, Fm, Fm_t, dt_u)
%INTERP_FORCE  Linear interpolation of the generalised force at time t.
%   On a uniform grid the bracketing index is found by arithmetic, which is
%   O(1); otherwise fall back to interp1 on the pre-transposed matrix.

    if t <= tg(1)
        F = Fm(:,1);
        return;
    elseif t >= tg(end)
        F = Fm(:,end);
        return;
    end

    if isfinite(dt_u)
        pos = (t - tg(1)) / dt_u;
        i0  = floor(pos) + 1;
        i0  = min(max(i0, 1), numel(tg) - 1);
        w   = (t - tg(i0)) / dt_u;
        F   = (1-w)*Fm(:,i0) + w*Fm(:,i0+1);
    else
        F = interp1(tg, Fm_t, t, 'linear').';
    end
end

% ------------------------------------------------------------------------
function wg_eff = kussner_convolve(t, wg, U, chord)
%KUSSNER_CONVOLVE  Exact Duhamel integral of the gust history with psi'(s).
%
%   Reduced time s = 2*U*t/chord (semi-chords travelled). With Jones's fit
%
%       psi'(s) = 0.065*exp(-0.13*s) + 0.5*exp(-s)
%
%   each exponential term y_j(s) = a_j * integral_0^s exp(-b_j*(s-sig))*w(sig) dsig
%   obeys  dy_j/ds = -b_j*y_j + a_j*w,  so the convolution can be advanced with
%   a recurrence that is exact whenever w is piecewise linear between samples.
%
%   This replaces a direct trapezoidal quadrature of the convolution integral.
%   That approach is O(n^2) and, more importantly, inaccurate: the fast
%   exp(-s) term decays over one semi-chord, so a grid coarser than that
%   over-integrates a convex kernel and the effective gust can exceed the true
%   one. The recurrence below is O(n) and grid-independent.

    s = 2 * U * t / chord;
    n = numel(s);
    wg_eff = zeros(1, n);
    if n < 2
        return;
    end

    a = [0.065, 0.5];
    b = [0.13,  1.0];
    y = [0, 0];

    for i = 2:n
        h = s(i) - s(i-1);
        if h <= 0
            wg_eff(i) = sum(y);
            continue;
        end
        w0 = wg(i-1);
        dw = wg(i) - w0;
        for j = 1:2
            E  = exp(-b(j)*h);
            I1 = (1 - E) / b(j);
            I2 = h*I1 - (1 - E*(1 + b(j)*h)) / b(j)^2;
            y(j) = E*y(j) + a(j) * (w0*I1 + (dw/h)*I2);
        end
        wg_eff(i) = sum(y);
    end
end
