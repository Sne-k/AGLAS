function P = plant_statespace(cfg, fem, modes, cs)
%PLANT_STATESPACE  Aeroservoelastic state-space model for control design.
%
%   P = PLANT_STATESPACE(cfg, fem, modes, cs)
%
%   Assembles the linear plant the gust load alleviation controller is designed
%   against. State vector, with n modes retained:
%
%       x = [ q ; qdot ; delta ; deltadot ]        (2n + 2 states)
%
%   q are the mass-normalised modal coordinates, delta the actual control
%   surface deflection. The actuator is a second-order lag driven by the
%   commanded deflection, so the controller cannot assume its command appears
%   on the surface instantly. That matters: an alleviation law that looks
%   excellent against an ideal actuator usually loses most of its benefit once
%   a realistic bandwidth and rate limit are in the loop.
%
%       xdot = A*x + B*delta_cmd + Bw*w_gust
%
%   Because the modes are mass-normalised, the structural block is exactly
%   [0 I; -diag(omega^2) -diag(2*zeta*omega)] with no generalised mass to carry.
%
%   Outputs
%     Measurements y (what the estimator actually sees):
%       1  root bending strain, as a strain gauge on the outer fibre
%       2  tip vertical acceleration
%     Performance outputs z (what the cost function penalises):
%       1  root bending moment
%       2  tip deflection
%
%   Tip acceleration is a function of qddot, so it carries direct feedthrough
%   from both the control command and the gust. Those feedthrough terms are
%   returned and must not be dropped; ignoring them biases the estimator.
%
%   Output fields
%     A, B, Bw                 state, control and disturbance matrices
%     C, D, Dw                 measurement matrices with feedthrough
%     Cz, Dz, Dzw              performance matrices
%     n_modes, n_states, idx   state bookkeeping
%     C_moment, C_tip, C_strain  modal-to-physical row vectors
%     kappa_root_modal         root curvature per unit modal coordinate

    nm    = modes.n_modes;
    omega = modes.omega(:);
    zeta  = cfg.damping.zeta;

    % --------------------------------------------- root curvature per mode
    % Exact Hermite second derivative at the start of element 1, which is the
    % clamped root. Uses the rotation DOFs as well as the translations.
    dx  = fem.dx;
    dpn = fem.dpn;
    B0  = [-6/dx^2, -4/dx, 6/dx^2, -2/dx];

    Phi_full = zeros(dpn*fem.n_nodes, nm);
    Phi_full(fem.free_dof, :) = modes.Phi;
    elem1_dofs = [1, 2, dpn+1, dpn+2];
    kappa_root_modal = B0 * Phi_full(elem1_dofs, :);        % [1 x nm]

    C_moment = cfg.section.EI      * kappa_root_modal;      % root moment
    C_strain = cfg.section.c_outer * kappa_root_modal;      % outer-fibre strain
    C_tip    = modes.Phi(fem.idx_w(end), :);                % tip deflection

    % ------------------------------------------------------ input vectors
    % Control: generalised force per radian of surface deflection.
    g_delta = cs.B_modal(:);

    % Gust: quasi-steady strip theory, generalised force per m/s of gust.
    k_gust   = 0.5 * cfg.flight.rho_air * cfg.flight.U_inf * ...
               cfg.geom.chord * cfg.aero.cl_alpha_2d;
    arm      = (cfg.geom.ea_frac - 0.25) * cfg.geom.chord;
    shape_L  = consistent_line_load(fem, 1, 0);
    if fem.torsion
        shape_M = consistent_line_load(fem, 0, 1);
    else
        shape_M = zeros(size(shape_L));
    end
    g_gust = k_gust * (modes.Phi.' * shape_L) + ...
             k_gust * arm * (modes.Phi.' * shape_M);
    g_gust = g_gust(:);

    % ------------------------------------- quasi-steady aeroelastic feedback
    % The wing's own motion generates aerodynamic forces: plunge velocity
    % changes local incidence by -wdot/U, and twist changes it by phi. Without
    % these terms the plant is a structure in still air with a gust bolted on,
    % its poles do not move with airspeed, and a controller could never be
    % tested off its design condition in any meaningful way. At 50 m/s the
    % aerodynamic damping raises the first bending damping ratio from the 0.02
    % structural value to about 0.30, which is the dominant effect and matches
    % a hand calculation of 0.5*rho*U*c*cl_alpha acting on the mode shape.
    %
    % A caution on what this model can and cannot predict. Being quasi-steady,
    % it omits the C(k) phase lag, which is the destabilising mechanism behind
    % classical flutter. It therefore loses stability at about 377 m/s in a
    % 6.5 Hz oscillatory mode, whereas the unsteady p-k solution in c4p puts
    % flutter at 314.8 m/s and 13.1 Hz. Quasi-steady theory over-predicting the
    % flutter speed is the expected result, not a discrepancy to reconcile.
    % This plant is for control design well below that boundary; use FLUTTER_PK
    % for anything touching the stability margin itself.
    %
    % The k -> 0 limit of the same Theodorsen strip theory used for flutter is
    % taken here, split into an aerodynamic stiffness (in phase with
    % displacement) and an aerodynamic damping (in phase with velocity), using
    % the identity i*k <-> (b/U)*d/dt.
    G  = modal_span_integrals(fem, modes.Phi);
    k0 = 1e-5;                              % numerically safe stand-in for k -> 0
    Q0 = aero_generalised(cfg, G, k0, cfg.flight.U_inf);
    b  = cfg.aero.semi_chord;
    K_aero = real(Q0);
    C_aero = (b / (k0 * cfg.flight.U_inf)) * imag(Q0);

    % ------------------------------------------------------ actuator block
    wa = 2*pi*cfg.control.act_freq_hz;
    za = cfg.control.act_zeta;

    % ------------------------------------------------------------ assembly
    Z  = zeros(nm);
    I  = eye(nm);
    zc = zeros(nm, 1);
    zr = zeros(1, nm);

    A_qq  = -(diag(omega.^2)      - K_aero);
    A_qqd = -(diag(2*zeta*omega)  - C_aero);

    A = [  Z,      I,      zc,      zc;
           A_qq,   A_qqd,  g_delta, zc;
           zr,     zr,     0,       1;
           zr,     zr,    -wa^2,   -2*za*wa ];

    B  = [zc; zc; 0; wa^2];
    Bw = [zc; g_gust; 0; 0];

    n_states = 2*nm + 2;

    % --------------------------------------------------------- measurements
    % Row 1: root strain, a function of q alone.
    C_eps = [C_strain, zr, 0, 0];

    % Row 2: tip acceleration = C_tip * qddot, and qddot is the second block
    % row of the state equation, so the feedthrough terms come with it.
    row_qddot_q     = C_tip * A_qq;
    row_qddot_qdot  = C_tip * A_qqd;
    row_qddot_delta =  C_tip * g_delta;
    C_acc = [row_qddot_q, row_qddot_qdot, row_qddot_delta, 0];
    D_acc = 0;                                % delta_cmd enters only via the actuator
    Dw_acc = C_tip * g_gust;                  % gust accelerates the tip at once

    C  = [C_eps; C_acc];
    D  = [0; D_acc];
    Dw = [0; Dw_acc];

    % --------------------------------------------------- performance outputs
    Cz  = [ C_moment, zr, 0, 0;
            C_tip,    zr, 0, 0 ];
    Dz  = [0; 0];
    Dzw = [0; 0];

    P.A = A;  P.B = B;  P.Bw = Bw;
    P.C = C;  P.D = D;  P.Dw = Dw;
    P.Cz = Cz; P.Dz = Dz; P.Dzw = Dzw;

    P.n_modes   = nm;
    P.n_states  = n_states;
    P.idx.q     = 1:nm;
    P.idx.qdot  = nm + (1:nm);
    P.idx.delta = 2*nm + 1;
    P.idx.ddot  = 2*nm + 2;

    P.C_moment = C_moment;
    P.C_tip    = C_tip;
    P.C_strain = C_strain;
    P.kappa_root_modal = kappa_root_modal;
    P.g_delta  = g_delta;
    P.g_gust   = g_gust;
    P.omega    = omega;
    P.zeta     = zeta;
    P.act      = struct('wn', wa, 'zeta', za);
    P.K_aero   = K_aero;
    P.C_aero   = C_aero;
    P.meas_names = {'root strain', 'tip acceleration'};
    P.perf_names = {'root bending moment', 'tip deflection'};
end
