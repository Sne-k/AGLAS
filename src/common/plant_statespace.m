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
    % it omits the C(k) phase lag and the apparent mass, and its stability
    % boundary is wrong: it loses stability at about 216 m/s, against the
    % 314.8 m/s p-k flutter speed. (An earlier version reported 377 m/s; that
    % figure came from the spurious torsion damping described below.) Use it
    % for control design well below that boundary. cfg.aero.plant_model =
    % 'unsteady' adds the circulation lag and reproduces the p-k boundary to
    % about 1 %; FLUTTER_PK remains the reference for the margin itself.
    %
    % The quasi-steady limit of the same Theodorsen strip theory used for
    % flutter is C(k) = 1, split into an aerodynamic stiffness (in phase with
    % displacement) and an aerodynamic damping (in phase with velocity) with
    % i*k <-> (b/U)*d/dt, plus the non-circulatory damping.
    %
    % This used to be obtained by evaluating AERO_GENERALISED at k = 1e-5 and
    % dividing the imaginary part by k. That is not a limit: Im C(k) ~ k*ln(k)
    % near zero, so the division leaves a term in ln(k0) that grows without
    % bound as k0 shrinks. At k0 = 1e-5 it overstated the torsion-mode
    % aerodynamic damping by a factor of about 27.
    G   = modal_span_integrals(fem, modes.Phi);
    rho = cfg.flight.rho_air;
    U   = cfg.flight.U_inf;
    b   = cfg.aero.semi_chord;
    a   = cfg.aero.a_ea;
    Kc   = rho * U * b * cfg.aero.cl_alpha_2d;
    Wl   = G.wp + b*(a + 0.5)*G.pp;
    E_q  = U * Wl;
    E_qd = -(G.ww + b*(a + 0.5)*G.pw) + b*(0.5 - a)*Wl;
    C_nc = pi*rho*U*b^2 * (G.wp - b*(0.5 - a)*G.pp);
    K_aero = Kc * E_q;
    C_aero = Kc * E_qd + C_nc;

    % ------------------------------------------------------ actuator block
    wa = 2*pi*cfg.control.act_freq_hz;
    za = cfg.control.act_zeta;

    if strcmpi(get_model(cfg), 'unsteady')
        P = unsteady_plant(cfg, modes, G, g_delta, g_gust, wa, za, ...
                           C_moment, C_tip, C_strain, kappa_root_modal);
        P.K_aero = K_aero;
        P.C_aero = C_aero;
        return;
    end

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
    P.model    = 'quasisteady';
    P.idx.aero = [];
    P.meas_names = {'root strain', 'tip acceleration'};
    P.perf_names = {'root bending moment', 'tip deflection'};
end

% ------------------------------------------------------------------------
function m = get_model(cfg)
    m = 'quasisteady';
    if isfield(cfg, 'aero') && isfield(cfg.aero, 'plant_model')
        m = cfg.aero.plant_model;
    end
end

% ------------------------------------------------------------------------
function P = unsteady_plant(cfg, modes, G, g_delta, g_gust, wa, za, ...
                            C_moment, C_tip, C_strain, kappa_root_modal)
%UNSTEADY_PLANT  Time-domain Theodorsen strip theory with indicial lags.
%
%   State vector, n modes:
%
%     x = [ q ; qdot ; delta ; deltadot ; z1 ; z2 ; zd ; zg ]
%
%   z1, z2 (n each) are the circulatory lag states of the wing's own motion,
%   zd (2) the lag states of the flap lift and zg (2) the Kussner lag of the
%   gust. Each lag realises an indicial function
%
%       phi(s) = 1 - A1*exp(-b1*s) - A2*exp(-b2*s),   s = U*t/b,
%
%   whose transfer function phi(0) + sum A_i*beta_i/(p + beta_i), with
%   beta_i = b_i*U/b, is the time-domain counterpart of C(k): it is 1 at
%   p = 0 and phi(0) as p -> infinity. The circulatory force on the wing is
%   that filter applied to the three-quarter-chord downwash projected on the
%   modes, which is exactly the C(k) part of AERO_GENERALISED rewritten with
%   i*k <-> (b/U)*d/dt, so the plant agrees with the p-k flutter model up to
%   the accuracy of the exponential fit.
%
%   The non-circulatory (apparent mass) forces are carried exactly: an added
%   mass matrix and an added damping matrix, both incompressible. They are a
%   small share of the total at the reduced frequencies of interest here.
%
%   With cfg.aero.compressible the lift slopes become cl_alpha/beta and the
%   lag rates b_i*beta^2, beta = sqrt(1 - M^2) (Leishman & Nguyen 1990).

    nm    = modes.n_modes;
    omega = modes.omega(:);
    zeta  = cfg.damping.zeta;

    rho = cfg.flight.rho_air;
    U   = cfg.flight.U_inf;
    b   = cfg.aero.semi_chord;
    a   = cfg.aero.a_ea;
    cla = cfg.aero.cl_alpha_2d;

    switch lower(cfg.aero.indicial)
        case 'jones',    A_ind = [0.165, 0.335]; b_ind = [0.0455, 0.300];
        case 'leishman', A_ind = [0.300, 0.700]; b_ind = [0.1400, 0.530];
        otherwise
            error('plant_statespace:indicial', ...
                  'cfg.aero.indicial must be ''jones'' or ''leishman''.');
    end

    beta_pg = 1;
    if cfg.aero.compressible
        M = U / cfg.flight.a_sound;
        if M >= 0.9
            error('plant_statespace:mach', ...
                  'Mach %.2f: Prandtl-Glauert compressibility is not valid this close to 1.', M);
        end
        beta_pg = sqrt(1 - M^2);
    end

    lam = b_ind * beta_pg^2 * U / b;          % lag rates [1/s]
    d0  = 1 - sum(A_ind);                     % phi(0), direct part
    Ab  = A_ind .* lam;

    % Kussner function for the gust, psi(s) = 1 - 0.5 e^(-0.13 s) - 0.5 e^(-s),
    % the same fit the analysis pipeline convolves in gust_modal_force.
    lam_g = [0.13, 1.0] * beta_pg^2 * U / b;
    Ab_g  = [0.5, 0.5] .* lam_g;

    % Circulatory force = Kc * filter( E_q*q + E_qd*qdot ).
    % Kc*E_q and Kc*E_qd are the C = 1 parts of AERO_GENERALISED.
    Kc   = rho * U * b * cla / beta_pg;
    Wl   = G.wp + b*(a + 0.5)*G.pp;           % lift weighting incl. moment arm
    E_q  = U * Wl;
    E_qd = -(G.ww + b*(a + 0.5)*G.pw) + b*(0.5 - a)*Wl;

    % Non-circulatory: apparent mass and damping (the k^2 and i*k terms of
    % AERO_GENERALISED with C = 0).
    M_app = pi*rho*b^2 * (G.ww + a*b*(G.wp + G.pw) + b^2*(1/8 + a^2)*G.pp);
    C_nc  = pi*rho*U*b^2 * (G.wp - b*(0.5 - a)*G.pp);

    g_delta = g_delta / beta_pg;
    g_gust  = g_gust  / beta_pg;

    Mi = inv(eye(nm) + M_app);

    Z  = zeros(nm);  I = eye(nm);
    zc = zeros(nm, 1);  z2 = zeros(nm, 2);  zr = zeros(1, nm);

    % qddot = Mi * ( Kq*q + Kqd*qdot + Kc*(Ab(1)*z1 + Ab(2)*z2)
    %                + g_delta*(d0*delta + Ab*zd) + g_gust*(Ab_g*zg) )
    Kq  = -diag(omega.^2)     + Kc*d0*E_q;
    Kqd = -diag(2*zeta*omega) + C_nc + Kc*d0*E_qd;

    row_q   = Mi*Kq;
    row_qd  = Mi*Kqd;
    row_del = Mi*g_delta*d0;
    row_z1  = Mi*Kc*Ab(1);
    row_z2  = Mi*Kc*Ab(2);
    row_zd  = Mi*g_delta*Ab;
    row_zg  = Mi*g_gust*Ab_g;

    A = [ Z,      I,      zc,      zc,   Z,              Z,              z2,              z2;
          row_q,  row_qd, row_del, zc,   row_z1,         row_z2,         row_zd,          row_zg;
          zr,     zr,     0,       1,    zr,             zr,             [0 0],           [0 0];
          zr,     zr,    -wa^2,   -2*za*wa, zr,          zr,             [0 0],           [0 0];
          E_q,    E_qd,   zc,      zc,  -lam(1)*I,       Z,              z2,              z2;
          E_q,    E_qd,   zc,      zc,   Z,             -lam(2)*I,       z2,              z2;
          zeros(2,nm), zeros(2,nm), [1;1], [0;0], zeros(2,nm), zeros(2,nm), -diag(lam),  zeros(2);
          zeros(2,nm), zeros(2,nm), [0;0], [0;0], zeros(2,nm), zeros(2,nm), zeros(2),   -diag(lam_g) ];

    ns = size(A, 1);
    B  = zeros(ns, 1);  B(2*nm + 2) = wa^2;
    Bw = zeros(ns, 1);  Bw(end-1:end) = 1;

    % Measurements: root strain from q; tip acceleration from the qddot row.
    % The gust reaches the forces only through the Kussner lag (psi(0) = 0),
    % so there is no gust feedthrough to the accelerometer.
    C_eps = [C_strain, zeros(1, ns - nm)];
    C_acc = C_tip * A(nm + (1:nm), :);
    C  = [C_eps; C_acc];
    D  = [0; 0];
    Dw = [0; 0];

    Cz  = [C_moment, zeros(1, ns - nm);
           C_tip,    zeros(1, ns - nm)];
    Dz  = [0; 0];
    Dzw = [0; 0];

    P.A = A;  P.B = B;  P.Bw = Bw;
    P.C = C;  P.D = D;  P.Dw = Dw;
    P.Cz = Cz; P.Dz = Dz; P.Dzw = Dzw;

    P.n_modes   = nm;
    P.n_states  = ns;
    P.idx.q     = 1:nm;
    P.idx.qdot  = nm + (1:nm);
    P.idx.delta = 2*nm + 1;
    P.idx.ddot  = 2*nm + 2;
    P.idx.aero  = (2*nm + 3):ns;

    P.C_moment = C_moment;
    P.C_tip    = C_tip;
    P.C_strain = C_strain;
    P.kappa_root_modal = kappa_root_modal;
    P.g_delta  = g_delta;
    P.g_gust   = g_gust;
    P.omega    = omega;
    P.zeta     = zeta;
    P.act      = struct('wn', wa, 'zeta', za);
    P.model    = 'unsteady';
    P.indicial = struct('A', A_ind, 'b', b_ind, 'beta', beta_pg, 'rates', lam);
    P.M_app    = M_app;
    P.meas_names = {'root strain', 'tip acceleration'};
    P.perf_names = {'root bending moment', 'tip deflection'};
end
