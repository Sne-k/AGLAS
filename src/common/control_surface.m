function cs = control_surface(cfg, fem, modes)
%CONTROL_SURFACE  Aerodynamic effectiveness of a trailing-edge control surface.
%
%   cs = CONTROL_SURFACE(cfg, fem, modes)
%
%   Adds the missing effector. Nothing in the open-loop model could produce a
%   control force, so gust load alleviation was not expressible: the plant had
%   a disturbance input and no control input.
%
%   Sectional aerodynamics (thin aerofoil theory with a plain flap)
%     With the chordwise coordinate written as x = (c/2)*(1 - cos(theta)), the
%     hinge at x_h maps to theta_h = acos(1 - 2*x_h/c), and a flap deflection
%     delta, positive trailing edge down, gives
%
%         dCl/ddelta      = 2*((pi - theta_h) + sin(theta_h))
%         dCm_c4/ddelta   = (1/4)*(sin(2*theta_h) - 2*sin(theta_h))
%
%     The first is the familiar result that a quarter-chord flap recovers about
%     60 % of the lift a whole-wing incidence change would give. The second is
%     negative: deflecting the flap down pitches the section nose down, which
%     twists the wing and partly cancels the lift it just created. That
%     aeroelastic loss of control effectiveness is the main reason an outboard
%     surface on a flexible wing underperforms its rigid estimate, so it has to
%     be in the model for the control design to mean anything.
%
%   Spanwise distribution
%     The surface covers a band of the semi-span. Its limits are snapped to
%     element boundaries so the work-equivalent nodal load vector stays exact
%     rather than needing partial-element integration.
%
%   Output fields
%     B_modal    generalised force per radian of deflection  [n_modes x 1]
%     dCl_ddelta, dCm_ddelta   sectional derivatives         [1/rad]
%     lift_per_rad, moment_per_rad   line loads              [N/m/rad, N.m/m/rad]
%     x_start, x_end, span_covered   actual extent after snapping   [m]
%     elements   indices of the elements the surface covers
%     effectiveness_ratio   dCl/ddelta divided by the 2D lift-curve slope

    c     = cfg.geom.chord;
    q_inf = 0.5 * cfg.flight.rho_air * cfg.flight.U_inf^2;

    % --- sectional derivatives -------------------------------------------
    Ef    = cfg.control.flap_chord_frac;
    x_h   = (1 - Ef) * c;                       % hinge, measured from the LE
    th_h  = acos(1 - 2*x_h/c);

    dCl = 2 * ((pi - th_h) + sin(th_h));
    dCm = 0.25 * (sin(2*th_h) - 2*sin(th_h));

    % --- line loads per radian -------------------------------------------
    lift_per_rad = q_inf * c * dCl;                       % [N/m/rad]
    % Moment about the elastic axis: the quarter-chord couple plus the flap
    % lift acting on its arm from the elastic axis.
    arm = (cfg.geom.ea_frac - 0.25) * c;
    moment_per_rad = q_inf * c^2 * dCm + lift_per_rad * arm;   % [N.m/m/rad]

    % --- spanwise extent, snapped to element boundaries -------------------
    L  = cfg.geom.L;
    dx = fem.dx;
    i_start = max(1,        round(cfg.control.span_start_frac * L / dx) + 1);
    i_end   = min(fem.n_el, round(cfg.control.span_end_frac   * L / dx));
    if i_end < i_start
        error('control_surface:emptySpan', ...
              'The control surface covers no whole element. Widen its span or refine the mesh.');
    end
    elements = i_start:i_end;

    x_start = (i_start - 1) * dx;
    x_end   = i_end * dx;

    % --- work-equivalent nodal load vector over those elements only -------
    f_full = zeros(fem.dpn * fem.n_nodes, 1);
    fe = lift_per_rad * dx * [1/2; dx/12; 1/2; -dx/12];
    if fem.torsion
        ft = moment_per_rad * dx * [1/2; 1/2];
    end
    for e = elements
        n1 = e; n2 = e + 1;
        bd = [fem.dpn*(n1-1) + (1:2), fem.dpn*(n2-1) + (1:2)];
        f_full(bd) = f_full(bd) + fe;
        if fem.torsion
            td = [fem.dpn*(n1-1) + 3, fem.dpn*(n2-1) + 3];
            f_full(td) = f_full(td) + ft;
        end
    end
    f_red = f_full(fem.free_dof);

    cs.B_modal            = modes.Phi.' * f_red;
    cs.f_nodal            = f_red;
    cs.dCl_ddelta         = dCl;
    cs.dCm_ddelta         = dCm;
    cs.lift_per_rad       = lift_per_rad;
    cs.moment_per_rad     = moment_per_rad;
    cs.arm                = arm;
    cs.x_start            = x_start;
    cs.x_end              = x_end;
    cs.span_covered       = x_end - x_start;
    cs.elements           = elements;
    cs.effectiveness_ratio = dCl / cfg.aero.cl_alpha_2d;
    cs.hinge_x            = x_h;
end
