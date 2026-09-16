function fem = beam_fem(cfg)
%BEAM_FEM  Assemble the cantilever wing finite element model.
%
%   fem = BEAM_FEM(cfg) builds the global mass and stiffness matrices for a
%   prismatic cantilever wing and applies the clamped-root boundary condition.
%
%   Element formulation
%     Bending  : 2-node Euler-Bernoulli beam, cubic Hermite interpolation,
%                consistent mass matrix. Exact for a prismatic beam.
%     Torsion  : 2-node Saint-Venant shaft, linear interpolation, consistent
%                mass matrix. Present only when cfg.fem.include_torsion is true.
%     Coupling : inertial only. Bending and torsion are elastically uncoupled
%                because twist is measured about the elastic axis, but they are
%                inertially coupled whenever the section CG is offset from that
%                axis. That static unbalance is what makes classical
%                bending-torsion flutter possible, and it is why the original
%                bending-only model could never produce a flutter speed.
%
%   Sign convention
%     w   vertical deflection, positive up                              [m]
%     th  bending rotation dw/dx                                        [rad]
%     ph  torsional twist about the elastic axis, positive nose-up      [rad]
%
%   With the CG a distance d aft of the elastic axis, a nose-up twist lowers
%   the CG, so the CG height is z_cg = w - d*ph and the kinetic energy carries
%   a negative cross term, -S*wdot*phdot with S = m*d.
%
%   DOF ordering is node-major: node n owns DOFs [w th ph] (or [w th]).
%
%   Output fields
%     M, K            reduced (constrained) global matrices
%     M_full, K_full  unconstrained global matrices
%     free_dof        indices of the retained DOFs in the full system
%     n_nodes, n_el, dx, dpn (DOFs per node)
%     x_nodes         nodal spanwise stations, including the root  [1 x n_nodes]
%     idx_w, idx_th, idx_ph  reduced-system indices of each DOF type
%
%   See also AGLAS_CONFIG, SECTION_PROPERTIES, MODAL_ANALYSIS.

    sec = cfg.section;
    L   = cfg.geom.L;
    n_el = cfg.fem.n_elements;
    dx  = L / n_el;
    n_nodes = n_el + 1;
    torsion = logical(cfg.fem.include_torsion);
    dpn = 2 + double(torsion);          % DOFs per node

    n_dof = dpn * n_nodes;
    K = zeros(n_dof);
    M = zeros(n_dof);

    % ------------------------------------------------- element sub-matrices
    EI = sec.EI;
    m  = sec.m_total;

    k_b = (EI/dx^3) * [ 12,     6*dx,   -12,     6*dx;
                        6*dx,   4*dx^2, -6*dx,   2*dx^2;
                       -12,    -6*dx,    12,    -6*dx;
                        6*dx,   2*dx^2, -6*dx,   4*dx^2 ];

    m_b = (m*dx/420) * [ 156,    22*dx,   54,    -13*dx;
                         22*dx,  4*dx^2,  13*dx, -3*dx^2;
                         54,     13*dx,   156,   -22*dx;
                        -13*dx, -3*dx^2, -22*dx,  4*dx^2 ];

    if torsion
        GJ = sec.GJ;
        It = sec.I_theta;
        S  = m * sec.d_ea_cg;           % static unbalance per unit span [kg]

        k_t = (GJ/dx)   * [ 1, -1; -1,  1 ];
        m_t = (It*dx/6) * [ 2,  1;  1,  2 ];

        % Inertial bending-torsion coupling: -S * integral( Nw' * Nph ) dx,
        % with Nw the cubic Hermite set and Nph the linear pair.
        m_bt = -S * dx * [  7/20,    3/20;
                            dx/20,   dx/30;
                            3/20,    7/20;
                           -dx/30,  -dx/20 ];
    end

    % ------------------------------------------------------------ assembly
    for e = 1:n_el
        n1 = e;  n2 = e + 1;
        bd = [dpn*(n1-1) + (1:2), dpn*(n2-1) + (1:2)];     % bending DOFs
        K(bd, bd) = K(bd, bd) + k_b;
        M(bd, bd) = M(bd, bd) + m_b;

        if torsion
            td = [dpn*(n1-1) + 3, dpn*(n2-1) + 3];         % torsion DOFs
            K(td, td) = K(td, td) + k_t;
            M(td, td) = M(td, td) + m_t;
            M(bd, td) = M(bd, td) + m_bt;
            M(td, bd) = M(td, bd) + m_bt.';
        end
    end

    % Enforce exact symmetry (guards against round-off asymmetry).
    K = (K + K.')/2;
    M = (M + M.')/2;

    % ------------------------------------------ clamped root at node 1, x = 0
    free_dof = (dpn + 1):n_dof;

    fem = struct();
    fem.K_full   = K;
    fem.M_full   = M;
    fem.K        = K(free_dof, free_dof);
    fem.M        = M(free_dof, free_dof);
    fem.free_dof = free_dof;
    fem.n_nodes  = n_nodes;
    fem.n_el     = n_el;
    fem.dx       = dx;
    fem.dpn      = dpn;
    fem.torsion  = torsion;
    fem.x_nodes  = linspace(0, L, n_nodes);
    fem.n_dof    = numel(free_dof);

    % Indices, within the reduced system, of each physical DOF type.
    % Reduced DOF j corresponds to full DOF j + dpn, i.e. nodes 2..n_nodes.
    fem.idx_w  = 1:dpn:fem.n_dof;
    fem.idx_th = 2:dpn:fem.n_dof;
    if torsion
        fem.idx_ph = 3:dpn:fem.n_dof;
    else
        fem.idx_ph = [];
    end
end
