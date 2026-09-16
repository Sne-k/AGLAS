function f = consistent_line_load(fem, q_lift, q_torque)
%CONSISTENT_LINE_LOAD  Work-equivalent nodal load vector for a line load.
%
%   f = CONSISTENT_LINE_LOAD(fem, q_lift)
%   f = CONSISTENT_LINE_LOAD(fem, q_lift, q_torque)
%
%   Converts a uniform distributed lift q_lift [N/m] (and, when the model has a
%   torsional DOF, a uniform distributed torque q_torque [N.m/m] about the
%   elastic axis) into the equivalent nodal load vector for the reduced
%   (clamped-root) system.
%
%   For a cubic Hermite beam element of length dx under a uniform transverse
%   load q, the work-equivalent element load vector is
%
%       f_e = q*dx * [ 1/2 ; dx/12 ; 1/2 ; -dx/12 ]
%
%   The dx/12 terms are nodal moments. They are not optional decoration: they
%   are what makes the discrete load do the same virtual work as the real
%   distributed load, and dropping them biases the deflection and, much more
%   strongly, the recovered bending moment.
%
%   For the linear torsion element the equivalent vector is
%
%       f_t = q_torque*dx * [ 1/2 ; 1/2 ]
%
%   The original code instead wrote force_distribution(1:2:end) = 1, i.e. a
%   unit point force on every translational DOF with no element length scaling
%   and no moment terms. That is not a discretisation of any distributed load:
%   it ignores dx entirely, so refining the mesh changes the total applied
%   force instead of converging.

    if nargin < 3 || isempty(q_torque)
        q_torque = 0;
    end

    dx  = fem.dx;
    dpn = fem.dpn;
    n_full = dpn * fem.n_nodes;
    f_full = zeros(n_full, 1);

    fe = q_lift * dx * [1/2; dx/12; 1/2; -dx/12];
    if fem.torsion
        ft = q_torque * dx * [1/2; 1/2];
    end

    for e = 1:fem.n_el
        n1 = e; n2 = e + 1;
        bd = [dpn*(n1-1) + (1:2), dpn*(n2-1) + (1:2)];
        f_full(bd) = f_full(bd) + fe;
        if fem.torsion
            td = [dpn*(n1-1) + 3, dpn*(n2-1) + 3];
            f_full(td) = f_full(td) + ft;
        end
    end

    f = f_full(fem.free_dof);
end
