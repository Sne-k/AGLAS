function A = span_integral_matrices(fem)
%SPAN_INTEGRAL_MATRICES  Span integrals of shape-function products.
%
%   A = SPAN_INTEGRAL_MATRICES(fem) returns the matrices needed to evaluate
%
%       A.ww(i,j) = integral over span of  phi_w,i(x) * phi_w,j(x)  dx
%       A.wp(i,j) = integral over span of  phi_w,i(x) * phi_phi,j(x) dx
%       A.pp(i,j) = integral over span of  phi_phi,i(x)*phi_phi,j(x) dx
%
%   as quadratic forms in the reduced DOF vector, so that for mode shapes Phi
%
%       integral(phi_w,i * phi_w,j) = Phi(:,i)' * A.ww * Phi(:,j)
%
%   These are exactly the consistent "mass" matrices of the element with unit
%   coefficients, so they are assembled from the same interpolation used by
%   BEAM_FEM. Strip-theory aerodynamics needs them to turn sectional lift and
%   moment coefficients into generalised aerodynamic forces.

    dx  = fem.dx;
    dpn = fem.dpn;
    n_full = dpn * fem.n_nodes;

    Aww = zeros(n_full);
    Awp = zeros(n_full);
    App = zeros(n_full);

    % integral of Nw' * Nw  (cubic Hermite, unit mass per length)
    e_ww = (dx/420) * [ 156,    22*dx,   54,    -13*dx;
                        22*dx,  4*dx^2,  13*dx, -3*dx^2;
                        54,     13*dx,   156,   -22*dx;
                       -13*dx, -3*dx^2, -22*dx,  4*dx^2 ];

    % integral of Nphi' * Nphi  (linear)
    e_pp = (dx/6) * [ 2, 1; 1, 2 ];

    % integral of Nw' * Nphi  (cubic Hermite x linear)
    e_wp = dx * [  7/20,    3/20;
                   dx/20,   dx/30;
                   3/20,    7/20;
                  -dx/30,  -dx/20 ];

    for e = 1:fem.n_el
        n1 = e; n2 = e + 1;
        bd = [dpn*(n1-1) + (1:2), dpn*(n2-1) + (1:2)];
        Aww(bd, bd) = Aww(bd, bd) + e_ww;
        if fem.torsion
            td = [dpn*(n1-1) + 3, dpn*(n2-1) + 3];
            App(td, td) = App(td, td) + e_pp;
            Awp(bd, td) = Awp(bd, td) + e_wp;
        end
    end

    fd = fem.free_dof;
    A.ww = Aww(fd, fd);
    A.wp = Awp(fd, fd);
    A.pp = App(fd, fd);
    A.pw = A.wp.';
end
