function [t_sol, q_sol, U_dof, gust_info] = solve_modal_response(cfg, fem, md, t_grid, wg, model, max_step)
%SOLVE_MODAL_RESPONSE  Integrate the modal equations under a gust.
%
%   [t_sol, q_sol, U_dof, gust_info] = SOLVE_MODAL_RESPONSE(cfg, fem, md, ...
%                                        t_grid, wg, model, max_step)
%
%   Solves, for mass-normalised mode shapes,
%
%       qddot_i + 2*zeta_i*omega_i*qdot_i + omega_i^2*q_i = F_i(t)
%
%   with F from GUST_MODAL_FORCE, and reconstructs the physical DOF history.
%
%   Mass normalisation is what makes this form valid. With peak-normalised
%   mode shapes the generalised mass is not unity, so neither diag(2*zeta*omega)
%   nor diag(omega^2) is the correct modal damping or stiffness.
%
%   Inputs
%     md        modal struct trimmed to the modes to retain
%     t_grid    output time grid                                    [s]
%     wg        gust velocity on t_grid                             [m/s]
%     model     'kussner' (default) or 'quasisteady'
%     max_step  optional ode45 MaxStep; defaults to one twentieth of the
%               shortest retained modal period, which keeps the solver from
%               stepping over the highest frequency it is meant to resolve
%
%   Outputs
%     t_sol     solution times                                      [s]
%     q_sol     generalised coordinates      [n_time x n_modes]
%     U_dof     reduced physical DOFs        [n_time x n_dof]
%     gust_info diagnostics from GUST_MODAL_FORCE

    if nargin < 6 || isempty(model),    model = 'kussner'; end
    if nargin < 7 || isempty(max_step)
        max_step = 1 / (20 * max(md.freq_hz));
    end

    n_modes = md.n_modes;
    [F_fun, gust_info] = gust_modal_force(cfg, fem, md, t_grid, wg, model);

    C_modal = 2 * cfg.damping.zeta * md.omega(:);
    K_modal = md.omega(:).^2;

    odefun = @(t, x) [ x(n_modes+1:end);
                       F_fun(t) - C_modal.*x(n_modes+1:end) - K_modal.*x(1:n_modes) ];

    opts = odeset('RelTol', cfg.time.rel_tol, ...
                  'AbsTol', cfg.time.abs_tol, ...
                  'MaxStep', max_step);

    [t_sol, x_sol] = ode45(odefun, t_grid, zeros(2*n_modes,1), opts);

    t_sol = t_sol(:);
    q_sol = x_sol(:, 1:n_modes);
    U_dof = (md.Phi * q_sol.').';
end
