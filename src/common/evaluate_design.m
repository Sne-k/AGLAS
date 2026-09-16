function out = evaluate_design(cfg, t_end)
%EVALUATE_DESIGN  Full re-analysis of one candidate wing design.
%
%   out = EVALUATE_DESIGN(cfg)
%   out = EVALUATE_DESIGN(cfg, t_end)
%
%   Runs the whole chain for the configuration given: section properties, FE
%   assembly, modal analysis, gust response and stress recovery, and returns
%   the quantities a sizing loop needs.
%
%   Every call recomputes the section from geometry rather than scaling a
%   previous answer. That matters because the factor of safety is not a simple
%   power of wall thickness. For a thin-walled box almost all the bending
%   material sits at a fixed distance from the neutral axis, so I grows close
%   to linearly with wall thickness and FoS does too; for a solid rectangle
%   I goes as h^3 while the outer fibre distance goes as h, so FoS goes as h^2.
%   Changing the thickness also changes the mass, hence the frequencies, hence
%   the dynamic amplification. Assuming any fixed scaling law, as the original
%   optimiser did with FoS_scaled = FoS * thickness_factor, is wrong for one
%   of those two sections and only approximately right for the other.
%
%   Output fields
%     min_FoS      worst factor of safety over span and time
%     peak_stress  peak outer-fibre stress                          [Pa]
%     mass         structural mass of the semi-span                 [kg]
%     area         cross-sectional area of material                 [m^2]
%     f1           fundamental frequency                            [Hz]
%     tip_defl     peak tip deflection                              [m]
%     x, FoS_span  spanwise station and worst-in-time FoS

    if nargin < 2 || isempty(t_end)
        t_end = max(4*cfg.gust.t_g, 4);
    end

    cfg.section = section_properties(cfg);

    fem   = beam_fem(cfg);
    modes = modal_analysis(fem, cfg.fem.n_modes);

    t_grid = 0:cfg.time.dt_out:t_end;
    wg     = gust_one_minus_cos(t_grid, cfg.gust.U_ds, cfg.gust.t_g);

    [~, ~, U_dof] = solve_modal_response(cfg, fem, modes, t_grid, wg, ...
                                         'kussner', cfg.gust.t_g/20);

    rec = recover_bending(cfg, fem, U_dof.');

    FoS_field = cfg.material.sigma_yield ./ max(abs(rec.stress), eps);

    out.min_FoS     = min(FoS_field(:));
    out.FoS_span    = min(FoS_field, [], 2);
    out.peak_stress = max(abs(rec.stress(:)));
    out.area        = cfg.section.A;
    out.mass        = cfg.section.m_struct * cfg.geom.L;
    out.mass_total  = cfg.section.m_total  * cfg.geom.L;
    out.f1          = modes.freq_hz(1);
    out.tip_defl    = max(abs(U_dof(:, fem.idx_w(end))));
    out.x           = rec.x(:);
    out.wall_t      = cfg.geom.wall_t;
end
