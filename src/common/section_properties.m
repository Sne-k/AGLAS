function sec = section_properties(cfg)
%SECTION_PROPERTIES  Cross-section stiffness and inertia properties.
%
%   sec = SECTION_PROPERTIES(cfg) returns area, second moment of area, torsion
%   constant and the derived rigidities and inertias for the section described
%   by cfg.geom.
%
%   Supported sections
%     'solid' : solid rectangle, width box_w x height box_h.
%     'box'   : thin-walled closed rectangular box, uniform wall thickness.
%               Torsion constant from the Bredt-Batho single-cell formula
%                   J = 4 * Am^2 / (contour integral of ds/t)
%               which for a constant wall thickness reduces to
%                   J = 4 * Am^2 * t / s_m
%               with Am the area enclosed by the wall mid-line and s_m its
%               perimeter.
%
%   Output fields
%     A        cross-sectional area of load-carrying material      [m^2]
%     I        second moment of area about the chordwise axis      [m^4]
%     J        Saint-Venant torsion constant                       [m^4]
%     c_outer  distance neutral axis -> extreme fibre              [m]
%     EI, GJ   bending and torsional rigidities                    [N m^2]
%     m_struct structural mass per unit span                       [kg/m]
%     m_total  structural + non-structural mass per unit span      [kg/m]
%     I_theta  polar mass moment of inertia per unit span about the
%              elastic axis                                        [kg m^2/m]

    g = cfg.geom;
    m = cfg.material;

    switch lower(g.section)
        case 'solid'
            sec.A = g.box_w * g.box_h;
            sec.I = g.box_w * g.box_h^3 / 12;
            % Saint-Venant torsion constant of a solid rectangle,
            % J = beta * a * b^3 with a >= b (Roark / Timoshenko series).
            a = max(g.box_w, g.box_h);
            b = min(g.box_w, g.box_h);
            r = b / a;
            beta = 1/3 - 0.21*r*(1 - r^4/12);
            sec.J = beta * a * b^3;

        case 'box'
            t = g.wall_t;
            if ~isfinite(t) || t <= 0
                error('section_properties:wallThickness', ...
                      'Box section requires a positive wall thickness.');
            end
            if 2*t >= min(g.box_w, g.box_h)
                error('section_properties:wallTooThick', ...
                      'Wall thickness %.4f m is too large for a %.3f x %.3f m box.', ...
                      t, g.box_w, g.box_h);
            end
            wi = g.box_w - 2*t;   hi = g.box_h - 2*t;       % inner dimensions
            sec.A = g.box_w*g.box_h - wi*hi;
            sec.I = (g.box_w*g.box_h^3 - wi*hi^3) / 12;

            Am  = (g.box_w - t) * (g.box_h - t);            % wall mid-line area
            s_m = 2 * ((g.box_w - t) + (g.box_h - t));      % wall mid-line perimeter
            sec.J = 4 * Am^2 * t / s_m;

        otherwise
            error('section_properties:unknownSection', ...
                  'Unknown section type "%s".', g.section);
    end

    sec.c_outer  = g.box_h / 2;
    sec.EI       = m.E * sec.I;
    sec.GJ       = m.G * sec.J;
    sec.m_struct = m.rho * sec.A;
    sec.m_total  = sec.m_struct + g.m_nonstruct;

    % Polar mass moment of inertia per unit span about the elastic axis.
    % The section mass is idealised as uniformly smeared over the chord, then
    % shifted to the elastic axis by the parallel-axis theorem.
    d_ea_cg  = (g.cg_frac - g.ea_frac) * g.chord;
    sec.I_theta = sec.m_total * (g.chord^2 / 12 + d_ea_cg^2);
    sec.d_ea_cg = d_ea_cg;
end
