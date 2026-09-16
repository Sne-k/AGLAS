function wg = gust_one_minus_cos(t, U_ds, t_g)
%GUST_ONE_MINUS_COS  Discrete "1-cosine" gust velocity time history.
%
%   wg = GUST_ONE_MINUS_COS(t, U_ds, t_g)
%
%   Implements the CS-25.341(a) / FAR 25.341(a) discrete gust shape
%
%       wg(t) = (U_ds/2) * (1 - cos(2*pi*t/t_g)),   0 <= t <= t_g
%       wg(t) = 0,                                   otherwise
%
%   over one complete cosine cycle: the velocity starts at zero, reaches the
%   peak U_ds at the midpoint t_g/2, and returns smoothly to zero at t_g.
%
%   The original code used cos(pi*t/t_g), which is only half a cycle. That
%   ramps monotonically up to U_ds and then drops discontinuously to zero at
%   t = t_g. The result is not a 1-cosine gust: it has a step in velocity, an
%   impulsive derivative, and it injects a spurious broadband transient that an
%   adaptive ODE solver will chatter on. It also doubles the impulse delivered
%   to the wing relative to the certification shape.
%
%   Inputs
%     t     time vector or scalar                    [s]
%     U_ds  design gust velocity (peak)              [m/s]
%     t_g   gust duration, one full cycle            [s]

    if t_g <= 0
        error('gust_one_minus_cos:badDuration', 't_g must be positive.');
    end
    wg = zeros(size(t));
    in = (t >= 0) & (t <= t_g);
    wg(in) = 0.5 * U_ds * (1 - cos(2*pi*t(in)/t_g));
end
