function r = probe_signal(cfg, t, wg, d, w_ind, C)
%PROBE_SIGNAL  What a gust sensor d metres ahead of the wing measures.
%
%   r = PROBE_SIGNAL(cfg, t, wg, d)
%   r = PROBE_SIGNAL(cfg, t, wg, d, w_ind, C)
%
%   wg is the gust that reaches the wing on the uniform time grid t. Under
%   Taylor's frozen-turbulence hypothesis the sensor sees the same gust
%   tau = d/U earlier, so with C = 0 (or no w_ind) the result is simply wg
%   advanced by tau.
%
%   Real turbulence evolves while it travels, so the gust that arrives is not
%   exactly the gust that was measured, and the loss grows with frequency and
%   with distance. With C > 0 the measurement is mixed, frequency by frequency,
%   with an independent turbulence realisation w_ind of the same spectrum:
%
%       R(f) = gamma(f) * W_ahead(f) + sqrt(1 - gamma(f)^2) * W_ind(f)
%       gamma(f) = exp(-C * |f| * d / U)
%
%   an exponential (Davenport-type) coherence. C is a decay constant; it is
%   swept rather than asserted, since its value for streamwise evolution over
%   tens of metres is not well established. The model ignores the lidar's
%   own spatial averaging along its beam and any spanwise gust variation.
%
%   Pass the result as ctrl.ff.ref to CLOSED_LOOP_SIM, with ff.probe_ahead = d
%   and enough FIR taps to span the preview.

    t  = t(:);  wg = wg(:);
    U  = cfg.flight.U_inf;
    r  = interp1(t, wg, t + d/U, 'linear', wg(end));

    if nargin < 6 || C == 0 || d == 0
        return;
    end

    n  = numel(t);
    dt = t(2) - t(1);
    f  = (0:n-1).' / (n*dt);
    f(f > 1/(2*dt)) = f(f > 1/(2*dt)) - 1/dt;        % two-sided frequencies
    g  = exp(-C * abs(f) * d / U);
    R  = g .* fft(r) + sqrt(1 - g.^2) .* fft(w_ind(:));
    r  = real(ifft(R));
end
