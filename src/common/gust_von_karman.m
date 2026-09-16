function [wg, t, info] = gust_von_karman(cfg, t_end, dt)
%GUST_VON_KARMAN  Continuous turbulence from the von Karman vertical PSD.
%
%   [wg, t, info] = GUST_VON_KARMAN(cfg, t_end, dt)
%
%   Synthesises a vertical gust velocity time history whose power spectral
%   density matches the von Karman model of MIL-F-8785C / MIL-HDBK-1797:
%
%                        L_w    1 + (8/3)*(1.339*L_w*Omega)^2
%     Phi_w(Omega) = s^2 ----- ---------------------------------
%                         pi   [1 + (1.339*L_w*Omega)^2]^(11/6)
%
%   with Omega the spatial frequency [rad/m], s the turbulence RMS [m/s] and
%   L_w the turbulence length scale [m]. Taylor's frozen-field hypothesis maps
%   spatial to temporal frequency, Omega = omega/U, and the temporal one-sided
%   PSD follows as S(omega) = Phi_w(omega/U)/U.
%
%   The realisation is built by random-phase Fourier synthesis: each spectral
%   line gets amplitude sqrt(2*S*domega) and a uniformly distributed phase, so
%   the sample variance converges to s^2 and the spectrum is correct by
%   construction.
%
%   Why not a Butterworth filter. The original code produced its "von Karman"
%   signal as second-order Butterworth-filtered white noise. A Butterworth
%   response rolls off in integer powers of frequency, whereas von Karman rolls
%   off as Omega^(-5/3) in the inertial subrange; no Butterworth section of any
%   order reproduces a fractional exponent. That version also had no RNG seed,
%   so it was not reproducible, and it needed butter() from the Signal
%   Processing Toolbox. This implementation is exact in the spectrum, seeded,
%   and uses only core MATLAB.
%
%   Outputs
%     wg    gust velocity time history                 [m/s]
%     t     time vector                                [s]
%     info  struct with the frequency grid, the target PSD and the achieved RMS

    if nargin < 2 || isempty(t_end), t_end = cfg.time.t_end;  end
    if nargin < 3 || isempty(dt),    dt    = cfg.time.dt_out; end

    U   = cfg.flight.U_inf;
    sig = cfg.turb.sigma;
    Lw  = cfg.turb.L_scale;

    t = 0:dt:t_end;
    N = numel(t);
    if mod(N, 2) ~= 0        % even length keeps the Hermitian bookkeeping simple
        N = N + 1;
        t = (0:N-1) * dt;
    end

    T      = N * dt;
    domega = 2*pi / T;
    kmax   = N/2 - 1;
    omega  = (1:kmax) * domega;          % strictly positive frequencies

    % Temporal one-sided PSD [ (m/s)^2 / (rad/s) ]
    Om  = omega / U;
    a   = (1.339 * Lw * Om).^2;
    Phi = sig^2 * (Lw/pi) * (1 + (8/3)*a) ./ (1 + a).^(11/6);
    S   = Phi / U;

    amp = sqrt(2 * S * domega);

    % Deterministic, reproducible phases. Seeds the global Mersenne Twister
    % stream; portable across MATLAB releases and GNU Octave.
    try
        rng(cfg.turb.seed, 'twister');
    catch
        rand('twister', cfg.turb.seed);       %#ok<RAND> legacy fallback
    end
    phase = 2*pi * rand(1, kmax);

    X = zeros(1, N);
    X(2:kmax+1) = (N/2) * amp .* exp(1i*phase);
    X(N:-1:N-kmax+1) = conj(X(2:kmax+1));      % enforce a real signal

    wg = real(ifft(X));
    wg = wg(:).';
    t  = t(:).';

    info.omega      = omega;
    info.S          = S;
    info.rms_target = sig;
    info.rms_actual = std(wg);
    info.dt         = dt;
end
