function C = theodorsen_Ck(k)
%THEODORSEN_CK  Theodorsen's circulation function C(k).
%
%   C = THEODORSEN_CK(k) evaluates
%
%       C(k) = H1^(2)(k) / ( H1^(2)(k) + 1i * H0^(2)(k) )
%
%   where Hn^(2) = Jn - 1i*Yn is the Hankel function of the second kind and k
%   is the reduced frequency k = omega*b/U with b the semi-chord.
%
%   Limits:  C(0) = 1 (quasi-steady),  C(k) -> 0.5 as k -> infinity.
%
%   Uses only besselj and bessely, which are core MATLAB; no toolbox needed.
%
%   The original flutter script hard-coded C_k = 0.5 and labelled it "C(k) for
%   k = 0.2 (real part)". The true value at k = 0.2 is 0.7276 - 0.1886i, so the
%   magnitude was wrong by 31 % and the phase lag, which is the entire physical
%   mechanism by which unsteady aerodynamics feeds energy into or out of the
%   structure, was discarded.

    k = double(k);
    C = zeros(size(k));

    zero_k = (k == 0);
    C(zero_k) = 1;

    kk = k(~zero_k);
    if ~isempty(kk)
        H0 = besselj(0, kk) - 1i*bessely(0, kk);
        H1 = besselj(1, kk) - 1i*bessely(1, kk);
        C(~zero_k) = H1 ./ (H1 + 1i*H0);
    end
end
