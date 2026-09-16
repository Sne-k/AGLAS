function [rc, ro, info] = ctrb_obsv_rank(A, B, C)
%CTRB_OBSV_RANK  Controllability and observability ranks, without any toolbox.
%
%   [rc, ro, info] = CTRB_OBSV_RANK(A, B, C)
%
%   Builds the Krylov matrices [B, A*B, ..., A^(n-1)*B] and its observability
%   dual and returns their numerical ranks. ctrb and obsv live in the Control
%   System Toolbox, which this project does not depend on.
%
%   The raw Krylov matrix is badly scaled when the plant mixes structural
%   frequencies with actuator dynamics, so each block is normalised before the
%   rank test. info also returns the Hautus minimum singular values, which are
%   a better-conditioned indicator of how close a mode is to being
%   uncontrollable or unobservable than the rank alone.

    n = size(A, 1);

    Cm = zeros(n, 0);  blk = B;
    for k = 1:n
        s = norm(blk, 'fro');
        if s > 0, Cm = [Cm, blk/s]; else, Cm = [Cm, blk]; end %#ok<AGROW>
        blk = A * blk;
    end
    rc = rank(Cm, 1e-10);

    ro = NaN;  Om = [];
    if nargin >= 3 && ~isempty(C)
        Om = zeros(0, n);  blk = C;
        for k = 1:n
            s = norm(blk, 'fro');
            if s > 0, Om = [Om; blk/s]; else, Om = [Om; blk]; end %#ok<AGROW>
            blk = blk * A;
        end
        ro = rank(Om, 1e-10);
    end

    % Hautus test at each open-loop eigenvalue.
    lam = eig(A);
    info.eig = lam;
    info.hautus_ctrb = zeros(numel(lam), 1);
    info.hautus_obsv = nan(numel(lam), 1);
    for k = 1:numel(lam)
        M = [lam(k)*eye(n) - A, B];
        info.hautus_ctrb(k) = min(svd(M)) / max(1, norm(A, 'fro'));
        if nargin >= 3 && ~isempty(C)
            N = [lam(k)*eye(n) - A; C];
            info.hautus_obsv(k) = min(svd(N)) / max(1, norm(A, 'fro'));
        end
    end
    info.n = n;
end
