function P = lyap_solve(A, Q)
%LYAP_SOLVE  Continuous Lyapunov equation A'*P + P*A = -Q, without a toolbox.
%
%   P = LYAP_SOLVE(A, Q)
%
%   Vectorising gives  (kron(I,A') + kron(A',I)) * vec(P) = -vec(Q),  a single
%   linear solve of size n^2. For the state dimensions here that is a trivial
%   cost and avoids depending on lyap from the Control System Toolbox.
%
%   A must be Hurwitz and Q symmetric positive semi-definite for P to be the
%   positive definite solution the Lyapunov argument needs.

    n = size(A, 1);
    if any(real(eig(A)) >= 0)
        error('lyap_solve:notHurwitz', ...
              'A must have all eigenvalues in the open left half plane.');
    end
    Q = (Q + Q.')/2;

    I = eye(n);
    M = kron(I, A.') + kron(A.', I);
    P = reshape(-(M \ Q(:)), n, n);
    P = (P + P.')/2;

    res = norm(A.'*P + P*A + Q, 'fro') / max(1, norm(Q, 'fro'));
    if res > 1e-8
        warning('lyap_solve:residual', ...
                'Lyapunov residual %.2e is larger than expected.', res);
    end
end
