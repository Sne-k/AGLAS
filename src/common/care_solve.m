function [P, K, info] = care_solve(A, B, Q, R)
%CARE_SOLVE  Continuous algebraic Riccati equation, without any toolbox.
%
%   [P, K, info] = CARE_SOLVE(A, B, Q, R)
%
%   Solves      A'*P + P*A - P*B*inv(R)*B'*P + Q = 0
%   and returns the optimal gain  K = inv(R)*B'*P,  so that u = -K*x minimises
%   the cost  integral( x'*Q*x + u'*R*u ) dt.
%
%   Method: the Hamiltonian matrix
%
%       H = [  A   -B*inv(R)*B'
%             -Q   -A'          ]
%
%   has eigenvalues symmetric about the imaginary axis. Taking the n
%   eigenvectors with negative real part and partitioning them as [X1; X2]
%   gives the stabilising solution P = X2 / X1. This is the standard
%   eigenvector method and needs only eig and a linear solve.
%
%   Written out rather than calling lqr or care because the Control System
%   Toolbox is not installed everywhere, and the rest of this project has no
%   toolbox dependencies. Keeping that property means the control design runs
%   on a bare MATLAB installation and under GNU Octave.
%
%   info returns the Riccati residual and the closed-loop eigenvalues, so the
%   caller can confirm the solution rather than trusting it.

    n = size(A, 1);
    if size(Q,1) ~= n || size(Q,2) ~= n
        error('care_solve:badQ', 'Q must be %dx%d.', n, n);
    end
    if rcond(R) < eps
        error('care_solve:singularR', 'R must be positive definite.');
    end

    G = B * (R \ B.');
    G = (G + G.')/2;
    Q = (Q + Q.')/2;

    H = [ A,  -G;
         -Q,  -A.' ];

    [V, D] = eig(H);
    lam = diag(D);

    stable = real(lam) < 0;
    if sum(stable) ~= n
        error('care_solve:noSplit', ...
            ['Found %d stable Hamiltonian eigenvalues, expected %d. The pair ' ...
             '(A,B) is probably not stabilisable, or Q is not positive ' ...
             'semi-definite on the unobservable subspace.'], sum(stable), n);
    end

    Vs = V(:, stable);
    X1 = Vs(1:n,     :);
    X2 = Vs(n+1:end, :);

    if rcond(X1) < 1e-14
        error('care_solve:illConditioned', ...
              'The stable invariant subspace is ill conditioned; rescale Q and R.');
    end

    P = real(X2 / X1);
    P = (P + P.')/2;

    K = R \ (B.' * P);

    res = A.'*P + P*A - P*G*P + Q;
    info.residual   = norm(res, 'fro') / max(1, norm(Q, 'fro'));
    info.closed_loop_eig = eig(A - B*K);
    info.stable     = all(real(info.closed_loop_eig) < 0);
    info.P_posdef   = all(eig(P) > -1e-9*max(1, norm(P)));
end
