function u_ad = aglas_u_adaptive(in)
%AGLAS_U_ADAPTIVE  Adaptive control increment, u_ad = theta' * phi.
%
%   u_ad = AGLAS_U_ADAPTIVE([theta; xhat])
%
%   Kept as its own block so the adaptive contribution can be logged and
%   inspected separately from the baseline LQG command. When diagnosing an
%   adaptive loop, the first question is always how much of the command the
%   adaptation is responsible for, and that is impossible to answer if the two
%   are summed inside one block.

    S  = evalin('base', 'AGLAS');
    ns = S.n_states;

    theta = in(1:ns);
    phi   = in(ns+1 : ns+ns);

    u_ad = theta.' * phi;
end
