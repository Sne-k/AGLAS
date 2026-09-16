function theta_dot = aglas_adapt_law(in)
%AGLAS_ADAPT_LAW  MRAC parameter update, for an Interpreted MATLAB Function block.
%
%   theta_dot = AGLAS_ADAPT_LAW([xhat; xm; theta])
%
%   Implements
%       theta_dot = -Gamma*( phi*(PB'*e)/(1 + phi'*phi) + sigma*theta )
%   with phi the estimated structural states and e = phi - xm the deviation of
%   the aircraft from its reference model.
%
%   A soft barrier replaces the hard projection used in the .m simulation: once
%   the parameter norm passes theta_max, the outward component of theta_dot is
%   removed. A hard reset would introduce a discontinuity that a variable-step
%   solver has to chase, and Simulink would take very small steps around it.

    S = evalin('base', 'AGLAS');

    ns = S.n_states;
    na = S.n_aug;

    xhat  = in(1:na);
    xm    = in(na+1 : na+ns);
    theta = in(na+ns+1 : na+ns+ns);

    phi = xhat(1:ns);
    e   = phi - xm;

    s_e = S.PB.' * e;
    if S.normalise
        s_e = s_e / (1 + phi.'*phi);
    end

    theta_dot = -S.Gamma * (phi*s_e + S.sigma*theta);

    nrm = norm(theta);
    if nrm > S.theta_max
        outward = (theta.'*theta_dot) / max(nrm^2, eps);
        if outward > 0
            theta_dot = theta_dot - outward*theta;
        end
    end
end
