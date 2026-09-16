function Q = aero_generalised(cfg, G, k, U)
%AERO_GENERALISED  Generalised aerodynamic force matrix, Theodorsen strip theory.
%
%   Q = AERO_GENERALISED(cfg, G, k, U)
%
%   Returns the complex generalised aerodynamic matrix Q(k) such that, for
%   harmonic motion of the modal coordinates q at reduced frequency k, the
%   generalised aerodynamic force is Q(k)*q.
%
%   G is the struct of modal span integrals from MODAL_SPAN_INTEGRALS:
%     G.ww(i,j) = integral phi_w,i * phi_w,j   dx
%     G.wp(i,j) = integral phi_w,i * phi_phi,j dx
%     G.pw(i,j) = integral phi_phi,i * phi_w,j dx
%     G.pp(i,j) = integral phi_phi,i * phi_phi,j dx
%
%   Sectional aerodynamics
%     Theodorsen's thin-aerofoil result for a section oscillating in plunge and
%     pitch about an elastic axis at a*b aft of mid-chord, written for vertical
%     displacement w positive up and twist alpha positive nose-up:
%
%       L     = pi*rho*b^2*(-wddot + U*alphadot - b*a*alphaddot)
%             + 2*pi*rho*U*b*C(k)*(-wdot + U*alpha + b*(1/2-a)*alphadot)
%
%       M_ea  = pi*rho*b^2*(-b*a*wddot - U*b*(1/2-a)*alphadot
%                           - b^2*(1/8+a^2)*alphaddot)
%             + 2*pi*rho*U*b^2*(a+1/2)*C(k)*(-wdot + U*alpha + b*(1/2-a)*alphadot)
%
%     Substituting harmonic motion and omega = k*U/b gives the four complex
%     coefficients below, each multiplying pi*rho*U^2.
%
%   Checks built into the formulation
%     At k = 0, C = 1, the lift reduces to L = pi*rho*U^2*c*alpha, which is the
%     steady result q_inf*c*cl_alpha*alpha, and M_ea reduces to L*(x_ea - c/4),
%     the steady lift acting at the quarter chord. As k -> infinity, C -> 1/2.
%
%   The original c4p.m used, in place of all of this,
%       Q_aero(i,j) = (0.5*rho_air*b^2*(2*pi)) * (L/(i+j-1));
%   a closed-form guess with no aerodynamic derivation, a hard-coded C(k)=0.5,
%   and the structural section width substituted for the aerodynamic semi-chord.
%   Being real and symmetric positive definite, it could only ever add damping,
%   so the damping of every mode stayed negative at every speed and the script
%   returned flutter_speed = NaN.

    rho = cfg.flight.rho_air;
    b   = cfg.aero.semi_chord;
    a   = cfg.aero.a_ea;

    C = theodorsen_Ck(k);

    Lw = k^2 - 2i*k*C;
    La = 1i*k + a*k^2 + 2*C + 2i*k*C*(0.5 - a);
    Mw = a*k^2 - 2i*(a + 0.5)*k*C;
    Ma = -1i*k*(0.5 - a) + (1/8 + a^2)*k^2 + 2*(a + 0.5)*C ...
         + 2i*(a + 0.5)*(0.5 - a)*k*C;

    Q = pi*rho*U^2 * ( Lw*G.ww + b*La*G.wp + b*Mw*G.pw + b^2*Ma*G.pp );
end
