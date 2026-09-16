function G = modal_span_integrals(fem, Phi)
%MODAL_SPAN_INTEGRALS  Project the span-integral matrices onto a modal basis.
%
%   G = MODAL_SPAN_INTEGRALS(fem, Phi) returns the modal span integrals used by
%   strip-theory aerodynamics, where Phi holds the mode shapes in reduced DOFs.
%
%   See also SPAN_INTEGRAL_MATRICES, AERO_GENERALISED.

    A = span_integral_matrices(fem);
    G.ww = Phi.' * A.ww * Phi;
    G.wp = Phi.' * A.wp * Phi;
    G.pw = Phi.' * A.pw * Phi;
    G.pp = Phi.' * A.pp * Phi;
end
