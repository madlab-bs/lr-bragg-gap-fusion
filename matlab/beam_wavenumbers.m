function [k1, k2, k1sq, k2sq, ks2] = beam_wavenumbers(f, p)
%BEAM_WAVENUMBERS  Timoshenko host wavenumbers, roots of Eq. (4.2) (Section 4.1).
%
%   k^4 - A w^2 k^2 + (B w^4 - C w^2) = 0,   A,B,C from params_beam,
%   solved as a quadratic in k^2 by the closed-form root formula.
%
%  Branch rule (e^{+i w t}, outgoing/decaying field e^{-i k |x|}, Im k <= 0 as in
%  Section 4.1):
%     k1 = +sqrt(k1sq)                          (always propagating)
%     k2 = -1i*sqrt(-k2sq)   below the shear cutoff   (k2 = -i kappa, decaying)
%     k2 = +sqrt(k2sq)       above the shear cutoff   (outgoing, by the same
%                                                      convention as k1)
%
%  Below w_c one root has k^2 > 0 and the other k^2 < 0; above w_c both are
%  positive, Eq. (4.3).  The discriminant is strictly positive (Section 4.1), so
%  k1sq and k2sq are always real and distinct.

  if nargin < 2 || isempty(p), p = params_beam(); end
  w  = 2*pi*f;
  disc = p.A^2*w.^4 - 4*(p.B*w.^4 - p.C*w.^2);   % = (rhoI/EI - rhoA/GA)^2 w^4 + 4 rhoA/EI w^2
  sq   = sqrt(max(disc, 0));
  k1sq = (p.A*w.^2 + sq)/2;
  k2sq = (p.A*w.^2 - sq)/2;
  ks2  = w.^2*p.rhoA/p.GA;                       % k_s^2 of Eq. (4.4)

  k1 = sqrt(k1sq);
  k2 = zeros(size(k2sq));
  neg = k2sq < 0;
  k2(neg)  = -1i*sqrt(-k2sq(neg));
  k2(~neg) = sqrt(k2sq(~neg));
end
