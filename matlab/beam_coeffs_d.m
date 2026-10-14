function [a1, a2, KT] = beam_coeffs_d(f, fr, eta, p)
%BEAM_COEFFS_D  Damped Green's-function coefficients and dynamic stiffness.
%
%  Same closed form as beam_coeffs.m, Eqs. (4.4) and (3.1), with the complex-modulus
%  substitutions of Section 6.3: complex GA and a complex resonator spring
%  K -> K(1+i eta_r):
%     a1 =  k1 (k2^2 - ks^2) / [ 2i GA(1+i eta_s) ks^2 (k2^2 - k1^2) ]
%     a2 = -k2 (k1^2 - ks^2) / [ 2i GA(1+i eta_s) ks^2 (k2^2 - k1^2) ]
%     K_T = Kc m w^2 / (m w^2 - Kc),   Kc = m(2 pi fr)^2 (1 + i eta_r)
%
%  eta = [eta_b, eta_s, eta_r].  With eta_r > 0 the K_T pole leaves the real axis
%  (Kc/m = wr^2(1+i eta_r)), so K_T is finite for every real f: damping
%  regularizes the resonance.  KT is Inf only at eta_r = 0, f = fr.

  if nargin < 4 || isempty(p),   p   = params_beam(); end
  if nargin < 3 || isempty(eta), eta = [0 0 0];       end

  GAc = p.GA*(1 + 1i*eta(2));
  [k1, k2, k1sq, k2sq, ks2] = beam_wavenumbers_d(f, eta, p);
  den = ks2.*(k2sq - k1sq);
  a1  =  1./(2i*GAc) .* k1.*(k2sq - ks2)./den;
  a2  = -1./(2i*GAc) .* k2.*(k1sq - ks2)./den;

  if nargout > 2
    w  = 2*pi*f;
    Kc = p.m*(2*pi*fr)^2*(1 + 1i*eta(3));
    KT = Kc*p.m*w.^2./(p.m*w.^2 - Kc);
  end
end
