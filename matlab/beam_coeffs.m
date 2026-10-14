function [a1, a2, KT] = beam_coeffs(f, fr, p)
%BEAM_COEFFS  Timoshenko beam: Green's-function coefficients and resonator stiffness.
%
%  Eq. (4.4), with ks^2 = w^2 rhoA/GA -- the 1/ks^2 factor is REQUIRED:
%     a1 =  k1 (k2^2 - ks^2) / [ 2i GA ks^2 (k2^2 - k1^2) ]
%     a2 = -k2 (k1^2 - ks^2) / [ 2i GA ks^2 (k2^2 - k1^2) ]
%
%  Resonator dynamic stiffness, Eq. (3.1):
%     K_T = K m w^2/(m w^2 - K),   K = m (2 pi fr)^2
%
%  KT is Inf at f = fr; callers must use the limit (see beam_cmu).

  if nargin < 3 || isempty(p), p = params_beam(); end
  [k1, k2, k1sq, k2sq, ks2] = beam_wavenumbers(f, p);
  den = ks2.*(k2sq - k1sq);
  a1  =  1./(2i*p.GA) .* k1.*(k2sq - ks2)./den;
  a2  = -1./(2i*p.GA) .* k2.*(k1sq - ks2)./den;

  if nargout > 2
    w  = 2*pi*f;
    K  = p.m*(2*pi*fr)^2;
    KT = K*p.m*w.^2./(p.m*w.^2 - K);
  end
end
