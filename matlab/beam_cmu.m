function [cmu, b1, b2] = beam_cmu(f, fr, p)
%BEAM_CMU  Timoshenko beam: the two Bloch branches cos(mu), Eq. (4.5).
%
%     cos^2(mu) - beta1 cos(mu) + beta2 = 0                        Eq. (4.5)
%     beta1 = (C1+C2) + i K_T (a1 S1 + a2 S2)                       Eq. (4.6)
%     beta2 = C1 C2   + i K_T (a1 S1 C2 + a2 S2 C1)                 Eq. (4.6)
%
%  beta1, beta2 are real for every real w, both below and above the shear
%  cutoff (Section 4.2), so the roots are real iff beta1^2 - 4 beta2 >= 0.
%  real() below only removes round-off imaginary parts.
%
%  VECTORIZED.  For f of length N, cmu is 2xN; entries are NaN where the roots
%  are complex (a total stop band).  Row 1 is the smaller root, row 2 the larger.
%
%  Roots are taken in Vieta form: one of them sits within ~1e-3 of +1 over
%  hundreds of Hz, and the differenced quadratic formula loses precision there.
%
%  At the resonator pole f = fr, K_T is infinite: one branch diverges and the
%  other tends to Psi of Eq. (4.9), Psi = (a1 S1 C2 + a2 S2 C1)/(a1 S1 + a2 S2).
%  Psi is returned in row 1 and NaN (the divergent branch) in row 2.

  if nargin < 3 || isempty(p), p = params_beam(); end
  f = f(:).';  N = numel(f);
  [k1, k2] = beam_wavenumbers(f, p);
  [a1, a2] = beam_coeffs(f, fr, p);
  C1 = cos(k1*p.L); C2 = cos(k2*p.L);
  S1 = sin(k1*p.L); S2 = sin(k2*p.L);

  w  = 2*pi*f;
  K  = p.m*(2*pi*fr)^2;
  KT = K*p.m*w.^2./(p.m*w.^2 - K);

  b1 = real((C1 + C2) + 1i*KT.*(a1.*S1 + a2.*S2));
  b2 = real(C1.*C2    + 1i*KT.*(a1.*S1.*C2 + a2.*S2.*C1));

  d   = b1.^2 - 4*b2;
  cmu = NaN(2, N);
  ok  = d >= 0;
  s   = sqrt(max(d, 0));
  r   = (b1 + sign0(b1).*s)/2;
  lo  = r; hi = b2./r;
  z   = (r == 0);
  lo(z) = -s(z)/2; hi(z) = s(z)/2;
  cmu(1, ok) = min(lo(ok), hi(ok));
  cmu(2, ok) = max(lo(ok), hi(ok));

  % ---- pole: replace by the analytic limit, Eq. (4.9) ---------------------
  at = abs(f - fr) <= p.tol_pole*max(fr, 1);
  if any(at)
    num = a1(at).*S1(at).*C2(at) + a2(at).*S2(at).*C1(at);
    den = a1(at).*S1(at) + a2(at).*S2(at);
    cmu(1, at) = real(num./den);
    cmu(2, at) = NaN;
    b1(at) = Inf; b2(at) = Inf;
  end
end

function s = sign0(x)
  s = sign(x); s(s == 0) = 1;
end
