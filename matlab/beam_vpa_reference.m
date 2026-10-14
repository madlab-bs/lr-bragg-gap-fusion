function V = beam_vpa_reference(R, p, D)
%BEAM_VPA_REFERENCE  Independent D-digit references for the validators.
%
%  V = beam_vpa_reference(R, p, D) recomputes, in vpa with the sandwich inputs (Table 1(a) of
%  Sharma & Sun, 2016) held as exact rationals, the quantities that beam_reference_values (R)
%  computes in double:
%     V.fB1     closed-form root of k1 L = pi (quadratic in w^2, flexural root, Vieta form)
%     V.frstar  closed form of beta1(f_B1) = 0, equivalently Eq. (5.2)
%     V.hi      upper edge of the coincidence band (cos mu = +1), so V.W = V.hi - V.fB1
%     V.lo_star lower edge of the mode-exchange band (branch degeneracy), V.Wstar = V.fB1 - V.lo_star
%     V.psiB    the pole limit Psi, Eq. (4.9), at f = f_r = f_B1 (expected -1, Eq. (4.10))
%  The double values in R are used only to bracket each edge (+-1e-7 Hz); the bracket must contain
%  a sign change of the vpa edge function, which is then bisected in vpa to 10^(8-D) Hz.  Nothing
%  here calls the double-precision dispersion routines.  Returned values are double.
%  Requires the Symbolic Math Toolbox.  Same algebra as verification/headline_vpa.m.

  if nargin < 2 || isempty(p), p = params_beam(); end
  if nargin < 3 || isempty(D), D = 40; end
  d0 = digits;  digits(D);  restore = onCleanup(@() digits(d0));
  P = struct();
  for fn = {'Ef','hf','rhof','rhoc','hc','Gc','b','L','m'}, P.(fn{1}) = sym(p.(fn{1})); end
  P.EI   = P.Ef*P.b*(P.hc^2*P.hf/2 + P.hc*P.hf^2);
  P.GA   = P.Gc*P.b*P.hc;
  P.rhoA = 2*P.rhof*P.b*P.hf + P.rhoc*P.b*P.hc;
  P.rhoI = P.rhof*P.b*(P.hc^2*P.hf/2 + P.hc*P.hf^2) + P.rhoc*P.b*P.hc^3/12;
  P.A = P.rhoI/P.EI + P.rhoA/P.GA;  P.B = P.rhoI*P.rhoA/(P.EI*P.GA);  P.C = P.rhoA/P.EI;

  k = sym(pi)/P.L;  bb = P.A*k^2 + P.C;
  fB = sqrt(vpa(2*k^4/(bb + sqrt(bb^2 - 4*P.B*k^4))))/(2*sym(pi));

  [~, a2, C2, S2] = coeffs(fB, P);
  KTs = 2*sin(waven2(fB, P)*P.L/2)^2/(1i*a2*S2);           % (1 - C2)/(i a2 S2)
  wB2 = (2*sym(pi)*fB)^2;
  Ks  = real(KTs)*P.m*wB2/(P.m*wB2 + real(KTs));
  frs = sqrt(Ks/P.m)/(2*sym(pi));

  hi = bisect(@(f) edge(f, fB, P, 'p'), R.coinc(2), D);
  lo = bisect(@(f) edge(f, frs, P, 'd'), R.mexch(1), D);
  [a1, a2, C2, S2, C1, S1] = coeffs(fB, P);
  psi = real((a1*S1*C2 + a2*S2*C1)/(a1*S1 + a2*S2));

  V = struct('D', D, 'fB1', double(fB), 'frstar', double(frs), 'hi', double(hi), ...
             'W', double(hi - fB), 'lo_star', double(lo), 'Wstar', double(fB - lo), ...
             'psiB', double(psi), 'psiB_err', double(abs(psi + 1)));
end

function r = bisect(fun, e0, D)
  a = vpa(sym(e0 - 1e-7, 'f'));  b = vpa(sym(e0 + 1e-7, 'f'));
  fa = sign(double(fun(a)));
  assert(fa ~= sign(double(fun(b))), 'beam_vpa_reference: no sign change within 1e-7 Hz of %.12f', e0);
  while (b - a) > sym(10)^(8 - D)
    m = (a + b)/2;
    if sign(double(fun(m))) == fa, a = m; else, b = m; end
  end
  r = (a + b)/2;
end

function [k1, k2, k1sq, k2sq, ks2] = waven(f, P)
  w = 2*sym(pi)*vpa(f);
  sq = sqrt(P.A^2*w^4 - 4*(P.B*w^4 - P.C*w^2));
  k1sq = (P.A*w^2 + sq)/2;  k2sq = (P.A*w^2 - sq)/2;  ks2 = w^2*P.rhoA/P.GA;
  k1 = kd(k1sq);  k2 = kd(k2sq);
end
function k2 = waven2(f, P), [~, k2] = waven(f, P); end
function k = kd(ksq), k = sqrt(ksq); if double(imag(k)) > 0, k = -k; end, end

function [a1, a2, C2, S2, C1, S1] = coeffs(f, P)
  [k1, k2, k1sq, k2sq, ks2] = waven(f, P);
  den = ks2*(k2sq - k1sq);
  a1 =  1/(2i*P.GA) * k1*(k2sq - ks2)/den;
  a2 = -1/(2i*P.GA) * k2*(k1sq - ks2)/den;
  C1 = cos(k1*P.L);  C2 = cos(k2*P.L);  S1 = sin(k1*P.L);  S2 = sin(k2*P.L);
end

function v = edge(f, fr, P, which)
  [a1, a2, C2, S2, C1, S1] = coeffs(f, P);
  w2 = (2*sym(pi)*vpa(f))^2;  K = P.m*(2*sym(pi)*fr)^2;  KT = K*P.m*w2/(P.m*w2 - K);
  b1 = real((C1 + C2) + 1i*KT*(a1*S1 + a2*S2));
  b2 = real(C1*C2 + 1i*KT*(a1*S1*C2 + a2*S2*C1));
  switch which
    case 'p', v = 1 - b1 + b2;
    otherwise, v = b1^2 - 4*b2;
  end
end
