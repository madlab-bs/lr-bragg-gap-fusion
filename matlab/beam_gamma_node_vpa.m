function [g, fB1str] = beam_gamma_node_vpa(eta_r, p, D)
%BEAM_GAMMA_NODE_VPA  gamma(f_B1) under pure resonator loss, f_r = f_B1, in D-digit vpa.
%
%  Requires the Symbolic Math Toolbox.
%
%  The damped chain beam_wavenumbers_d -> beam_coeffs_d -> beam_cmu_d quadratic, ported to
%  vpa with NO nodal snap and NO pole substitution (eta_r > 0, so K_T is finite).  The
%  sandwich inputs of params_beam (Table 1(a) of Sharma & Sun, 2016) are held as exact
%  rationals.  f_B1 is the closed-form root of k1 L = pi: with k = pi/L the characteristic
%  equation (4.2), k^4 - A w^2 k^2 + B w^4 - C w^2 = 0, is a quadratic in w^2, whose
%  flexural (lower) root is taken in Vieta form.  The same port, its convergence in D, and
%  four controls are in verification/gamma_fB1_vpa.m (output gamma_fB1_vpa_output.txt).
%
%  The returned values are the D-digit round-off floor of an identity whose exact value
%  is 0 (Section 6.3); they scale as 10^(-D/2) and are upper bounds, returned as double.

  if nargin < 2 || isempty(p), p = params_beam(); end
  if nargin < 3 || isempty(D), D = 100; end
  assert(all(eta_r > 0), 'eta_r must be positive: the K_T pole is not handled here');
  d0 = digits;  digits(D);  restore = onCleanup(@() digits(d0));

  P = struct();
  for fn = {'Ef','hf','rhof','rhoc','hc','Gc','b','L','m'}
    P.(fn{1}) = sym(p.(fn{1}));                  % 'r' flag: 1.2e-3 -> 3/2500
  end
  EI   = P.Ef*P.b*(P.hc^2*P.hf/2 + P.hc*P.hf^2);
  GA   = P.Gc*P.b*P.hc;
  rhoA = 2*P.rhof*P.b*P.hf + P.rhoc*P.b*P.hc;
  rhoI = P.rhof*P.b*(P.hc^2*P.hf/2 + P.hc*P.hf^2) + P.rhoc*P.b*P.hc^3/12;
  A = rhoI/EI + rhoA/GA;  B = rhoI*rhoA/(EI*GA);  C = rhoA/EI;

  % f_B1 from k1 L = pi (closed form)
  kB = sym(pi)/P.L;  bb = A*kB^2 + C;
  fB = sqrt(vpa(2*kB^4/(bb + sqrt(bb^2 - 4*B*kB^4))))/(2*sym(pi));
  fB1str = char(vpa(fB, 30));

  % lossless host: wavenumbers and Green's-function coefficients at f_B1 (eta_b = eta_s = 0)
  w = 2*sym(pi)*fB;
  sq = sqrt(A^2*w^4 - 4*(B*w^4 - C*w^2));
  k1sq = (A*w^2 + sq)/2;  k2sq = (A*w^2 - sq)/2;  ks2 = w^2*rhoA/GA;
  k1 = kdecay(k1sq);  k2 = kdecay(k2sq);
  den = ks2*(k2sq - k1sq);
  a1 =  1/(2i*GA) * k1*(k2sq - ks2)/den;
  a2 = -1/(2i*GA) * k2*(k1sq - ks2)/den;
  C1 = cos(k1*P.L);  C2 = cos(k2*P.L);  S1 = sin(k1*P.L);  S2 = sin(k2*P.L);

  g = zeros(size(eta_r));
  for k = 1:numel(eta_r)
    Kc = P.m*w^2*(1 + 1i*sym(eta_r(k), 'f'));    % f_r = f_B1
    KT = Kc*P.m*w^2/(P.m*w^2 - Kc);
    b1 = (C1 + C2) + 1i*KT*(a1*S1 + a2*S2);
    b2 = C1*C2     + 1i*KT*(a1*S1*C2 + a2*S2*C1);
    sqd = sqrt(b1^2 - 4*b2);
    sgn = 1;  if double(real(conj(b1)*sqd)) < 0, sgn = -1; end
    rbig = (b1 + sgn*sqd)/2;  rsmall = b2/rbig;
    g(k) = double(min(abs(imag(acos(rsmall))), abs(imag(acos(rbig)))));
  end
end

function k = kdecay(ksq)
  k = sqrt(ksq);  if double(imag(k)) > 0, k = -k; end
end
