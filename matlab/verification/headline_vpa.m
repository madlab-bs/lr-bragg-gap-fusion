function R = headline_vpa(D)
%HEADLINE_VPA  Headline values recomputed in D-digit vpa, independently of the
%  double-precision production routines (no gap_edges, narrow_band, bisection tolerances,
%  snaps or pole windows).  Writes only headline_vpa_output_<D>.txt and .mat here.
%  Run from matlab/verification/:  headline_vpa(50)   (Symbolic Math Toolbox; about two
%  minutes).  Run it before reference_values_check, which reads headline_vpa_results_50.mat.
%
%  Lossless Timoshenko beam, exact rational inputs from Table 1(a) of ref. [44].
%    cos^2 mu - b1 cos mu + b2 = 0,  b1 = C1+C2 + i KT (a1 S1 + a2 S2),
%                                    b2 = C1 C2 + i KT (a1 S1 C2 + a2 S2 C1)    Eqs. (4.5)-(4.6)
%  Band-edge functions:  Em(f) = 1 + b1 + b2   (a branch at cos mu = -1, mu = pi)
%                        Ep(f) = 1 - b1 + b2   (a branch at cos mu = +1, mu = 0)
%                        Dl(f) = b1^2 - 4 b2   (branch degeneracy, Eq. (5.6))
%  Unit-cell config (b) (sliding, mu = 0): 1 - i KT (a1 cot h1 + a2 cot h2) = 0     Eq. (5.1)(b)
%  Each root is bracketed from a coarse double scan of the SAME vpa functions, then bisected in
%  vpa to below 10^(-D+8) Hz, and the bracketing function's identity is reported.

  if nargin < 1 || isempty(D), D = 50; end
  here = fileparts(mfilename('fullpath'));
  addpath(fullfile(here, '..', '..', 'matlab'));
  d0 = digits;  digits(D);  restore = onCleanup(@() digits(d0));
  out = fullfile(here, sprintf('headline_vpa_output_%d.txt', D));
  fid = fopen(out, 'w');  cl = onCleanup(@() fclose(fid));
  pr = @(varargin) both(fid, varargin{:});
  pb = params_beam();  P = symparams(pb);
  pr('headline_vpa  D = %d digits  %s\n\n', D, char(datetime('now')));

  % ---- Bragg frequencies, shear cutoff, evanescent floor -------------------------------
  fB1 = fBragg(P, 1);  fB2 = fBragg(P, 2);
  fc  = sqrt(P.GA/P.rhoI)/(2*sym(pi));
  [k1, k2] = waven(fB1, P);
  kL  = abs(imag(k2))*P.L;
  pr('f_B1  = %s Hz   (|sin k1L| = %.2g)\n', s(fB1), dbl(abs(sin(k1*P.L))));
  pr('f_B2  = %s Hz   (|sin k1L| = %.2g)\n', s(fB2), dbl(abs(sin(waven(fB2, P)*P.L))));
  pr('f_c   = %s Hz\n', s(fc));
  pr('kappaL= %s   (acosh(cos k2L) check: %.2g)\n', s(kL), dbl(abs(acosh(real(cos(k2*P.L))) - kL)));

  % ---- residual band at exact coincidence, f_r = f_B1 ----------------------------------
  pr('\n-- f_r = f_B1 --\n');
  Em0 = dbl(psiB(fB1, fB1, P));
  pr('f = f_r = f_B1 is the K_T pole; finite branch Psi, Eq. (4.9), = %s (expect -1)\n', s(psiB(fB1, fB1, P)));
  [hi, id] = edge_root(fB1, P, fB1 + 0.005, fB1 + 0.08, D);
  W = hi - fB1;
  pr('upper edge  = %s Hz  [%s]\n', s(hi), id);
  pr('W           = %s Hz\n', s(W));

  % ---- mode exchange f_r*: beta1(f_B1) = 0  <=>  Eq. (5.2), 1 = i K_T a2 cot h2 -------
  pr('\n-- mode exchange --\n');
  [~, ~, a2, C2, S2] = coeffs(fB1, P);
  KTs = (1 - C2)/(1i*a2*S2);
  wB2 = (2*sym(pi)*fB1)^2;
  Ks  = KTs*P.m*wB2/(P.m*wB2 + KTs);             % K_T = K m w^2/(m w^2 - K) solved for K
  frs = sqrt(real(Ks)/P.m)/(2*sym(pi));
  pr('f_r*        = %s Hz   (Im K* / Re K* = %.2g)\n', s(frs), dbl(imag(Ks)/real(Ks)));
  h2 = waven2(fB1, P)*P.L/2;
  pr('Eq. (5.2) residual 1 - i K_T a2 cot h2 at f_r*: %.2g\n', ...
     dbl(abs(1 - 1i*KTof(fB1, frs, P)*a2*cot(h2))));
  pr('second root at f_B1, f_r*: C2 + i K_T a2 S2 = %s (expect +1)\n', ...
     s(real(C2 + 1i*KTof(fB1, frs, P)*a2*S2)));
  % double gap_edges only BRACKETS the pass band containing f_B1 - 0.01 Hz; vpa identifies the
  % vanishing edge function in each bracket and bisects it.  (A first-sign-change rule is wrong
  % here: a branch degeneracy at ~636.4082 Hz lies inside the stop band, where both roots > 1.)
  [lo, hiS, idl, idh, dlo, dhi] = band_containing(frs, fB1 - sym(1)/100, pb, P, D);
  Ws = hiS - lo;
  pr('pass band   = [%s, %s] Hz  [%s | %s]\n', s(lo), s(hiS), idl, idh);
  pr('upper edge - f_B1 = %.3g Hz;  vpa minus double edges: %.2g, %.2g Hz\n', dbl(hiS - fB1), dlo, dhi);
  pr('W*          = %s Hz\n', s(Ws));
  pr('W(f_r) of the same band near f_r* (is f_r* the minimizer?):\n');
  for dl = [-1e-2 -1e-3 -1e-4 1e-4 1e-3 1e-2]
    [l2, h2b] = band_containing(frs + sym(dl, 'f'), fB1 - sym(1)/100, pb, P, D);
    pr('  f_r* %+7.0e Hz : W - W* = %+.6e Hz\n', dl, dbl(h2b - l2 - Ws));
  end

  % ---- degeneracy edge and sliding-cell mode at f_r = 200 Hz ---------------------------
  pr('\n-- f_r = 200 Hz --\n');
  fr = sym(200);
  fdeg = root_of(@(f) disc(f, fr, P), sym(364.2559), sym(364.25605), D);
  pr('degeneracy  = %s Hz\n', s(fdeg));
  fsl  = root_of(@(f) cfgb(f, fr, P), sym(364.2565), sym(364.2573), D);
  pr('config (b)  = %s Hz\n', s(fsl));
  pr('offset      = %s Hz\n', s(fsl - fdeg, 6));

  % ---- rod (exact) --------------------------------------------------------------------
  prd = params_rod();
  th = sym(prd.m)/(sym(prd.rhoA)*sym(prd.L));
  G1 = sym(pi)^2*th/4;
  pr('\n-- rod --\ntheta_c = %s = %s\nGamma_1 = %s\nc = %s m/s, f_B1 = %s Hz\n', char(th), s(th), s(G1), ...
     s(sqrt(sym(prd.EA)/sym(prd.rhoA))), s(sqrt(sym(prd.EA)/sym(prd.rhoA))/(2*sym(prd.L))));

  R = struct('D', D, 'fB1', s(fB1), 'fB2', s(fB2), 'fc', s(fc), 'kappaL', s(kL), 'W', s(W), ...
    'hi', s(hi), 'frstar', s(frs), 'Wstar', s(Ws), 'lo_star', s(lo), 'hi_star', s(hiS), 'fdeg', s(fdeg), ...
    'fslide', s(fsl), 'theta_c', s(th), 'Gamma1', s(G1), 'Em0', Em0);
  save(fullfile(here, sprintf('headline_vpa_results_%d.mat', D)), '-struct', 'R');
  [~, on, oe] = fileparts(out);  pr('\nwrote %s\n', [on oe]);   % file name only (no absolute path)
end

% ======================================================================================
function P = symparams(p)
  for fn = {'Ef','hf','rhof','rhoc','hc','Gc','b','L','m'}, P.(fn{1}) = sym(p.(fn{1})); end
  P.EI   = P.Ef*P.b*(P.hc^2*P.hf/2 + P.hc*P.hf^2);
  P.GA   = P.Gc*P.b*P.hc;
  P.rhoA = 2*P.rhof*P.b*P.hf + P.rhoc*P.b*P.hc;
  P.rhoI = P.rhof*P.b*(P.hc^2*P.hf/2 + P.hc*P.hf^2) + P.rhoc*P.b*P.hc^3/12;
  P.A = P.rhoI/P.EI + P.rhoA/P.GA;  P.B = P.rhoI*P.rhoA/(P.EI*P.GA);  P.C = P.rhoA/P.EI;
end

function fB = fBragg(P, n)
%  k = n pi/L in k^4 - A w^2 k^2 + B w^4 - C w^2 = 0, flexural (lower) root in w^2
  k = n*sym(pi)/P.L;  bb = P.A*k^2 + P.C;
  fB = sqrt(vpa(2*k^4/(bb + sqrt(bb^2 - 4*P.B*k^4))))/(2*sym(pi));
end

function [k1, k2, k1sq, k2sq, ks2] = waven(f, P)
  w = 2*sym(pi)*vpa(f);
  sq = sqrt(P.A^2*w^4 - 4*(P.B*w^4 - P.C*w^2));
  k1sq = (P.A*w^2 + sq)/2;  k2sq = (P.A*w^2 - sq)/2;  ks2 = w^2*P.rhoA/P.GA;
  k1 = kd(k1sq);  k2 = kd(k2sq);
end
function k2 = waven2(f, P), [~, k2] = waven(f, P); end
function k = kd(ksq), k = sqrt(ksq); if double(imag(k)) > 0, k = -k; end, end

function [a1, C1, a2, C2, S2, S1] = coeffs(f, P)
  [k1, k2, k1sq, k2sq, ks2] = waven(f, P);
  den = ks2*(k2sq - k1sq);
  a1 =  1/(2i*P.GA) * k1*(k2sq - ks2)/den;
  a2 = -1/(2i*P.GA) * k2*(k1sq - ks2)/den;
  C1 = cos(k1*P.L); C2 = cos(k2*P.L); S1 = sin(k1*P.L); S2 = sin(k2*P.L);
end

function KT = KTof(f, fr, P)
  w2 = (2*sym(pi)*vpa(f))^2;  K = P.m*(2*sym(pi)*fr)^2;
  KT = K*P.m*w2/(P.m*w2 - K);
end

function [b1, b2] = betas(f, fr, P)
  [a1, C1, a2, C2, S2, S1] = coeffs(f, P);  KT = KTof(f, fr, P);
  b1 = real((C1 + C2) + 1i*KT*(a1*S1 + a2*S2));
  b2 = real(C1*C2 + 1i*KT*(a1*S1*C2 + a2*S2*C1));
end

function v = edge(f, fr, P, which)
  [b1, b2] = betas(f, fr, P);
  switch which, case 'm', v = 1 + b1 + b2; case 'p', v = 1 - b1 + b2; case 'd', v = b1^2 - 4*b2; end
end
function v = disc(f, fr, P), v = edge(f, fr, P, 'd'); end

function v = cfgb(f, fr, P)
  [a1, ~, a2] = coeffs(f, P);  [k1, k2] = waven(f, P);
  v = real(1 - 1i*KTof(f, fr, P)*(a1*cot(k1*P.L/2) + a2*cot(k2*P.L/2)));
end

function ps = psiB(f, fr, P) %#ok<INUSD>
  [a1, C1, a2, C2, S2, S1] = coeffs(f, P);
  ps = real((a1*S1*C2 + a2*S2*C1)/(a1*S1 + a2*S2));
end

function [r, id] = edge_root(fr, P, a, b, D)
%  first band edge in (a, b): scan the three edge functions, bisect the first sign change
  fg = linspace(double(a), double(b), 400);  best = Inf;  id = '';
  for w = {'m','p','d'}
    v = arrayfun(@(f) dbl(edge(sym(f,'f'), fr, P, w{1})), fg);
    j = find(sign(v(1:end-1)) ~= sign(v(2:end)), 1);
    if ~isempty(j) && fg(j) < best, best = fg(j); jw = j; id = w{1}; end
  end
  assert(isfinite(best), 'no band edge found in the window');
  r = root_of(@(f) edge(f, fr, P, id), sym(fg(jw),'f'), sym(fg(jw+1),'f'), D);
  switch id
    case 'm', id = 'cos mu = -1 (mu = pi)';
    case 'p', id = 'cos mu = +1 (mu = 0)';
    otherwise, id = 'branch degeneracy';
  end
end

function [lo, hi, idl, idh, dlo, dhi] = band_containing(fr, fprobe, pb, P, D)
%  topology from the double production routine, values from vpa
  fp = double(fprobe);  frd = double(fr);
  e = gap_edges(@(f) beam_cmu(f, frd, pb), fp + [-0.06 0.06], ...
                struct('n', 240000, 'tol', 1e-12, 'pole', frd));
  ilo = find(e < fp, 1, 'last');  ihi = find(e > fp, 1, 'first');
  assert(~isempty(ilo) && ~isempty(ihi), 'probe is not inside an interior band');
  c = beam_cmu(fp, frd, pb);  assert(any(c >= -1 & c <= 1), 'probe is not in a pass band');
  [lo, idl] = refine(e(ilo), fr, P, D);  [hi, idh] = refine(e(ihi), fr, P, D);
  dlo = double(lo) - e(ilo);  dhi = double(hi) - e(ihi);
end

function [r, id] = refine(e0, fr, P, D)
  h = 1e-7;  a = sym(e0 - h, 'f');  b = sym(e0 + h, 'f');
  for w = {'m','p','d'}
    if sign(dbl(edge(a, fr, P, w{1}))) ~= sign(dbl(edge(b, fr, P, w{1})))
      r = root_of(@(f) edge(f, fr, P, w{1}), a, b, D);
      switch w{1}
        case 'm', id = 'cos mu = -1';  case 'p', id = 'cos mu = +1';  otherwise, id = 'degeneracy';
      end
      return
    end
  end
  error('no edge function changes sign within 1e-7 Hz of the double edge %.12f', e0);
end

function r = root_of(fun, a, b, D)
%  bisection in vpa to |b - a| < 10^(-D+8) (relative to 1 kHz scale), then one secant step
  a = vpa(a);  b = vpa(b);  fa = dbl(fun(a));  fb = dbl(fun(b));
  assert(sign(fa) ~= sign(fb), 'root not bracketed');
  tol = sym(10)^(-D + 8);
  while (b - a) > tol
    m = (a + b)/2;  fm = fun(m);
    if sign(dbl(fm)) == sign(fa), a = m; fa = dbl(fm); else, b = m; end
  end
  r = (a + b)/2;
end

function x = dbl(v), x = double(v); end
function t = s(v, n), if nargin < 2, n = 15; end, t = char(vpa(v, n)); end
function both(fid, varargin), fprintf(varargin{:}); fprintf(fid, varargin{:}); end
