function R = gamma_fB1_vpa(Dlist)
%GAMMA_FB1_VPA  gamma(f_B1) under pure resonator loss, in extended precision.
%
%  Checks the resonator-loss exception of Section 6.3: with loss only in the resonator
%  (eta_r > 0) and f_r = f_B1, the nodal Bragg state is unaffected, so gamma(f_B1) = 0
%  exactly, and any nonzero value in double precision is round-off.  The script writes
%  only gamma_fB1_vpa_output.txt and gamma_fB1_vpa_results.mat (in this folder).
%  Supports the resonator-only curve of manuscript Fig. 9(b).
%  Run from matlab/verification/:  gamma_fB1_vpa   (Symbolic Math Toolbox; under a minute).
%
%  The damped Timoshenko chain (beam_wavenumbers_d, beam_coeffs_d, the quadratic of
%  beam_cmu_d) is ported line by line to a single code path, gam_chain, that runs in
%  double or in vpa.  It has NO nodal snap and NO pole substitution (eta_r > 0 here,
%  so K_T is finite).  Parameters are the sandwich inputs of Table 1(a) of ref. [44]
%  (as in params_beam), held as exact rationals.
%
%  f_B1 (k1 L = pi) is computed two ways at each precision:
%    (i)  closed form: with k = pi/L the characteristic equation (4.2),
%         k^4 - A w^2 k^2 + B w^4 - C w^2 = 0, is a quadratic in w^2;
%    (ii) vpasolve on k1(f) L = pi, k1 from the same branch formula as beam_wavenumbers_d.
%
%  Controls (each must behave differently from the pinned state, or the port is suspect):
%    C1  double run of gam_chain reproduces beam_cmu_d off the node (generic f, mixed loss);
%    C2  host loss eta_s = 1e-4, f_r = f_B1 + 0.1 Hz (off the K_T pole): gamma(f_B1) is
%        O(1e-2) and independent of digits;
%    C2b the same at f_r = f_B1 through the pole limit Psi of Eq. (4.9), the branch that
%        beam_cmu_d substitutes there; compared with beam_cmu_d for information only;
%    C3  off-node f = f_B1 + 1 mHz, eta_r = 1e-3: gamma finite and independent of digits;
%    C4  detuned resonator f_r = f_B1 + 0.1 Hz, eta_r = 1e-3: gamma(f_B1) still falls
%        with digits (the nodal solution of Section 6.1 holds for any complex K_T, not
%        only at coincidence).

  if nargin < 1 || isempty(Dlist), Dlist = [32 50 75 100]; end
  here = fileparts(mfilename('fullpath'));
  addpath(fullfile(here, '..', '..', 'matlab'));
  p  = params_beam();
  ee = logspace(-6, -1, 51);                     % Fig. 9(b) loss factors
  out = fullfile(here, 'gamma_fB1_vpa_output.txt');
  fid = fopen(out, 'w');  cleanup = onCleanup(@() fclose(fid));
  pr = @(varargin) fprintf_both(fid, varargin{:});
  dig0 = digits;

  vs = ver('symbolic');
  pr('gamma_fB1_vpa  %s  MATLAB %s  Symbolic %s\n', char(datetime('now')), version, vs.Version);

  % ---- exact parameters ----------------------------------------------------
  P = sym_params(p);
  pr('\nExact (rational) Table 1(a) inputs and the derived constants:\n');
  for fn = {'Ef','hf','rhof','rhoc','hc','Gc','b','L','m'}
    pr('  %-5s = %-14s  |double - params_beam| = %g\n', fn{1}, char(P.(fn{1})), ...
       abs(double(P.(fn{1})) - p.(fn{1})));
  end
  for fn = {'EI','GA','rhoA','rhoI'}
    pr('  %-5s = %s  (rel. diff. from params_beam %.2g)\n', fn{1}, char(P.(fn{1})), ...
       abs(double(P.(fn{1})) - p.(fn{1}))/abs(p.(fn{1})));
  end

  % ---- double precision, no snap: the round-off-limited reference values ----
  fBd = beam_bragg_frequency(p);
  gd  = zeros(1, numel(ee));
  for k = 1:numel(ee), gd(k) = gam_chain(fBd, fBd, [0 0 ee(k)], p); end
  [~, ~, gsnap] = beam_cmu_d(fBd, fBd, [0 0 1e-3], p);
  pr('\nDouble precision, f_B1 = %.12f Hz (beam_bragg_frequency):\n', fBd);
  pr('  sin(k1 L) = %.3g;  unsnapped gamma(f_B1) in [%.3g, %.3g];  beam_cmu_d (snapped) at eta_r=1e-3: %g\n', ...
     sin(real(beam_wavenumbers(fBd, p))*p.L), min(gd), max(gd), gsnap);

  % ---- C1: port reproduces production code in double, off the node ---------
  fg = fBd + [-0.02 0.004 0.03 0.5];  etag = [1e-4 2e-5 3e-3];
  [~, ~, gref] = beam_cmu_d(fg, fBd, etag, p);
  gport = arrayfun(@(f) gam_chain(f, fBd, etag, p), fg);
  c1 = max(abs(gport - gref)./abs(gref));
  pr('\nC1 port vs beam_cmu_d (double, eta=[1e-4 2e-5 3e-3], 4 frequencies): max rel. diff = %.2g\n', c1);

  % ---- extended precision --------------------------------------------------
  nD = numel(Dlist);
  G = zeros(nD, numel(ee));  fBv = cell(1, nD);  fBgap = zeros(1, nD);  S1v = zeros(1, nD);
  c2 = zeros(1, nD);  c2b = zeros(1, nD);  c3 = zeros(1, nD);  c4 = zeros(1, nD);
  for i = 1:nD
    D = Dlist(i);  digits(D);
    [fcf, fns] = fB1_vpa(P);
    fBgap(i) = double(abs(fcf - fns)/fcf);
    fB = fcf;  fBv{i} = fB;
    k1 = kdecay(k1sq_of(fB, [0 0 0], P));
    S1v(i) = double(abs(sin(k1*P.L)));
    for k = 1:numel(ee)
      G(i, k) = double(gam_chain(fB, fB, [0 0 sym(ee(k), 'f')], P));
    end
    c2(i) = double(gam_chain(fB, fB + sym(1)/10, [0 sym(1e-4,'f') 0], P));
    c2b(i) = double(gam_chain(fB, fB, [0 sym(1e-4,'f') 0], P));
    c3(i) = double(gam_chain(fB + sym(1)/1000, fB, [0 0 sym(1e-3,'f')], P));
    c4(i) = double(gam_chain(fB, fB + sym(1)/10, [0 0 sym(1e-3,'f')], P));
    pr('\n%d digits:  f_B1 = %s Hz\n', D, char(vpa(fB, min(D, 40))));
    pr('  closed form vs vpasolve: rel. diff = %.3g;  |sin(k1 L)| = %.3g\n', fBgap(i), S1v(i));
    pr('  gamma(f_B1), pure eta_r, f_r = f_B1:  max over 51 eta = %.3g,  min = %.3g\n', ...
       max(G(i,:)), min(G(i,:)));
    pr('  C2 eta_s=1e-4, f_r=f_B1+0.1 Hz: gamma(f_B1) = %.12g\n', c2(i));
    pr('  C2b eta_s=1e-4, f_r=f_B1 (Psi limit): gamma(f_B1) = %.12g\n', c2b(i));
    pr('  C3 f_B1+1 mHz, eta_r=1e-3 gamma = %.12g\n', c3(i));
    pr('  C4 f_r=f_B1+0.1 Hz, eta_r=1e-3 gamma(f_B1) = %.3g\n', c4(i));
  end
  digits(dig0);

  pr('\nf_B1 double (beam_bragg_frequency) minus %d-digit value: %.3g Hz\n', ...
     Dlist(end), double(sym(fBd, 'f') - fBv{end}));

  pr('\nPer-eta table (gamma(f_B1), pure resonator loss, f_r = f_B1):\n');
  hdr = compose('%d digits', Dlist);
  pr('  %-12s %-12s', 'eta_r', 'double');  pr(' %-12s', hdr{:});  pr('\n');
  for k = 1:numel(ee)
    pr('  %-12.4e %-12.4e', ee(k), gd(k));  pr(' %-12.4e', G(:,k));  pr('\n');
  end

  % ---- convergence verdicts --------------------------------------------------
  pr('\nConvergence (max over eta of gamma, and log10 ratio to 10^(-D/2)):\n');
  for i = 1:nD
    pr('  D=%3d  max gamma = %.3g   log10(max gamma) + D/2 = %+.2f\n', Dlist(i), max(G(i,:)), ...
       log10(max(G(i,:))) + Dlist(i)/2);
  end
  mono = all(all(diff(log10(G), 1, 1) < -5));
  ok2  = max(abs(c2 - c2(end)))/c2(end) < 1e-12 && c2(end) > 1e-3;
  ok3  = max(abs(c3 - c3(end)))/c3(end) < 1e-12 && c3(end) > 1e-6;
  ok4  = all(diff(log10(c4)) < -5);
  pr('\nVERDICTS\n');
  pr('  [%s] gamma(f_B1) falls with digits at every eta (each step > 5 decades)\n', pf(mono));
  pr('  [%s] C1 port reproduces beam_cmu_d off the node (< 1e-10)\n', pf(c1 < 1e-10));
  pr('  [%s] C2 host loss: gamma(f_B1) = %.10g, digit-independent\n', pf(ok2), c2(end));
  [~, ~, g2b] = beam_cmu_d(fBd, fBd, [0 1e-4 0], p);
  pr('  [info] C2b Psi limit at coincidence: vpa %.10g, beam_cmu_d %.10g (rel. diff %.2g)\n', ...
     c2b(end), g2b, abs(c2b(end) - g2b)/g2b);
  pr('  [%s] C3 off-node: gamma = %.10g, digit-independent\n', pf(ok3), c3(end));
  pr('  [%s] C4 detuned f_r: gamma(f_B1) falls with digits\n', pf(ok4));
  pr('  [%s] f_B1 closed form = vpasolve to < 10^(5-D) at every D\n', ...
     pf(all(fBgap < 10.^(5 - Dlist))));

  R = struct('Dlist', Dlist, 'eta', ee, 'gamma_vpa', G, 'gamma_double_nosnap', gd, ...
             'fB1_double', fBd, 'fB1_vpa_str', {cellfun(@(x) char(vpa(x, 40)), fBv, 'uni', 0)}, ...
             'fB1_route_gap', fBgap, 'absS1', S1v, 'c1', c1, 'c2', c2, 'c2b', c2b, 'c3', c3, 'c4', c4, ...
             'verdicts', [mono, c1 < 1e-10, ok2, ok3, ok4]);
  save(fullfile(here, 'gamma_fB1_vpa_results.mat'), '-struct', 'R');
  [~, on, oe] = fileparts(out);  pr('\nwrote %s\n', [on oe]);   % file name only (no absolute path)
end

% =========================================================================
function g = gam_chain(f, fr, eta, P)
%  beam_wavenumbers_d + beam_coeffs_d + beam_cmu_d quadratic, for scalar f.
%  Works for double or sym (vpa) inputs.  No snap, no pole substitution.
  if isa(f, 'sym'), I = sym(1i); twopi = 2*sym(pi); else, I = 1i; twopi = 2*pi; end
  EIc = P.EI*(1 + I*eta(1));  GAc = P.GA*(1 + I*eta(2));
  w   = twopi*f;
  [k1sq, k2sq, ks2] = k1sq_of(f, eta, P);
  k1 = kdecay(k1sq);  k2 = kdecay(k2sq);
  den = ks2*(k2sq - k1sq);
  a1  =  1/(2*I*GAc) * k1*(k2sq - ks2)/den;
  a2  = -1/(2*I*GAc) * k2*(k1sq - ks2)/den;
  C1 = cos(k1*P.L);  C2 = cos(k2*P.L);  S1 = sin(k1*P.L);  S2 = sin(k2*P.L);
  if eta(3) == 0 && isequal(f, fr)
    % exact K_T pole: one branch diverges, the finite one is Psi, Eq. (4.9), as beam_cmu_d
    psi = (a1*S1*C2 + a2*S2*C1)/(a1*S1 + a2*S2);
    g = abs(imag(acos(psi)));  return
  end
  Kc  = P.m*(twopi*fr)^2*(1 + I*eta(3));
  KT  = Kc*P.m*w^2/(P.m*w^2 - Kc);
  b1 = (C1 + C2) + I*KT*(a1*S1 + a2*S2);
  b2 = C1*C2     + I*KT*(a1*S1*C2 + a2*S2*C1);
  d  = b1^2 - 4*b2;  sq = sqrt(d);
  sgn = 1;  if double(real(conj(b1)*sq)) < 0, sgn = -1; end
  rbig = (b1 + sgn*sq)/2;
  if double(abs(rbig)) == 0, rsmall = (b1 - sgn*sq)/2; else, rsmall = b2/rbig; end
  g = min(abs(imag(acos(rsmall))), abs(imag(acos(rbig))));
end

function [k1sq, k2sq, ks2] = k1sq_of(f, eta, P)
  if isa(f, 'sym'), I = sym(1i); twopi = 2*sym(pi); else, I = 1i; twopi = 2*pi; end
  EIc = P.EI*(1 + I*eta(1));  GAc = P.GA*(1 + I*eta(2));
  Ac = P.rhoI/EIc + P.rhoA/GAc;  Bc = P.rhoI*P.rhoA/(EIc*GAc);  Cc = P.rhoA/EIc;
  w = twopi*f;
  disc = Ac^2*w^4 - 4*(Bc*w^4 - Cc*w^2);  sq = sqrt(disc);
  k1sq = (Ac*w^2 + sq)/2;  k2sq = (Ac*w^2 - sq)/2;  ks2 = w^2*P.rhoA/GAc;
end

function k = kdecay(ksq)
  k = sqrt(ksq);  if double(imag(k)) > 0, k = -k; end
end

function [fcf, fns] = fB1_vpa(P)
%  (i) closed form.  With k = pi/L:  B W^2 - (A k^2 + C) W + k^4 = 0,  W = w^2.
%  The flexural (lower) root, in cancellation-free Vieta form.
  A = P.rhoI/P.EI + P.rhoA/P.GA;  B = P.rhoI*P.rhoA/(P.EI*P.GA);  C = P.rhoA/P.EI;
  k = sym(pi)/P.L;  bb = A*k^2 + C;
  W = vpa(2*k^4/(bb + sqrt(bb^2 - 4*B*k^4)));
  fcf = sqrt(W)/(2*sym(pi));
  % (ii) vpasolve on k1(f) L = pi along the lossless branch of beam_wavenumbers_d
  syms ff positive
  w = 2*sym(pi)*ff;
  k1sq = (A*w^2 + sqrt(A^2*w^4 - 4*(B*w^4 - C*w^2)))/2;
  fns = vpasolve(k1sq*P.L^2 == sym(pi)^2, ff, [636 637]);
  assert(isscalar(fns), 'vpasolve did not return a unique root in [636,637] Hz');
end

function P = sym_params(p)
%  Table 1(a) inputs of ref. [44] as exact rationals; derived constants as in params_beam.
  for fn = {'Ef','hf','rhof','rhoc','hc','Gc','b','L','m'}
    P.(fn{1}) = sym(p.(fn{1}));                 % 'r' flag: 1.2e-3 -> 3/2500
  end
  P.EI   = P.Ef*P.b*(P.hc^2*P.hf/2 + P.hc*P.hf^2);
  P.GA   = P.Gc*P.b*P.hc;
  P.rhoA = 2*P.rhof*P.b*P.hf + P.rhoc*P.b*P.hc;
  P.rhoI = P.rhof*P.b*(P.hc^2*P.hf/2 + P.hc*P.hf^2) + P.rhoc*P.b*P.hc^3/12;
end

function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end

function fprintf_both(fid, varargin)
  fprintf(varargin{:});  fprintf(fid, varargin{:});
end
