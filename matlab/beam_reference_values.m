function R = beam_reference_values(p)
%BEAM_REFERENCE_VALUES  Computed reference frequencies of the Timoshenko beam.
%
%  The figure scripts read these computed values (f_r*, the residual pass-band edges, and
%  the f_r = 200 Hz degeneracy and sliding-cell frequencies) rather than typed literals, so
%  every figure uses one computed source.
%
%   R.fB1      k1 L = pi, Eq. (4.7)                         beam_bragg_frequency
%   R.frstar   mode exchange: beta1(f_B1) = 0, equivalently Eq. (5.2), 1 = i K_T a2 cot(alpha_2),
%              solved in closed form for K and f_r (no root search)
%   R.coinc    [lo hi W], residual pass band at f_r = f_B1: lo = f_B1 (pole limit Psi = -1,
%              Eqs. (4.9)-(4.10)); hi refined on its own edge function
%   R.mexch    [lo hi W], pass band adjacent to f_B1 at f_r = f_r*: hi = f_B1; lo refined
%   R.deg200   branch-degeneracy edge at f_r = 200 Hz (beta1^2 = 4 beta2, Eq. (5.6))
%   R.slide200 configuration (b) (sliding cell, mu = 0) at f_r = 200 Hz, Eq. (5.1)(b)
%   R.edge_id  which edge function vanishes at each refined edge
%
%  Band topology (which interval propagates) comes from gap_edges; each edge is then solved by
%  fzero to machine precision on the edge function that changes sign there (1 + b1 + b2 for
%  cos mu = -1, 1 - b1 + b2 for cos mu = +1, b1^2 - 4 b2 for a degeneracy).  Every root window is
%  asserted to contain exactly one sign change.  Cross-checked against 50-digit vpa in
%  verification/reference_values_check.m.

  if nargin < 1 || isempty(p), p = params_beam(); end
  R.fB1 = beam_bragg_frequency(p);
  fB = R.fB1;

  % ---- mode exchange, closed form ------------------------------------------------------
  [k1, k2] = beam_wavenumbers(fB, p);  [~, a2] = beam_coeffs(fB, fB, p);
  S2 = sin(k2*p.L);
  omc2 = 2*sin(k2*p.L/2)^2;                        % 1 - cos(k2 L), without cancellation
  KTs = omc2/(1i*a2*S2);                           % beta1(f_B1) = 0 with S1 = 0, C1 = -1
  assert(abs(imag(KTs)) <= 1e-12*abs(KTs), 'K_T* is not real');
  KTs = real(KTs);  wB2 = (2*pi*fB)^2;
  Ks  = KTs*p.m*wB2/(p.m*wB2 + KTs);               % K_T = K m w^2/(m w^2 - K), solved for K
  R.frstar = sqrt(Ks/p.m)/(2*pi);
  R.kL1_residual = abs(sin(k1*p.L));               % provenance: |sin k1 L| at f_B1

  % ---- residual pass band at coincidence -----------------------------------------------
  %  lo = f_B1 exactly: the K_T pole, where the finite branch is Psi = -1, Eqs. (4.9)-(4.10)
  [lo, hi, ida, idb] = band_containing(fB, fB + 0.01, p);
  assert(lo == fB, 'coincidence band does not start at f_B1');
  R.coinc = [lo hi hi - lo];  R.edge_id.coinc = {ida, idb};

  % ---- pass band adjacent to f_B1 at the mode exchange ---------------------------------
  %  hi = f_B1 exactly: S1 = 0, C1 = -1 make cos mu = -1 a root for every f_r, Eq. (4.7).  At f_r*
  %  the second root also equals +1 there (the exchange), so no single edge function is refined.
  [lo, hi, ida, idb] = band_containing(R.frstar, fB - 0.01, p);
  assert(hi == fB, 'mode-exchange band does not end at f_B1');
  R.mexch = [lo hi hi - lo];  R.edge_id.mexch = {ida, idb};

  % ---- f_r = 200 Hz: degeneracy edge and sliding-cell mode -----------------------------
  u = unitcell_modes('beam', 200, [300 400], p);
  assert(numel(u.config_b) == 1, 'expected one config-(b) mode in 300-400 Hz at f_r = 200 Hz');
  w = u.config_b + [-0.01 0.01];
  R.slide200 = solve_unique(@(f) cfg_b(f, 200, p), w);
  R.deg200   = solve_unique(@(f) edge_fun(f, 200, p, 'd'), w);
end

% ======================================================================================
function [lo, hi, ida, idb] = band_containing(fr, fprobe, p)
%  lo and hi are refined unless they are the pole or the exact k1 L = pi edge (handled by caller)
  e = gap_edges(@(f) beam_cmu(f, fr, p), fprobe + [-0.06 0.06], ...
                struct('n', 240000, 'tol', 1e-12, 'pole', fr));
  c = beam_cmu(fprobe, fr, p);
  assert(any(c >= -1 & c <= 1), 'probe frequency is not in a pass band');
  ilo = find(e < fprobe, 1, 'last');  ihi = find(e > fprobe, 1, 'first');
  assert(~isempty(ilo) && ~isempty(ihi), 'probe is not inside an interior band');
  [lo, ida] = refine(e(ilo), fr, p);  [hi, idb] = refine(e(ihi), fr, p);
end

function [r, id] = refine(e0, fr, p)
%  solve the edge function that changes sign within +-1e-7 Hz of the bracketing edge.  Two edges
%  are exact and are not refined: k1 L = pi (cos mu = -1 for every f_r, Eq. (4.7)) and the K_T pole.
  h = 1e-7;  fB = beam_bragg_frequency(p);
  if abs(e0 - fB) <= 1e-9, r = fB; id = 'cos mu = -1 (exact at k1 L = pi)'; return, end
  if abs(e0 - fr) <= 2*h,  r = e0; id = 'pole'; return, end
  ids = {'m','p','d'};  names = {'cos mu = -1','cos mu = +1','degeneracy'};
  hit = false(1, 3);
  for j = 1:3
    hit(j) = sign(edge_fun(e0 - h, fr, p, ids{j})) ~= sign(edge_fun(e0 + h, fr, p, ids{j}));
  end
  assert(nnz(hit) == 1, 'edge at %.12f Hz: %d edge functions change sign', e0, nnz(hit));
  r = solve_unique(@(f) edge_fun(f, fr, p, ids{hit}), e0 + [-h h]);  id = names{hit};
end

function r = solve_unique(fun, w)
%  exactly one sign change on a 2001-point scan of w, then fzero to machine precision
  g = linspace(w(1), w(2), 2001);  v = arrayfun(fun, g);
  k = find(sign(v(1:end-1)) ~= sign(v(2:end)));
  assert(numel(k) == 1, 'window [%.9f, %.9f] has %d sign changes, not 1', w(1), w(2), numel(k));
  r = fzero(fun, g([k k+1]), optimset('TolX', eps));
end

function v = edge_fun(f, fr, p, which)
  [~, b1, b2] = beam_cmu(f, fr, p);
  switch which
    case 'm', v = 1 + b1 + b2;
    case 'p', v = 1 - b1 + b2;
    otherwise, v = b1^2 - 4*b2;
  end
end

function v = cfg_b(f, fr, p)
%  Eq. (5.1)(b), the normalized configuration-(b) condition, real on both sides of f_c
  [k1, k2] = beam_wavenumbers(f, p);  [a1, a2, KT] = beam_coeffs(f, fr, p);
  v = real(1 - 1i*KT*(a1/tan(k1*p.L/2) + a2/tan(k2*p.L/2)));
end
