function out = unitcell_modes(model, fr, frange, p, which)
%UNITCELL_MODES  Unit-cell eigenfrequencies: beam Section 5.1, rod Section 3.2.
%
%   out = unitcell_modes(model, fr, [f0 f1], p, which)
%   model : 'beam' or 'rod'
%   which : 'all' (default) | 'a' -- config (a) only (it is resonator-free, so
%           this path skips the two transcendental scans entirely)
%
%  Three configurations, all solved analytically (no finite elements, no mesh),
%  with h1 = k1 L/2 and h2 = k2 L/2 (alpha_1 and alpha_2 of Eq. (5.1)):
%
%   (a) simply supported cell, resonators at the ENDS.  Resonator-independent.
%       beam: k1 L = n pi   Eq. (5.1)(a)        rod: k L = n pi   Eq. (3.3)
%   (b) sliding cell, resonator at the CENTER.  Bloch mu = 0.
%       beam: 1 - i K_T (a1 cot h1 + a2 cot h2) = 0   Eq. (5.1)(b)
%       rod : tan(kL/2) = K_T/(2 EA k)                Eq. (3.4), first relation
%   (c) simply supported cell, resonator at the CENTER.  Bloch mu = pi.
%       beam: 1 + i K_T (a1 tan h1 + a2 tan h2) = 0   Eq. (5.1)(c)
%       rod : tan(kL/2) = -2 EA k/K_T                 Eq. (3.4), second relation
%
%  The NORMALIZED forms (5.1)(b)/(c) are used rather than the raw half-angle
%  braces of the mu = 0 / mu = pi edge factorizations (Section 4.3), which carry
%  the extra factors sin h1 sin h2 / cos h1 cos h2.  Reason: below the shear
%  cutoff h2 is imaginary, and the raw
%  H brace  sin(h1)sin(h2) - i K_T(...)  is then purely IMAGINARY -- taking its
%  real part returns identically zero and every config-(b) root is lost.  The
%  normalized forms are real on both sides of the cutoff (a1 is imaginary and
%  a2 real below f_c; both imaginary above).  Dividing out sin/cos only removes
%  the h = n*pi/2 roots, which belong to config (a) and are found separately.

  if nargin < 5 || isempty(which), which = 'all'; end
  if nargin < 4 || isempty(p)
    if strcmpi(model,'beam'), p = params_beam(); else, p = params_rod(); end
  end
  if nargin < 3 || isempty(frange), frange = [1, 1500]; end

  out.model    = model;
  out.config_a = local_config_a(model, frange, p);
  if strcmpi(which,'a'), out.config_b = []; out.config_c = []; return; end
  out.config_b = local_roots(@(f) braceH(model, f, fr, p), frange, fr, p);
  out.config_c = local_roots(@(f) braceG(model, f, fr, p), frange, fr, p);
end

% ---------------------------------------------------------------------------
function fa = local_config_a(model, frange, p)
  fa = []; n = 1;
  while n <= 50
    if strcmpi(model,'beam')
      g = @(f) real(beam_wavenumbers(f, p))*p.L - n*pi;
    else
      g = @(f) 2*pi*f/p.c*p.L - n*pi;
    end
    fn = local_bisect_mono(g, 1e-6, 1e6, 1e-11);
    if isempty(fn), break; end
    if fn > frange(2), break; end
    if fn >= frange(1), fa(end+1) = fn; end %#ok<AGROW>
    n = n + 1;
  end
end

% ---- brace functions (vectorized, real-valued) ----------------------------
function v = braceG(model, f, fr, p)          % config (c), mu = pi
  f = f(:).';
  if strcmpi(model,'beam')
    [k1, k2] = beam_wavenumbers(f, p);
    [a1, a2, KT] = beam_coeffs(f, fr, p);
    h1 = k1*p.L/2; h2 = k2*p.L/2;
    v  = real(1 + 1i*KT.*(a1.*tan(h1) + a2.*tan(h2)));        % Eq. (5.1)(c)
  else
    k = 2*pi*f/p.c; KT = KTrod(f, fr, p);
    v = 1 + KT./(2*p.EA*k).*tan(k*p.L/2);                     % Eq. (3.4), mu = pi
  end
end

function v = braceH(model, f, fr, p)          % config (b), mu = 0
  f = f(:).';
  if strcmpi(model,'beam')
    [k1, k2] = beam_wavenumbers(f, p);
    [a1, a2, KT] = beam_coeffs(f, fr, p);
    h1 = k1*p.L/2; h2 = k2*p.L/2;
    v  = real(1 - 1i*KT.*(a1./tan(h1) + a2./tan(h2)));        % Eq. (5.1)(b)
  else
    k = 2*pi*f/p.c; KT = KTrod(f, fr, p);
    v = 1 - KT./(2*p.EA*k)./tan(k*p.L/2);                     % Eq. (3.4), mu = 0
  end
end

function KT = KTrod(f, fr, p)
  w = 2*pi*f; K = p.m*(2*pi*fr)^2;
  KT = K*p.m*w.^2./(p.m*w.^2 - K);
end

% ---- vectorized scan + scalar bisection -----------------------------------
function r = local_roots(g, frange, fr, p)
  n  = 60000;
  fs = linspace(frange(1), frange(2), n);
  v  = g(fs);
  bad = ~isfinite(v) | (abs(fs - fr) < 1e-6*max(fr,1));
  v(bad) = NaN;
  %  A sign change across a tan/cot POLE is not a root: reject it by requiring
  %  the function to be small on at least one side of the crossing.
  sgn = sign(v);
  small = min(abs(v(1:end-1)), abs(v(2:end)));
  idx = find(sgn(1:end-1).*sgn(2:end) < 0 & abs(diff(v)) < 1e3 & small < 1e2);
  r = zeros(1, numel(idx));
  for j = 1:numel(idx)
    i = idx(j);
    r(j) = local_bisect(g, fs(i), fs(i+1), p.tol_bisect);
  end
  r = r(isfinite(r));
end

function x = local_bisect(g, lo, hi, tol)
  a = g(lo);
  while (hi - lo) > tol
    mid = 0.5*(lo + hi);
    if g(mid)*a > 0, lo = mid; else, hi = mid; end
  end
  x = 0.5*(lo + hi);
end

function x = local_bisect_mono(g, lo, hi, tol)
  if g(lo) > 0 || g(hi) < 0, x = []; return; end
  while (hi - lo) > tol
    mid = 0.5*(lo + hi);
    if g(mid) < 0, lo = mid; else, hi = mid; end
  end
  x = 0.5*(lo + hi);
end
