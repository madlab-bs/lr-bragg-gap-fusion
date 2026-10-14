function [fpk, Tpk, info] = transmission_peaks(fr, eta, N, band, p, nscan)
%TRANSMISSION_PEAKS  Resolved transmission maxima inside a pass band.
%
%  A plain sampled maximum (tens of points) steps over the lossless Fabry-Perot resonances,
%  which are only 1e-9 to 1e-8 Hz wide, so every local maximum is refined here.  Used by
%  damping_study and the finite-transmission figure (fig10_finite_transmission).
%
%  [fpk, Tpk] = transmission_peaks(fr, eta, N, [lo hi], p, nscan)
%  1. Scan the closed band on nscan points (default 20(N+1), at least 2000).  A sample placed
%     exactly on the K_T pole (f = fr, eta_r = 0) is moved inside the band by 1e-9 of its width.
%  2. Every local maximum of the samples, including a maximum at either band edge, is refined
%     by nested zooming (41 points, window shrinking ~20x per pass) restricted to the band,
%     until the window is below 1e-12 f or the peak value changes by < 1e-12 relative.
%  Validity (verification/transmission_convergence.m): with loss (eta_s >= 1e-6) the
%  resonances are broadened, the peak count is stable and the band maximum agrees to < 1e-7 dB
%  under a doubled or 10x scan.  WITHOUT loss the resonances are 1e-9 to 1e-8 Hz wide, and those
%  near the band edge are narrower than double precision resolves in f: the set of peaks is NOT
%  converged and must not be plotted.  Only the lossless MAXIMUM is reliable, because the value
%  found, T = 1.000000, attains the passivity bound T <= 1.

  if nargin < 6 || isempty(nscan), nscan = max(2000, 20*(N + 1)); end
  lo = band(1);  hi = band(2);  span = hi - lo;
  T = @(f) beam_transmission(f, fr, eta, N, p);
  g = linspace(lo, hi, nscan);
  pole = eta(3) == 0 & abs(g - fr) <= 1e-12*max(fr, 1);
  g(pole) = g(pole) + sign(mean(band) - fr + eps)*1e-9*span;
  t = T(g);
  assert(all(isfinite(t)), 'non-finite transmission on the scan');
  k = 1 + find(t(2:end-1) > t(1:end-2) & t(2:end-1) >= t(3:end));
  if t(1) > t(2), k = [1 k]; end
  if t(end) > t(end-1), k = [k numel(g)]; end
  fpk = zeros(1, numel(k));  Tpk = fpk;
  for q = 1:numel(k)
    a = g(max(k(q) - 1, 1));  b = g(min(k(q) + 1, numel(g)));  best = -Inf;
    for it = 1:60
      z = linspace(a, b, 41);  tz = T(z);  [m, j] = max(tz);
      conv = (b - a) < 1e-12*z(j) || abs(m - best) <= 1e-12*m;
      best = m;  fbest = z(j);
      a = z(max(j - 1, 1));  b = z(min(j + 1, 41));
      if conv, break; end
    end
    fpk(q) = fbest;  Tpk(q) = best;
  end
  info = struct('nscan', nscan, 'npeaks', numel(k), 'band', band);
end
