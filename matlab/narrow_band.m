function [lo, hi, W] = narrow_band(cmufun, frange, opts)
%NARROW_BAND  Narrowest INTERIOR passband of a dispersion relation on a window.
%
%  "Interior" means bounded by band edges on both sides -- the intervals that
%  touch the ends of the window are excluded, since their width is an artifact
%  of the window.  Returns W = 0 when no interior passband exists (fully merged
%  gaps).  Used by compute_tuning_tolerance, fig06_rod_beam_coincidence and
%  validate_figure_program.

  if nargin < 3, opts = struct(); end
  if ~isfield(opts,'n'),    opts.n    = 20000; end
  if ~isfield(opts,'tol'),  opts.tol  = 1e-7;  end
  if ~isfield(opts,'pole'), opts.pole = [];    end

  e = gap_edges(cmufun, frange, opts);
  lo = NaN; hi = NaN; W = 0;
  if numel(e) < 2, return; end
  best = Inf;
  for i = 1:numel(e)-1
    mid = 0.5*(e(i) + e(i+1));
    if ~isempty(opts.pole) && abs(mid - opts.pole) < 1e-6*max(opts.pole,1)
      mid = mid + 1e-5*max(opts.pole,1);
    end
    c = cmufun(mid);
    if ~isempty(c) && any(c(:) >= -1 & c(:) <= 1)
      if (e(i+1) - e(i)) < best
        best = e(i+1) - e(i); lo = e(i); hi = e(i+1);
      end
    end
  end
  if isfinite(best), W = best; end
end
