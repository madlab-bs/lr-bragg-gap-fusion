function [edges, gaps, bands] = gap_edges(cmufun, frange, opts)
%GAP_EDGES  Generic band-edge finder.
%
%   [edges, gaps, bands] = gap_edges(cmufun, [f0 f1], opts)
%
%  cmufun(fvec) must return an (nbranch x numel(fvec)) matrix of REAL cos(mu)
%  values, NaN where a branch has no real root.  A frequency is a passband iff
%  at least one entry lies in [-1, 1] -- this is the ONLY propagation test used
%  for the lossless band structures.  Gaps are never inferred from |Im mu| of a
%  selected root.
%
%  opts fields (all optional):
%     .n      scan points                          (default 40000)
%     .tol    bisection tolerance in Hz            (default 1e-9)
%     .pole   frequency of the K_T pole            (default [])
%     .branch use ONLY this branch for the test    (above the shear cutoff
%             some published band edges bound a single branch of Eq. (4.5);
%             Section 4.5)
%
%  The scan is vectorized; only the bisection is scalar.

  if nargin < 3, opts = struct(); end
  if ~isfield(opts,'n'),      opts.n      = 40000; end
  if ~isfield(opts,'tol'),    opts.tol    = 1e-9;  end
  if ~isfield(opts,'pole'),   opts.pole   = [];    end
  if ~isfield(opts,'branch'), opts.branch = [];    end

  fs = linspace(frange(1), frange(2), opts.n);
  pv = prop_vec(cmufun, fs, opts.branch);

  idx = find(diff(pv) ~= 0);
  edges = zeros(1, numel(idx));
  for j = 1:numel(idx)
    i = idx(j);
    lo = fs(i); hi = fs(i+1); plo = pv(i);
    while (hi - lo) > opts.tol
      mid = 0.5*(lo + hi);
      if prop_vec(cmufun, mid, opts.branch) == plo, lo = mid; else, hi = mid; end
    end
    edges(j) = 0.5*(lo + hi);
  end

  cuts  = [frange(1), edges, frange(2)];
  bands = []; gaps = [];
  for i = 1:numel(cuts)-1
    mid = 0.5*(cuts(i) + cuts(i+1));
    if ~isempty(opts.pole) && abs(mid - opts.pole) < 1e-6*max(opts.pole,1)
      mid = mid + 1e-5*max(opts.pole,1);
    end
    if prop_vec(cmufun, mid, opts.branch)
      bands(end+1,:) = [cuts(i), cuts(i+1)]; %#ok<AGROW>
    else
      gaps(end+1,:)  = [cuts(i), cuts(i+1)]; %#ok<AGROW>
    end
  end
end

function tf = prop_vec(cmufun, fs, branch)
  c = cmufun(fs);
  if isempty(c), tf = false(1, numel(fs)); return; end
  if size(c,2) ~= numel(fs), c = c(:).'; end
  if ~isempty(branch)
    if size(c,1) < branch, tf = false(1, numel(fs)); return; end
    c = c(branch, :);
  end
  ok = (c >= -1) & (c <= 1) & ~isnan(c);
  tf = any(ok, 1);
end
