function M = damping_metrics(fr, p, eta_vec)
%DAMPING_METRICS  Resolvability metrics of the residual passband under damping.
%
%  For a given resonator tuning fr, returns the lossless passband, the evanescent-
%  channel floor kappa*L (the value gamma_min saturates to, Eq. (6.9)), and the
%  least-decaying-wave decay gamma_min swept over three hysteretic-damping
%  directions -- eta_s (shear/core), eta_r (resonator), eta_core=eta_s=eta_r --
%  together with the adjacent stop-band shoulder gamma_wall and the resolvability
%  threshold eta_crit (Eq. (6.10), Section 6.3) under two criteria:
%     C1 (evanescent floor): gamma_min = 1/2 * kappa L
%     C2 (local wall)      : gamma_min = 1/2 * gamma_wall(+-W, same eta)
%
%  gamma_min = min over the two Bloch branches of |Im mu| = decay/cell of the
%  least-decaying wave (Section 6.3). Pure eta_r is exceptional: the nodal
%  Bragg solution is exact for complex K_T, so its continuum minimum is zero.

  if nargin < 2 || isempty(p), p = params_beam(); end
  if nargin < 3 || isempty(eta_vec), eta_vec = logspace(-6, 0, 61); end
  fB1 = beam_bragg_frequency(p);

  [lo, hi, W] = lossless_pb(fr, p);
  M.fr = fr; M.lo = lo; M.hi = hi; M.W = W;

  % evanescent-channel decay per cell at f_B1  (= saturation floor of gamma_min)
  [~, k2] = beam_wavenumbers(fB1, p);
  M.kappaL = abs(imag(k2))*p.L;

  ne = numel(eta_vec);  M.eta = eta_vec;
  M.gmin_s=zeros(1,ne); M.gmin_r=zeros(1,ne); M.gmin_core=zeros(1,ne);
  M.gwall_s=zeros(1,ne); M.gwall_r=zeros(1,ne); M.gwall_core=zeros(1,ne);
  for k = 1:ne
    e = eta_vec(k);
    dirs = {[0 e 0],[0 0 e],[0 e e]};
    [M.gmin_s(k),   M.gwall_s(k)]    = gpair(lo,hi,W,fr,dirs{1},p);
    [~,              M.gwall_r(k)]    = gpair(lo,hi,W,fr,dirs{2},p);
    M.gmin_r(k) = 0;  % nodal solution (Section 6.3): gamma(f_B1)=0 for any complex K_T
    [M.gmin_core(k),M.gwall_core(k)] = gpair(lo,hi,W,fr,dirs{3},p);
  end

  half = 0.5*M.kappaL;
  M.ncrit_floor = [first_cross(eta_vec,M.gmin_s,half), ...
                   first_cross(eta_vec,M.gmin_r,half), ...
                   first_cross(eta_vec,M.gmin_core,half)];
  M.ncrit_wall  = [first_cross(eta_vec, M.gmin_s   ./M.gwall_s,    0.5), ...
                   first_cross(eta_vec, M.gmin_r   ./M.gwall_r,    0.5), ...
                   first_cross(eta_vec, M.gmin_core./M.gwall_core, 0.5)];
  M.resonator_pin = true;
end

% ---------- helpers ----------
function [gm, gw] = gpair(lo, hi, W, fr, eta, p)
  d = hi - lo;
  f = linspace(lo, hi, 300);  f = f(f>lo+1e-4*d & f<hi-1e-4*d);
  f = f(abs(f-fr) > 1e-5*max(d,1e-4));
  gm = min(gam1(f, fr, eta, p));
  gL = gam1(lo-W, fr, eta, p);  gH = gam1(hi+W, fr, eta, p);
  gw = min(gL, gH);                        % shallower adjacent shoulder
end

function g = gam1(f, fr, eta, p), [~,~,g] = beam_cmu_d(f, fr, eta, p); end

function [lo,hi,W] = lossless_pb(fr, p)
  fB1 = beam_bragg_frequency(p); a = fB1-0.3; b = fB1+0.3;
  e = gap_edges(@(f) beam_cmu(f, fr, p), [a b], ...
                struct('n',60000,'tol',1e-10,'pole',fr));
  cuts=[a e b]; lo=NaN; hi=NaN; W=0; best=Inf;
  for i=1:numel(cuts)-1
    mid=0.5*(cuts(i)+cuts(i+1)); if abs(mid-fr)<1e-6*max(fr,1), mid=mid+1e-5; end
    c=beam_cmu(mid,fr,p);
    if any(c>=-1&c<=1)&&i>1&&i<numel(cuts)-1
      if (cuts(i+1)-cuts(i))<best, best=cuts(i+1)-cuts(i); lo=cuts(i); hi=cuts(i+1); end
    end
  end
  if isfinite(best), W=best; end
end

function xc = first_cross(x, y, yt)
  xc = NaN;
  for i=1:numel(x)-1
    if y(i) <= yt && y(i+1) > yt && y(i)~=y(i+1)
      t=(yt-y(i))/(y(i+1)-y(i));
      xc=10^(log10(x(i))+t*(log10(x(i+1))-log10(x(i)))); return
    end
  end
end
