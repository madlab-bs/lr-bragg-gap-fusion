function [cmu, info] = rod_cmu(f, fr, p)
%ROD_CMU  cos(mu) = F(w) for the resonator-loaded rod, Eq. (3.2) (Section 3.1).
%
%     cos(mu) = cos(k L) + K_T(w)/(2 EA k) * sin(k L),     k = w/c_0
%
%  f may be a vector.  Returns cmu of the same size.  A frequency propagates iff
%  |cmu| <= 1 (the relation is real for all real w, so there is no root
%  selection in this model).
%
%  The pole of K_T at w = w_r is handled by the analytic limit (Section 3.4):
%    - if sin(k_r L) ~= 0 the limit is +-Inf  (attenuating);
%    - if sin(k_r L)  = 0 (coincidence) the finite residue is
%      cos(k_r L)*(1 + K L/(4 EA)) = (-1)^n (1 + Gamma_n), Eq. (3.8).

  if nargin < 3 || isempty(p), p = params_rod(); end
  f = f(:).'; cmu = zeros(size(f));
  K  = p.m*(2*pi*fr)^2;
  kr = 2*pi*fr/p.c;
  gamma = K*p.L/(4*p.EA);                       % Gamma_n = K L/(4 EA), Eq. (3.7)

  for i = 1:numel(f)
    w = 2*pi*f(i); k = w/p.c; s = sin(k*p.L); cs = cos(k*p.L);
    if abs(f(i) - fr) <= p.tol_pole*max(fr,1)
      if abs(sin(kr*p.L)) <= 1e-9
        cmu(i) = cos(kr*p.L)*(1 + gamma);       % coincidence, Eq. (3.8)
      else
        cmu(i) = Inf;                           % attenuating
      end
    else
      KT = K*p.m*w^2/(p.m*w^2 - K);             % Eq. (3.1)
      cmu(i) = cs + KT/(2*p.EA*k)*s;
    end
  end

  if nargout > 1
    info.K = K; info.gamma = gamma; info.fB = p.c/(2*p.L);
    info.KT_at_fB = K*p.m*(2*pi*info.fB)^2/(p.m*(2*pi*info.fB)^2 - K);
    info.flip = info.KT_at_fB < 0;              % inverted (gap below f_B1), Eq. (3.6)
  end
end
