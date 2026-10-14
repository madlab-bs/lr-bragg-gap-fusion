function p = params_chain()
%PARAMS_CHAIN  Parameters for the discrete mass-in-mass chain (Section 2).
%  Kept in its own file so that no function carries a physical constant.  The
%  defaults (M = Kb = 1, m = 0.5, w_r = 1) are a reference sanity case for the
%  chain band edges of Section 2.2, Eq. (2.5).
  p.M   = 1;      % main mass
  p.Kb  = 1;      % chain spring
  p.m   = 0.5;    % resonator mass
  p.wr  = 1;      % resonator frequency sqrt(kr/m)
  p.a   = 1;      % lattice spacing
  % Reference sanity values for the chain (Section 2.2, Eq. (2.5)): acoustic
  % top sqrt(Lambda_-), optical bottom w_0 = w_r*sqrt(1+m/M) (Eq. (2.2)) and
  % spectrum top sqrt(Lambda_+), with Lambda_-+ the roots of Eq. (2.4).
  p.sanity_acoustic_top = 0.928731;
  p.sanity_optical_bot  = 1.2247448714;   % w_0 = w_r*sqrt(1+m/M), not w_r
  p.sanity_spectrum_top = 2.153476;
end
