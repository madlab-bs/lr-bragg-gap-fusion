function p = params_rod()
%PARAMS_ROD  Reference parameters for the resonator-loaded rod (Section 3).
%  Reference case of Section 3.3.  EA/rhoA = 10000 exactly, so c_0 = 100 m/s and
%  f_B1 = 1000 Hz exactly for L = 0.05 m.

  p.EA   = 1543;       % axial rigidity                       [N]
  p.rhoA = 0.1543;     % mass per unit length                 [kg/m]
  p.L    = 0.05;       % resonator spacing                    [m]
  p.m    = 0.025;      % resonator mass                       [kg]

  p.c   = sqrt(p.EA/p.rhoA);
  p.fB1 = p.c/(2*p.L);
  p.theta_c = p.m/(p.rhoA*p.L);                % mass ratio theta_c, Eq. (3.7)
  p.gamma1  = pi^2*p.theta_c/4;                % coupling parameter Gamma_1, Eq. (3.7) ff.

  p.fr_cases = [200, 700, 1000, 1400];         % tunings of the rod table, Section 3.3

  % Reference band-edge tables [Hz], printed for comparison only (never targets)
  p.spec_200  = [173.6, 391.1, 1000.1, 1118.4, 2000.1, 2063.3];
  p.spec_700  = [305.6, 907.1, 1000.1, 1677.3, 2000.1, 2522.0];
  p.spec_1000 = [320.4, 1841.1, 2000.1, 2724.0];
  p.spec_1400 = [328.0, 1000.1, 1051.7, 1940.8, 2000.1, 2863.2];

  p.tol_pole   = 1e-11;
  p.tol_bisect = 1e-9;
end
