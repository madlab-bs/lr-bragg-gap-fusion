function p = params_beam()
%PARAMS_BEAM  Reference parameters for the resonator-loaded Timoshenko beam (Section 4).
%
%  All Timoshenko constants are DERIVED from the sandwich geometry (Table 1(a) of
%  Sharma & Sun, JSV 364 (2016) 133-146) through Eq. (4.1), rather than copied
%  from the rounded Table 1(b) figures of that paper.
%
%  NOTE: the shear rigidity is GA = Gc*b*hc, Eq. (4.1).  The 2016 paper's printed
%  Eq. (2), GA = Gc*b*(hc + 2*hf), gives 685 N and contradicts its own
%  Table 1(b) value of 625 N.  It is a misprint.
%
%  Every physical value enters through this file or params_rod.m /
%  params_chain.m; no other function carries a physical constant of its own.

  % ---- sandwich geometry and materials (paper Table 1(a)) -------------------
  p.Ef = 85.42e9;      % facesheet Young's modulus            [Pa]
  p.hf = 1.2e-3;       % facesheet thickness                  [m]
  p.rhof = 1570;       % facesheet density                    [kg/m^3]
  p.rhoc = 96.11;      % core density                         [kg/m^3]
  p.hc = 25e-3;        % core thickness                       [m]
  p.Gc = 1.0e6;        % core shear modulus                   [Pa]
  p.b  = 25e-3;        % width                                [m]

  % ---- equivalent Timoshenko constants, Eq. (4.1) --------------------------
  p.EI   = p.Ef*p.b*(0.5*p.hc^2*p.hf + p.hc*p.hf^2);
  p.GA   = p.Gc*p.b*p.hc;                      % NOT Gc*b*(hc+2hf), see above
  p.rhoA = 2*p.rhof*p.b*p.hf + p.rhoc*p.b*p.hc;
  p.rhoI = p.rhof*p.b*(0.5*p.hc^2*p.hf + p.hc*p.hf^2) + p.rhoc*p.b*p.hc^3/12;

  % ---- periodicity and resonator ------------------------------------------
  p.L = 0.05;          % resonator spacing                    [m]
  p.m = 0.025;         % resonator mass                       [kg]

  % ---- hysteretic damping loss factors (Section 6.3) ----------------------
  %  EI->EI(1+i eta_b), GA->GA(1+i eta_s), K->K(1+i eta_r) under e^{+iwt}.
  %  Default lossless; the damped functions (*_d.m) take eta=[eta_b,eta_s,eta_r]
  %  explicitly, so these do not affect the lossless suite.  Realistic ranges for
  %  a sandwich beam: viscoelastic core eta_s,eta_r ~ 0.1-1; metal facesheet
  %  eta_b ~ 1e-3-1e-2.  Band-structure resolvability threshold
  %  eta_crit ~ 3e-5 to 6e-5, Eq. (6.10).
  p.eta_b = 0;  p.eta_s = 0;  p.eta_r = 0;
  p.eta   = [p.eta_b, p.eta_s, p.eta_r];

  % ---- derived groups ------------------------------------------------------
  p.A  = p.rhoI/p.EI + p.rhoA/p.GA;            % k^4 - A w^2 k^2 + B w^4 - C w^2 = 0, Eq. (4.2)
  p.B  = p.rhoI*p.rhoA/(p.EI*p.GA);
  p.C  = p.rhoA/p.EI;
  p.fc = sqrt(p.GA/p.rhoI)/(2*pi);             % shear cutoff w_c^2 = GA/rhoI, Section 4.1
  p.theta_c = p.m/(p.rhoA*p.L);                % resonator / cell mass ratio

  % ---- cases used by the figures and by validate_all -----------------------
  p.fr_cases     = [200, 850];                 % tunings of Sharma & Sun (2016), Section 4.5
  p.fr_coincide  = [];                         % filled by beam_bragg_frequency
  p.fr_spec_coin = 636.44;                     % tuning of Fig. 8(c) of Sharma & Sun (2016)

  % ---- rounded Table 1(b) values, kept ONLY for cross-checking -------------
  p.rounded = struct('EI',877.6905,'GA',625,'rhoA',0.1543,'rhoI',1.926e-5);

  % ---- 2016 published band edges [Hz], for COMPARISON only (never targets) --
  p.published_200 = [147.84, 364.36, 635.68, 800.42, 1270.82, 1368.38];
  p.published_850 = [208.2, 636.4, 665.8, 1228];

  % ---- numerical tolerances ------------------------------------------------
  p.tol_pole   = 1e-11;   % relative half-width of the K_T pole window
  p.tol_bisect = 1e-9;    % Hz, band-edge bisection tolerance
  p.tol_imag   = 1e-8;    % |Im| below which a root counts as real
end
