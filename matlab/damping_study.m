function damping_study()
%DAMPING_STUDY  Prints the damping results of Section 6.3.
%
%  Does the razor-thin residual passband (W ~ 0.037 Hz at f_B1 = 636.45 Hz) survive
%  realistic structural / viscoelastic damping?  Reports eta_crit (band structure)
%  and the transmission corollary, at coincidence and at the mode-exchange point,
%  and contrasts the rod (which merges -- no passband to protect).
%  Run from matlab/:  damping_study   (prints only; its output is recorded in
%  verification/damping_study_output.txt).
  more off
  p  = params_beam();
  pr = params_rod();
  R = beam_reference_values(p);  fB1 = R.fB1;  frstar = R.frstar;   % computed, not typed
  %  Every printed result below is computed, not typed.

  fprintf('\n########################################################################\n');
  fprintf('#   DAMPING STUDY -- residual passband under hysteretic damping         #\n');
  fprintf('#   e^{+iwt}: EI->EI(1+i eta_b), GA->GA(1+i eta_s), K->K(1+i eta_r)      #\n');
  fprintf('########################################################################\n');

  Mc = damping_metrics(fB1,    p);
  Mx = damping_metrics(frstar, p);
  fprintf('\nEvanescent-channel floor  kappa*L (at f_B1) = %.7f  (= acosh(cosh kappa L))\n', Mc.kappaL);
  fprintf('gamma_min saturates to this value under heavy host damping: the least-\n');
  fprintf('decaying Bloch wave becomes the bare evanescent Timoshenko near-field.\n');

  for MM = {Mc, Mx}
    M = MM{1};
    if abs(M.fr-fB1)<1e-6, nm=sprintf('COINCIDENCE   f_r = f_B1 = %.4f Hz (W = %.4f Hz)', fB1, R.coinc(3));
    else                    nm=sprintf('MODE-EXCHANGE f_r = f_r* = %.4f Hz (W = %.4f Hz)', frstar, R.mexch(3)); end
    fprintf('\n===== %s =====\n', nm);
    fprintf('  lossless passband : %.7f - %.7f Hz  (W = %.7f Hz)\n', M.lo, M.hi, M.W);
    fprintf('  eta_crit  (resolvable while below):        eta_s     eta_s=eta_r\n');
    fprintf('    C1  gamma_min = 1/2 kappaL            : %10.3e  %10.3e\n', M.ncrit_floor([1 3]));
    fprintf('    C2  gamma_min = 1/2 wall(+-W)         : %10.3e  %10.3e\n', M.ncrit_wall([1 3]));
  end
  gpin = beam_gamma_node_vpa([1e-6 1e-3 1e-1], p, 100);
  fprintf('\n  pure eta_r: no continuum eta_crit.  gamma(f_B1), 100-digit vpa, eta_r = 1e-6/1e-3/1e-1:\n');
  fprintf('    %.1e  %.1e  %.1e   (exact value 0: nodal pinning, Section 6.3)\n', gpin);
  ec = [Mc.ncrit_floor([1 3]) Mc.ncrit_wall([1 3]) Mx.ncrit_floor([1 3]) Mx.ncrit_wall([1 3])];  % both criteria

  % ---- transmission corollary (finite N, see beam_transmission) -----------
  fprintf('\n===== TRANSMISSION corollary (finite N, lossless leads) =====\n');
  fprintf('  contrast C_T = resolved band maximum (transmission_peaks) over the median\n');
  fprintf('  transmittance 100-250 mHz above f_B1 (the Fig. 10 convention), in dB:\n');
  fr = fB1; band = R.coinc(1:2); ffl = linspace(fB1+0.10, fB1+0.25, 24);
  Ns = [20 50 100 200];
  fprintf('    %-12s', 'eta_s \\ N'); for N=Ns, fprintf('%8d', N); end; fprintf('\n');
  for es = [0 1e-5 3e-5 1e-4 3e-4 1e-3]
    fprintf('    %-12.0e', es);
    for N=Ns
      [~,Tp] = transmission_peaks(fr,[0 es 0],N,band,p);  Tpk = max(Tp);
      Tfl = median(beam_transmission(ffl,fr,[0 es 0],N,p));
      fprintf('%8.1f', 10*log10(Tpk/Tfl));
    end
    fprintf('\n');
  end
  fprintf('  (eta_s = 0: the band maximum is T = 1, the passivity bound, attained.)\n');
  fprintf('  eta_s at which C_T = 3 dB (solved):');
  for N = Ns
    CT = @(le) 10*log10(max(pk(fr,10^le,N,band,p))/median(beam_transmission(ffl,fr,[0 10^le 0],N,p))) - 3;
    fprintf('  N=%d: %.2e', N, 10^fzero(CT, [-5 -2], optimset('TolX',1e-6)));
  end
  fprintf('\n');

  % ---- rod contrast --------------------------------------------------------
  fprintf('\n===== ROD contrast =====\n');
  npb = sum(abs(rod_cmu(995:0.002:1005, pr.fB1, pr)) <= 1);
  fprintf('  rod at coincidence: %d propagating points in 995-1005 Hz -> NO passband.\n', npb);
  fprintf('  |cos mu| = 1+gamma_1 = %.4f > 1, Eq. (3.8): the gaps FUSE, so damping is a\n', 1+pr.gamma1);
  fprintf('  non-question for the rod -- there is nothing to smear out.  The survival\n');
  fprintf('  of a passband is a two-channel (count) effect; damping tests whether it\n');
  fprintf('  survives OBSERVABLY, and for the beam the answer is set by eta_crit above.\n');

  % ---- physical verdict ----------------------------------------------------
  fprintf('\n===== PHYSICAL INTERPRETATION =====\n');
  fprintf('  eta_crit = %.1e - %.1e for host/combined loss (computed above).  Pure resonator loss\n', min(ec), max(ec));
  fprintf('  has no intrinsic continuum threshold because the Bragg state is pinned.\n');
  fprintf('  Realistic host loss\n');
  fprintf('  factors for a sandwich beam:\n');
  fprintf('    viscoelastic polymer core   eta ~ 0.1 - 1      (3-4 orders ABOVE eta_crit)\n');
  fprintf('    structural foam core        eta ~ 0.01 - 0.1   (2-3 orders ABOVE)\n');
  fprintf('    metal facesheets (bending)  eta_b ~ 1e-3 - 1e-2\n');
  fprintf('    even an all-metal analogue  eta ~ 1e-3         (>1 order ABOVE)\n');
  fprintf('  => the residual passband is NOT observable for any real sandwich beam.\n');
  fprintf('  It is a property of the CONSERVATIVE complex band structure (the absence of\n');
  fprintf('  merging is an exact statement about the lossless dispersion relation);\n');
  fprintf('  under realistic damping the least-decaying wave decays at ~kappaL/cell and\n');
  fprintf('  the passband is indistinguishable from its flanks.  Observing it would need\n');
  fprintf('  eta below about %.0e, unattainable for a viscoelastic-core sandwich.\n', min(ec));
  fprintf('\n(Figures: fig09_damping_attenuation.m and fig10_finite_transmission.m.)\n');
end

function T = pk(fr, e, N, band, p), [~, T] = transmission_peaks(fr, [0 e 0], N, band, p); end
