function nfail = validate_all()
%VALIDATE_ALL  Validation suite of the manuscript's models (35 checks).
%
%  NFAIL = VALIDATE_ALL() runs criteria 0-9 (35 checks) for the discrete chain
%  (Section 2), the rod (Section 3), the Timoshenko beam (Sections 4-5) and the
%  damped beam (Section 6.3), and returns the number of failed checks.
%  Run from matlab/:  validate_all   (expected summary: 35 passed, 0 failed;
%  needs the Symbolic Math Toolbox).
%
%  It also prints the values of the manuscript's tables: criterion 1 prints the
%  rod stop bands of Table 1 (f_r = 200, 700, 1000, 1400 Hz), and criterion 4
%  prints the beam band edges of Table 2 (f_r = 850 Hz, including the LR cut-off
%  above the shear cutoff).
%
%  Every pass/fail test is built so that it can fail: a computed quantity is
%  compared with an independent route (the edge function solved directly by
%  fzero, the 2x2 generalized eigenproblem of the chain cell, the 4x4 Timoshenko
%  transfer matrix (criteria 8, 9), the config-(a) eigenproblem with a physically
%  present resonator (criterion 6), or a 40-digit variable-precision evaluation),
%  never with a re-evaluation of the same formula.
%
%  Values published by Sharma and Sun, JSV 364 (2016) 133-146 (manuscript
%  ref. [44]) and other externally quoted reference values are printed alongside
%  as comparisons only.  Nothing here is tuned to make a check pass: where this
%  model and the 2016 paper disagree, the model's value is the target and the
%  difference is printed.
%
%  Criterion 7 also tests three misprints in the 2016 paper (its Eqs. (2), (6)
%  and (8)) through their numerical consequences.  Criterion 9 tests the
%  hysteretic-damping extension of Section 6.3: continuity with the lossless
%  code at zero loss, the damped quadratic against the damped 4x4 transfer
%  matrix, the sign of the loss factor (passivity), the evanescent floor
%  gamma_min -> kappa L of Eq. (6.9), energy conservation T + R = 1, and the
%  collapse of the residual-pass-band transmission peak under host loss.
%
%  Some printed check labels carry the authors' internal derivation and record
%  identifiers; the comment beside each check gives the corresponding manuscript
%  equation or section.

  more off
  nfail = 0; npass = 0;
  pb = params_beam(); pr = params_rod(); pc = params_chain();

  hdr('LR-Bragg validation suite');
  fprintf('  engine    : %s\n', version());
  fprintf('  beam      : EI=%.6f  GA=%.4f  rhoA=%.10f  rhoI=%.12e\n', pb.EI, pb.GA, pb.rhoA, pb.rhoI);
  fprintf('  shear cutoff f_c = %.6f Hz\n', pb.fc);
  fprintf('  rod       : c=%.6f m/s  f_B1=%.6f Hz  theta_c=%.6f  gamma1=%.8f\n', ...
          pr.c, pr.fB1, pr.theta_c, pr.gamma1);

  % =====================================================================
  hdr('CRITERION 0 -- mass-in-mass chain (Section 2)');
  ch = chain_dispersion(pc);
  %  Chain band edges (Section 2.2, Eqs. (2.4)-(2.5)).  Pass/fail is against an INDEPENDENT route,
  %  the 2x2 generalized eigenproblem of the cell at mu = pi and mu = 0; the reference sanity values
  %  quoted for the chain (pc.sanity_*) are printed as comparisons only.  The optical-branch bottom
  %  is the zero of the effective mass, w_0 = w_r*sqrt(1+m/M) (Eqs. (2.2), (2.5)), not w_r.
  wpi = chain_eig(pc, pi);  w0 = chain_eig(pc, 0);
  [npass,nfail] = chk(npass,nfail,'A0 acoustic branch top', ...
      abs(ch.w_acoustic_top - wpi(1)) < 1e-12*wpi(1), ...
      sprintf('%.10f vs eig %.10f   [reference value %.6f, d=%+.1e]', ch.w_acoustic_top, wpi(1), ...
              pc.sanity_acoustic_top, ch.w_acoustic_top - pc.sanity_acoustic_top));
  [npass,nfail] = chk(npass,nfail,'A0 spectrum top', ...
      abs(ch.w_spectrum_top - wpi(2)) < 1e-12*wpi(2), ...
      sprintf('%.10f vs eig %.10f   [reference value %.6f, d=%+.1e]', ch.w_spectrum_top, wpi(2), ...
              pc.sanity_spectrum_top, ch.w_spectrum_top - pc.sanity_spectrum_top));
  [npass,nfail] = chk(npass,nfail,'A0 optical branch bottom = wr*sqrt(1+m/M)', ...
      abs(ch.w_optical_bot - w0(2)) < 1e-12*w0(2) && w0(1) < 1e-7, ...
      sprintf('%.10f vs eig %.10f  [a quoted value of 1.0 is incorrect: the mu = 0 optical mode is w_r*sqrt(1+m/M)]', ch.w_optical_bot, w0(2)));
  %  One gap only, containing the resonator pole (Sections 2.2-2.3): the chain has no Bragg gap.
  [npass,nfail] = chk(npass,nfail,'chain has exactly ONE internal gap (Section 2.2)', ...
      ch.gap(1) < pc.wr && pc.wr < ch.gap(2), ...
      sprintf('gap (%.6f, %.6f) contains the pole wr=%.4f', ch.gap(1), ch.gap(2), pc.wr));

  % =====================================================================
  hdr('CRITERION 1 -- rod band edges: independent root check; reference table printed (comparison)');
  %  Rod stop bands at the four tunings of manuscript Table 1, from Eq. (3.2).  pr.spec_* are
  %  externally quoted 1-decimal reference edges, printed for comparison only.
  spec = {pr.spec_200, pr.spec_700, pr.spec_1000, pr.spec_1400};
  worst1 = 0;
  for i = 1:numel(pr.fr_cases)
    fr = pr.fr_cases(i);
    e  = gap_edges(@(f) rod_cmu(f, fr, pr), [1, 2999], ...
                   struct('n', 60000, 'tol', pr.tol_bisect, 'pole', fr));
    t  = spec{i};
    if numel(e) == numel(t), d = max(abs(e - t)); else, d = NaN; end
    worst1 = max(worst1, edge_root_err(@(f) rod_cmu(f, fr, pr) - 1, @(f) rod_cmu(f, fr, pr) + 1, [], e, fr));
    fprintf('    f_r=%-6g computed: %s\n', fr, fmt(e));
    fprintf('             ref. cmp: %s   max|d|=%.4f Hz  (comparison only)\n', fmt(t), d);
  end
  %  Pass/fail does not use the reference values pr.spec_*.  Every edge found by gap_edges (bisection on the propagation indicator) must be a root of
  %  cos(mu) = +1 or -1 of the rod relation (3.2), solved directly.
  [npass,nfail] = chk(npass,nfail,'rod edges are roots of cos(mu) = +-1 (reference table printed as comparison)', ...
      worst1 < 1e-8, sprintf('worst |edge - root| = %.2e Hz over 4 tunings', worst1));

  % Eq. (3.3), |cos mu| = 1 at kL = n*pi for every resonator, presumes f_r ~= n*f_B
  % (resonator pole away from the Bragg frequency).  The excluded case is coincidence,
  % where Eq. (3.8) applies instead and is tested under criterion 2.  Coincident pairs
  % are therefore skipped here, and reported so that the exclusion is visible.
  fB = pr.fB1;
  gmax = 0; nskip = 0;
  for fr = [200 700 1000 1400]
    for n = 1:3
      if abs(fr - n*fB) < 1e-9
        nskip = nskip + 1;
        fprintf('    skipped f_r=%g, n=%d : coincidence, hypothesis sin(k_r L)~=0 violated;\n', fr, n);
        fprintf('            there cos(mu) = %.6f = -(1+gamma1), tested under criterion 2\n', rod_cmu(n*fB, fr, pr));
        continue
      end
      gmax = max(gmax, abs(abs(rod_cmu(n*fB, fr, pr)) - 1));
    end
  end
  [npass,nfail] = chk(npass,nfail,'rod: kL=n*pi is a band edge for EVERY non-coincident resonator, Eq. (3.3)', ...
      gmax < 1e-12, sprintf('max||cos mu|-1| = %.3e over 11 (f_r,n) pairs; %d coincident pair skipped', gmax, nskip));

  %  Flip criterion, Eq. (3.6): K_T(f_B1) > 0 for f_r < f_B1 (gap above f_B1),
  %  K_T(f_B1) < 0 for f_r > f_B1 (gap below), with K_T from Eq. (3.1).
  KTfB = @(fr) pr.m*(2*pi*fr)^2*pr.m*(2*pi*fB)^2/(pr.m*(2*pi*fB)^2 - pr.m*(2*pi*fr)^2);
  [npass,nfail] = chk(npass,nfail,'rod flip criterion sign(K_T(f_B1)), Eq. (3.6)', ...
      KTfB(200) > 0 && KTfB(700) > 0 && KTfB(1400) < 0, ...
      sprintf('K_T(f_B1) = %+.4e / %+.4e / %+.4e for f_r=200/700/1400', KTfB(200), KTfB(700), KTfB(1400)));

  % =====================================================================
  hdr('CRITERION 2 -- rod at coincidence: no passband in 995-1005 Hz at 1e-3 Hz');
  fs = 995:0.001:1005;
  cm = rod_cmu(fs, pr.fB1, pr);
  npb = sum(abs(cm) <= 1);
  [npass,nfail] = chk(npass,nfail,'rod merged: zero propagating points in 995-1005 Hz', ...
      npb == 0, sprintf('%d of %d sample points propagate', npb, numel(fs)));
  %  Coincidence value, Eq. (3.8): cos(mu) = -(1+Gamma_1) at f = f_r = f_B1 (Gamma_1 = pr.gamma1,
  %  Eq. (3.7)).  It is tested as a LIMIT: rod_cmu at f = f_r = f_B1 returns the
  %  substituted formula, which is algebraically -(1+Gamma_1), so comparing that value would be
  %  circular.  The unsubstituted relation at f_B1(1 +- d) must converge to -(1+Gamma_1) linearly
  %  in d from both sides.
  rl = @(d) rod_cmu(pr.fB1*(1+d), pr.fB1, pr) + (1 + pr.gamma1);
  Fe = [max(abs([rl(1e-5) rl(-1e-5)])), max(abs([rl(1e-6) rl(-1e-6)])), max(abs([rl(1e-7) rl(-1e-7)]))];
  [npass,nfail] = chk(npass,nfail,'rod: cos(mu) at coincidence equals -(1+gamma1), Eq. (3.8)', ...
      Fe(3) < 1e-5 && all(Fe(2:3)./Fe(1:2) > 0.05 & Fe(2:3)./Fe(1:2) < 0.2), ...
      sprintf('|cos(mu)+1+gamma1| at d = 1e-5, 1e-6, 1e-7: %.2e, %.2e, %.2e (linear limit)', Fe));

  % =====================================================================
  hdr('CRITERION 3 -- beam, f_r = 200 Hz  (model values are the target; 2016 values are comparisons)');
  e200 = gap_edges(@(f) beam_cmu(f, 200, pb), [40, 906], ...
                   struct('n', 60000, 'tol', pb.tol_bisect, 'pole', 200));
  fprintf('    computed : %s\n', fmt(e200));
  fprintf('    2016 cmp : %s  (comparison only)\n', fmt(pb.published_200(1:4)));
  fprintf('    delta    : %s\n', fmt(abs(e200 - pb.published_200(1:4))));
  %  The target is this model's own edges, not the 2016 values (printed above for comparison).
  %  (i) REGRESSION against the model's 4-decimal reference edges (tgt3), at the rounding half-unit
  %  5e-5 Hz: this detects any unintended change of parameters or code.  A tolerance of order
  %  0.5 Hz would not; e.g. it would accept the rounded section constants of Table 1(b) of the 2016
  %  paper, which move f_B1 by 0.064 Hz.  (ii) CROSS-ROUTE: the degeneracy edge (Eq. (5.6), Section 5.4) equals the
  %  discriminant root beta1^2 = 4*beta2 solved independently (beam_reference_values).
  tgt3 = [147.9747, 364.2560, 636.4464, 799.8652];
  [Rr, ~, rok3, rmsg3] = safe_refs(pb, false);
  [npass,nfail] = chk(npass,nfail,'beam f_r=200 edges reproduce the derived values', ...
      rok3 && numel(e200) == 4 && max(abs(e200 - tgt3)) < 5e-5 && abs(e200(2) - Rr.deg200) < 1e-8, ...
      sprintf('max|d| = %.2e Hz vs record; degeneracy vs discriminant root %.2e Hz %s', ...
              max(abs(e200 - tgt3(1:min(4,numel(e200))))), abs(e200(min(2,end)) - Rr.deg200), rmsg3));

  % =====================================================================
  hdr('CRITERION 4 -- beam, f_r = 850 Hz (published 208.2 / 636.4 / 665.8)');
  e850 = gap_edges(@(f) beam_cmu(f, 850, pb), [40, 906], ...
                   struct('n', 60000, 'tol', pb.tol_bisect, 'pole', 850));
  fprintf('    computed : %s\n', fmt(e850));
  d4 = NaN;
  if numel(e850) == 3, d4 = max(abs(e850 - pb.published_850(1:3))); end
  fprintf('    2016 cmp : %s   max|d| = %.4f Hz  (comparison only)\n', fmt(pb.published_850(1:3)), d4);
  %  Pass/fail is on the model itself: each edge must be a root of its own edge function
  %  (1 - beta1 + beta2, 1 + beta1 + beta2, or the discriminant beta1^2 - 4*beta2; Eqs. (4.5)-(4.6),
  %  (5.6)).  The 2016 values (manuscript Table 2) are printed as comparisons only.
  w4 = edge_root_err(@(f) bfun(f, 850, pb, 'p'), @(f) bfun(f, 850, pb, 'm'), @(f) bfun(f, 850, pb, 'd'), e850, 850);
  [npass,nfail] = chk(npass,nfail,'beam f_r=850 edges are roots of their edge functions (2016 printed as comparison)', ...
      numel(e850) == 3 && w4 < 1e-8, sprintf('worst |edge - root| = %.2e Hz', w4));

  fLR = branch_cross(@(f) beam_cmu(f, 850, pb), -1, 1, [1150, 1300], pb.tol_bisect);
  %  LR cut-off at f_r = 850 Hz.  It lies above the shear cutoff f_c, where both wavenumbers are real
  %  and gaps are branch-wise, so it is located as a crossing of one Bloch branch through
  %  cos(mu) = -1 (Section 4.5).  The crossing must be a root of 1 + b1 + b2, and an eigenvalue of
  %  the independent 4x4 transfer matrix must give cos(mu) = -1 there; 1228 Hz (2016) is printed as
  %  a comparison.
  okLR = false;  tmLR = NaN;  rLR = NaN;
  if ~isempty(fLR)
    rLR = fzero(@(f) bfun(f, 850, pb, 'm'), fLR + [-1e-6 1e-6], optimset('TolX', eps));
    tmLR = min(abs(beam_bloch_tm(fLR, 850, [0 0 0], pb) + 1));
    okLR = fLR > pb.fc && abs(fLR - rLR) < 1e-8 && tmLR < 1e-6;
  end
  [npass,nfail] = chk(npass,nfail,'beam f_r=850 LR cut-off above the shear cutoff (single-branch edge, Section 4.5)', okLR, ...
      sprintf('branch crossing cos(mu)=-1 at %.4f Hz; |f - root| = %.1e Hz; 4x4 |cos mu + 1| = %.1e  [2016 cmp: 1228]', ...
              fLR, abs(fLR - rLR), tmLR));

  % =====================================================================
  hdr('CRITERION 5 -- beam keeps a passband at coincidence');
  %  Section 4.4, Eqs. (4.9)-(4.11): at f_r = f_B1 the beam retains a pass band whose lower edge
  %  is the geometric edge k1 L = pi.
  fB1b = bragg_frequency(pb);
  fprintf('    geometric Bragg cut-on (k1 L = pi): %.9f Hz\n', fB1b);
  ec = gap_edges(@(f) beam_cmu(f, fB1b, pb), [636.15, 636.95], ...
                 struct('n', 40000, 'tol', pb.tol_bisect, 'pole', fB1b));
  [pbLo, pbHi, W] = narrowest_band(@(f) beam_cmu(f, fB1b, pb), ec, [636.15 636.95]);
  fprintf('    exact params, f_r = f_B1 : passband %.7f - %.7f Hz (width %.7f)\n', pbLo, pbHi, W);
  [npass,nfail] = chk(npass,nfail,'beam passband survives at coincidence (width > 0)', ...
      W > 1e-4 && abs(pbLo - fB1b) < 1e-8, ...
      sprintf('width = %.7f Hz; lower edge - geometric edge = %.1e Hz (asserted < 1e-8)', W, pbLo - fB1b));

  % Comparison with an externally quoted residual-band window (636.382-636.431 Hz, not a manuscript
  % value), reproduced here with the ROUNDED Table 1(b) constants of the 2016 paper.
  pr2 = pb; pr2.EI = pb.rounded.EI; pr2.GA = pb.rounded.GA;
  pr2.rhoA = pb.rounded.rhoA; pr2.rhoI = pb.rounded.rhoI;
  pr2.A = pr2.rhoI/pr2.EI + pr2.rhoA/pr2.GA;
  pr2.B = pr2.rhoI*pr2.rhoA/(pr2.EI*pr2.GA); pr2.C = pr2.rhoA/pr2.EI;
  er = gap_edges(@(f) beam_cmu(f, pb.fr_spec_coin, pr2), [636.20, 636.70], ...
                 struct('n', 40000, 'tol', pb.tol_bisect, 'pole', pb.fr_spec_coin));
  [rLo, rHi, Wr] = narrowest_band(@(f) beam_cmu(f, pb.fr_spec_coin, pr2), er, [636.20 636.70]);
  fprintf('    rounded params, f_r = %.2f : passband %.6f - %.6f Hz (width %.6f)\n', ...
          pb.fr_spec_coin, rLo, rHi, Wr);
  %  Comparison only, no pass/fail (the window is an external number).  The printed shift from the
  %  exact-constant band shows the sensitivity of f_B1 and of the residual band to rounding the
  %  section constants; the model uses constants derived from the 2016 Table 1(a) geometry.
  fprintf('    [comparison] quoted rounded-constant window 636.382 - 636.431: computed %.6f - %.6f (shift from exact band %+.4f Hz)\n', ...
          rLo, rHi, rLo - pbLo);

  % =====================================================================
  hdr('CRITERION 6 -- unit-cell config (a) = Bragg cut-on and is resonator-free');
  uc = unitcell_modes('beam', 200, [100, 1400], pb);
  fa = uc.config_a(1);
  fprintf('    config (a) first eigenfrequency : %.9f Hz\n', fa);
  fprintf('    Bragg cut-on from beam_cmu      : %.9f Hz\n', e200(3));
  fprintf('    2016 Table 2 comparison         : 635.68 Hz  (delta %.4f Hz; comparison only)\n', abs(fa - 635.68));
  [npass,nfail] = chk(npass,nfail,'config (a) equals the computed Bragg cut-on (<0.01 Hz)', ...
      abs(fa - e200(3)) < 0.01, sprintf('|d| = %.3e Hz', abs(fa - e200(3))));

  % INDEPENDENT ROUTE.  Config (a) (Eq. (5.1)(a), Section 5.1) as a genuine eigenproblem of the
  % 4x4 transfer matrix, with a PHYSICALLY PRESENT resonator (m,K) at each end, not assumed to
  % sit at a node.  Solving real(beam_wavenumbers)*L = n*pi instead would be circular: that
  % function never receives m, K or f_r, so resonator independence would hold by construction.
  % Here the eigenfrequency must be invariant to m and f_r and equal the Bragg cut-on.
  refA = NaN; worstA = 0;
  for mm = [0.005 0.025 0.100 0.400]
    for frx = [100 400 fB1b 900 2000]
      faR = bisect_real(@(f) configA_det(f, mm, frx, pb), 636.0, 636.9, 1e-12);
      if isnan(refA), refA = faR; end
      worstA = max(worstA, abs(faR - refA));
    end
  end
  fprintf('    config (a), physical resonator (m,K) at each end : %.9f Hz\n', refA);
  [npass,nfail] = chk(npass,nfail,'config (a) resonator-free: physical (m,K) eigenproblem, invariant to m and f_r', ...
      worstA < 1e-8 && abs(refA - e200(3)) < 1e-6, ...
      sprintf('%.9f Hz; spread %.2e Hz over m in [0.005,0.4] (80x) and f_r in [100,2000]', refA, worstA));

  % Eq. (4.7): cos(mu) = -1 at k1 L = pi for every resonator.  Away from the
  % pole, machine-precision exactness is expected.
  gmaxb = 0;
  for frx = [50 200 850 2000]
    c = beam_cmu(fB1b, frx, pb);
    gmaxb = max(gmaxb, min(abs(c + 1)));
  end
  [npass,nfail] = chk(npass,nfail,'beam: cos(mu)=-1 is an EXACT root at k1L=pi, any resonator, Eq. (4.7)', ...
      gmaxb < 1e-12, sprintf('max|cos(mu)+1| = %.3e over f_r = 50..2000 Hz', gmaxb));

  % Psi, the finite root at the resonator pole (Eq. (4.9)), is the value beam_cmu substitutes at
  % f = f_r.  Comparing |Psi+1| with a bound built from the same S1, a1/a2 and S2 would be an
  % identity, so instead (i) Psi must be the LIMIT of the unsubstituted quadratic, converging
  % linearly from f_r(1 +- d), and (ii) at coincidence an independent 40-digit evaluation must
  % give Psi = -1 (Eq. (4.10)).
  [Rc, Vc, rokc, rmsgc] = safe_refs(pb, true);
  lim_ok = rokc;  lim_txt = rmsgc;
  frl = [200 850 1100 fB1b];  if isfinite(Rc.frstar), frl = [200 850 1100 Rc.frstar fB1b]; end
  for frx = frl
    cP = beam_cmu(frx, frx, pb);  psi = cP(1);
    %  At exact coincidence the pole limit and k1 L -> pi interact: a1 S1/(a2 S2) is small only
    %  within ~1e-5 Hz of f_B1, so the finite root sweeps from +1 to -1 there.  Linear
    %  convergence to Psi = -1 is observed only for d <= 1e-10, hence the smaller d at coincidence.
    ed = zeros(1, 2);  dd = [1e-6 1e-7];  if frx == fB1b, dd = [1e-10 1e-11]; end
    for q = 1:2
      cp = beam_cmu(frx*(1+dd(q)), frx, pb);  cm = beam_cmu(frx*(1-dd(q)), frx, pb);
      ed(q) = max(min(abs(cp - psi)), min(abs(cm - psi)));    % each side must converge
    end
    rt = ed(2)/ed(1);
    %  discriminating part: linear convergence (ratio ~0.1 per decade of d); magnitude sanity on
    %  the smaller d only, since the approach constant at coincidence is ~1e8 (scale ~1e-5 Hz)
    lim_ok = lim_ok && ed(2) < 1e-2 && (ed(2) < 1e-12 || (rt > 0.05 && rt < 0.2));
    lim_txt = [lim_txt sprintf(' %g:%.1e/%.1e', round(frx,3), ed)]; %#ok<AGROW>
  end
  [npass,nfail] = chk(npass,nfail,'beam: Psi is the pole limit, and Psi = -1 at coincidence, Eqs. (4.9)-(4.10)', ...
      lim_ok && Vc.psiB_err < 1e-30, ...
      sprintf('|Psi - c(f_r(1+-d))| at d=1e-6/1e-7 (coincidence 1e-10/1e-11) [f_r:err]%s; vpa |Psi+1| = %.1e', lim_txt, Vc.psiB_err));

  % =====================================================================
  hdr('CRITERION 7 -- boxed equations and the paper typo ledger');

  % Equation numbers "Eq.(2)", "Eq.(6)", "Eq.(8)" in the TYPO labels below are those of the 2016
  % paper (Sharma and Sun, ref. [44]), not of the manuscript.
  %
  % (7a) characteristic equation (4.2) satisfied by the computed wavenumbers
  % Relative residual = |sum of terms| / |largest term|.  Scaling by |k|^4 alone
  % is wrong for the evanescent root, whose k is O(1) while the omega^2 terms are
  % O(10^3-10^7).
  r7 = 0;
  for f = [100 400 636.4464381 900 1300]
    [k1,k2] = beam_wavenumbers(f, pb); w = 2*pi*f;
    for k = [k1 k2]
      t = [k^4, -pb.A*w^2*k^2, pb.B*w^4, -pb.C*w^2];
      r7 = max(r7, abs(sum(t))/max(abs(t)));
    end
  end
  [npass,nfail] = chk(npass,nfail,'characteristic equation Eq. (4.2) satisfied', ...
      r7 < 1e-12, sprintf('max relative residual %.3e (both roots, 100-1300 Hz)', r7));

  % (7b) a1,a2 of Eq. (4.4) against the two boundary conditions at the drive point, solved independently
  r7b = 0;
  for f = [200 636.4464381 850 1100]
    [k1,k2,k1s,k2s,ks2] = beam_wavenumbers(f, pb);
    [a1,a2] = beam_coeffs(f, 200, pb);
    Mbc = [ 1i*(k1s-ks2)/k1, 1i*(k2s-ks2)/k2 ; -1i*pb.GA*k1, -1i*pb.GA*k2 ];
    Bnum = Mbc \ [0; -0.5];                      % P = 1, shear jump V(0+) = -P/2 (Section 4.2)
    r7b = max(r7b, max(abs(Bnum(:) - [a1;a2])));
  end
  [npass,nfail] = chk(npass,nfail,'coefficients a1,a2 of Eq. (4.4) equal the boundary-condition solution', ...
      r7b < 1e-18, sprintf('max|difference| %.3e', r7b));

  % (7c) the 1/ks^2 factor in a1,a2 (Eq. (4.4)) is load-bearing
  eno = gap_edges(@(f) beam_cmu_noks2(f, 200, pb), [40, 906], ...
                  struct('n', 20000, 'tol', 1e-6, 'pole', 200));
  [npass,nfail] = chk(npass,nfail,'the 1/ks^2 factor is required (LR gap vanishes without it)', ...
      numel(eno) ~= 4, sprintf('without 1/ks^2 the edges become %s', fmt(eno)));

  % (7d) misprint in Eq. (2) of the 2016 paper: GA = Gc*b*(hc + 2*hf) evaluates to 685,
  % contradicting that paper's own Table 1(b) value 625 = Gc*b*hc.  Beyond the arithmetic,
  % the CONSEQUENCE is tested: GA = 685 misplaces the Bragg cut-on far from any published
  % value, while GA = 625 lands on 636.446 Hz.
  GA_printed = pb.Gc*pb.b*(pb.hc + 2*pb.hf);
  pbw = pb; pbw.GA = GA_printed;
  pbw.A = pbw.rhoI/pbw.EI + pbw.rhoA/pbw.GA;
  pbw.B = pbw.rhoI*pbw.rhoA/(pbw.EI*pbw.GA); pbw.C = pbw.rhoA/pbw.EI;
  fB_wrong = bragg_frequency(pbw); fB_right = bragg_frequency(pb);
  [npass,nfail] = chk(npass,nfail,'2016 paper, Eq. (2) misprint: printed GA=685 misplaces the Bragg cut-on', ...
      abs(GA_printed - 685) < 1e-9 && abs(fB_wrong - 636.4464) > 20 && abs(fB_right - 636.4464) < 1e-3, ...
      sprintf('GA=685 -> cut-on %.2f Hz (>20 Hz off, matches no published value); GA=625 -> %.4f Hz', fB_wrong, fB_right));

  % (7e) misprint in Eq. (6) of the 2016 paper: the rotary-inertia term needs rhoI, not rhoA
  %      (rhoI [kg m] and rhoA [kg/m] differ in units and by four orders of magnitude)
  [npass,nfail] = chk(npass,nfail,'2016 paper, Eq. (6) misprint: rotary inertia must be rhoI', ...
      abs(pb.rhoI/pb.rhoA - 1.2485e-4) < 1e-6, ...
      sprintf('rhoI/rhoA = %.6e (they differ by 4 orders of magnitude; not interchangeable)', pb.rhoI/pb.rhoA));

  % (7f) misprint in Eq. (8) of the 2016 paper: the printed omega^2 on the B term (it should
  % be omega^4, as in manuscript Eq. (4.2)) removes the shear cutoff.  The true cutoff mode
  % (k=0, f=f_c) satisfies the CORRECTED characteristic equation but leaves a large residual
  % in the printed one.  (Comparing sqrt(C/B) with sqrt(GA/rhoI) instead would be an identity
  % that cannot fail.)
  w_c = 2*pi*pb.fc;
  res_correct  = pb.B*w_c^4 - pb.C*w_c^2;        % k=0 in the corrected quartic  (~0)
  res_misprint = pb.B*w_c^2 - pb.C*w_c^2;        % k=0 in the printed quartic    (nonzero)
  [npass,nfail] = chk(npass,nfail,'2016 paper, Eq. (8) misprint: the printed omega^2 term fails the shear-cutoff mode', ...
      abs(res_correct) < 1e-6 && abs(res_misprint) > 1, ...
      sprintf('at (k=0, f_c): corrected residual %.3e (~0) vs printed %.3e (no cutoff mode)', res_correct, res_misprint));

  % (7g) intentionally absent.  Eq. (14) of the 2016 paper does carry the 1/ks^2 factor, so it
  % is not a misprint; the factor's numerical role is exercised by (7c) above.

  % (7h) cross-check of beam_cmu's Vieta root extraction (used for numerical stability, since
  % one root lies within ~1e-3 of +1) against LAPACK eig.
  % NOTE: this is NOT a test of the 2x2 reduction (5.4): only g^2 = pq enters det(M), so
  % eig(M) equals the roots of the quadratic (4.5) for ANY choice of the square root g.
  % What it checks is that the Vieta form agrees with eig.
  r7h = 0;
  for fr = [200 636.4 850]
    for f = [300 700 1100]
      [k1,k2] = beam_wavenumbers(f, pb);
      [a1,a2,KT] = beam_coeffs(f, fr, pb);
      pp = 1i*KT*a1*sin(k1*pb.L); qq = 1i*KT*a2*sin(k2*pb.L);
      Mm = [cos(k1*pb.L)+pp, sqrt(pp*qq); sqrt(pp*qq), cos(k2*pb.L)+qq];
      ev = sort(real(eig(Mm))); bl = beam_cmu(f, fr, pb);
      if numel(bl) == 2, r7h = max(r7h, max(abs(ev(:) - bl(:)))); end
    end
  end
  [npass,nfail] = chk(npass,nfail,'Vieta root-extraction (beam_cmu) agrees with LAPACK eig of the 2x2', ...
      r7h < 1e-12, sprintf('max|eig(M) - cos(mu)| = %.3e', r7h));

  % (7i) unit-cell eigenfrequencies (5.1)(a),(c) coincide with the corresponding band edges
  %  NOTE: this check must REQUIRE non-empty families.  A max([...]) over possibly-empty
  %  mins would pass silently if a family returned no roots at all.
  ucb = unitcell_modes('beam', 200, [100, 900], pb);
  okfam = ~isempty(ucb.config_a) && ~isempty(ucb.config_b) && ~isempty(ucb.config_c);
  if okfam
    d7i = max([min(abs(ucb.config_c - e200(1))), min(abs(ucb.config_a - e200(3))), ...
               min(abs(ucb.config_c - e200(4)))]);
  else
    d7i = Inf;
  end
  [npass,nfail] = chk(npass,nfail,'unit-cell configs (a),(c) equal their band edges exactly (Section 5.1)', ...
      okfam && d7i < 1e-6, sprintf('families non-empty: %d ; max|d| = %.3e Hz', okfam, d7i));

  %  DEGENERACY EDGES (Section 5.4, Eq. (5.6)).  In a two-wave-type host a passband can also
  %  end where the two Bloch branches meet, beta1^2 = 4 beta2 -- there both cos(mu) turn
  %  complex without either reaching +-1.  Mead's bounding-frequency theorem gives only the
  %  mu = 0 and mu = pi edges, so it does NOT capture these.  The rod relation, linear in
  %  cos(mu), has no discriminant and no such edge.  Classify every computed edge by which of
  %  the three conditions actually vanishes, and require a clear margin over the runner-up.
  fprintf('    edge classification (H: cos mu=+1, G: cos mu=-1, D: branch degeneracy)\n');
  ndeg = 0; nmu = 0; cmarg = Inf; ndeg200 = 0; fdeg200 = NaN;
  for fr = [200 850]
    ee = gap_edges(@(f) beam_cmu(f, fr, pb), [40, 906], ...
                   struct('n', 60000, 'tol', 1e-11, 'pole', fr));
    for i = 1:numel(ee)
      [~, bb1, bb2] = beam_cmu(ee(i), fr, pb);
      vals = [abs(1-bb1+bb2), abs(1+bb1+bb2), abs(bb1^2-4*bb2)];
      [vs, ord] = sort(vals); wq = ord(1); nm = {'H (mu=0)','G (mu=pi)','D (DEGENERACY)'};
      cmarg = min(cmarg, vs(2)/max(vs(1), realmin));          % classification margin (asserted)
      if wq == 3 && fr == 200, ndeg200 = ndeg200 + 1; fdeg200 = ee(i); end
      if wq == 3, ndeg = ndeg + 1; else, nmu = nmu + 1; end
      fprintf('      f_r=%3g  %14.9f Hz  ->  %s\n', fr, ee(i), nm{wq});
    end
  end
  [npass,nfail] = chk(npass,nfail,'every band edge is classified as mu=0, mu=pi, or a degeneracy', ...
      ndeg + nmu == 7 && cmarg > 1e3 && ndeg200 == 1 && abs(fdeg200 - Rc.deg200) < 1e-8, ...
      sprintf('%d Bloch-phase + %d degeneracy edges = 7; worst margin %.1e; f_r=200 degeneracy at %.9f Hz', ...
              nmu, ndeg, cmarg, fdeg200));

  %  Section 5.4: the f_r = 200 Hz degeneracy edge lies ~1e-3 Hz below the sliding-cell
  %  (config (b)) eigenfrequency.
  dgap = min(abs(ucb.config_b - e200(2)));
  [npass,nfail] = chk(npass,nfail,'config (b) sits just above the degeneracy edge it approximates', ...
      dgap < 1e-2 && dgap > 1e-9, ...
      sprintf('config (b) = %.9f Hz, passband edge = %.9f Hz, offset = %.2e Hz (a two-channel effect)', ...
              ucb.config_b(1), e200(2), dgap));

  %  Rod: the end-resonator family (Eq. (3.3)) and the two center-resonator families (Eq. (3.4))
  %  together reproduce every band edge, since the rod has no degeneracy edges.
  ucr = unitcell_modes('rod', 200, [1, 2999], pr);
  er  = gap_edges(@(f) rod_cmu(f, 200, pr), [1, 2999], ...
                  struct('n', 60000, 'tol', pr.tol_bisect, 'pole', 200));
  allr = sort([ucr.config_a, ucr.config_b, ucr.config_c]);
  d7j = Inf;
  if numel(allr) == numel(er), d7j = max(abs(allr - er)); end
  [npass,nfail] = chk(npass,nfail,'rod unit cell: the three families reproduce EVERY band edge', ...
      d7j < 1e-6, sprintf('%d families vs %d edges; max|d| = %.3e Hz', numel(allr), numel(er), d7j));

  % =====================================================================
  hdr('CRITERION 8 (independent) -- 4x4 transfer matrix vs the phased-array quadratic');
  % A completely independent route through the beam model: state y=[v,phi,M,V], one-cell
  % field matrix T=expm(A L) from the Timoshenko field equations (Section 4.1), point
  % resonator as the jump V(+)=V(-)+K_T v, and cos(mu)=1/2(lam+1/lam) from
  % eig(T_res T_field).  No Green's function, no a1/a2, no lattice sum, no "1 = -K_T S" sign
  % convention.  Agreement validates Eq. (4.4) including 1/ks^2 and the sign of a2, the shear
  % sign convention V = GA(v' + phi) with V(0+) = -P/2, and the e^{+i w t} wavenumber branch
  % rule Im k <= 0 (Sections 4.1-4.2) at once.
  worstTM = 0; nTM = 0;
  for fr = [200, fB1b, 850, 1500]
    for f = [80, 150, 300, 500, 700, 906.626210, 1000, 1300, 1800]
      if abs(f - fr) < 1e-6, continue; end
      w  = 2*pi*f;
      T  = expm(beam_stateA(w, pb)*pb.L);
      wr = 2*pi*fr; KTv = pb.m*w^2*wr^2/(w^2 - wr^2);
      R  = eye(4); R(4,1) = KTv;
      lam = eig(R*T);
      tm  = 0.5*(lam + 1./lam);
      [~, bb1, bb2] = beam_cmu(f, fr, pb);
      pa  = roots([1, -bb1, bb2]);
      dd  = 0;
      for it = 1:numel(tm), dd = max(dd, min(abs(tm(it) - pa))); end
      worstTM = max(worstTM, dd); nTM = nTM + 1;
    end
  end
  [npass,nfail] = chk(npass,nfail,'INDEPENDENT: 4x4 transfer-matrix cos(mu) = phased-array quadratic', ...
      worstTM < 1e-8, ...
      sprintf('max mismatch %.3e over %d (f,f_r) pairs, incl. f_B1 and the shear cutoff', worstTM, nTM));

  % =====================================================================
  hdr('CRITERION 9 -- hysteretic damping (Section 6.3)');
  %  Section 6.3: complex moduli EI(1+i eta_b), GA(1+i eta_s), K(1+i eta_r) under e^{+i w t};
  %  eta = [eta_b eta_s eta_r] is the loss-factor vector passed to the damped functions.
  fB1d = bragg_frequency(pb);

  % (9a) eta->0 CONTINUITY: the damped code reproduces the lossless beam_cmu
  % bit-for-bit wherever the lossless roots are real.  A bug in any damped
  % function (wrong complex branch, dropped factor) breaks this at O(1).
  wre = 0; wim = 0; np9 = 0;
  for fr = [200 fB1d 850 1500]
    for f = [120 150 300 500 700 950 1100]
      if abs(f - fr) < 1e-6, continue; end
      cl = beam_cmu(f, fr, pb);  cd = beam_cmu_d(f, fr, [0 0 0], pb);
      if all(~isnan(cl))
        wre = max(wre, max(abs(sort(cl(:)) - sort(real(cd(:))))));
        wim = max(wim, max(abs(imag(cd(:)))));  np9 = np9 + 1;
      end
    end
  end
  [npass,nfail] = chk(npass,nfail,'damped code at eta=0 reproduces lossless beam_cmu bit-for-bit', ...
      wre == 0 && wim < 1e-12, sprintf('worst |Re diff| = %.1e, worst |Im| = %.1e over %d pts', wre, wim, np9));

  % (9b) TWO INDEPENDENT damped routes agree on gamma: the phased-array complex
  % quadratic (beam_cmu_d) vs the 4x4 transfer matrix (beam_bloch_tm).  Damped
  % continuation of criterion 8; validates the complex beta1,beta2 machinery.
  wg = 0; ng = 0;
  DE = {[0 1e-3 0],[0 0 1e-3],[1e-3 1e-2 1e-2],[0 0.1 0.1]};
  for ie = 1:numel(DE)
    for fr = [200 fB1d 850]
      for f = [150 500 636.47 950 1200]
        if abs(f - fr) < 1e-6, continue; end
        [~,~,gq] = beam_cmu_d(f, fr, DE{ie}, pb);
        [~,~,gt] = beam_bloch_tm(f, fr, DE{ie}, pb);
        wg = max(wg, abs(gq - gt));  ng = ng + 1;
      end
    end
  end
  [npass,nfail] = chk(npass,nfail,'damped quadratic and 4x4 transfer matrix agree on gamma', ...
      wg < 1e-6, sprintf('max |gamma_quad - gamma_TM| = %.3e over %d damped (eta,f,f_r) pts', wg, ng));

  % (9c) PASSIVITY / sign (mutation).  Under e^{+iwt} with (1+i eta), eta_s>0 must
  % ATTENUATE a propagating wave (|t|^2 < lossless) and eta_s<0 must AMPLIFY it.
  % This pins the damping sign convention of Section 6.3; flipping it flips the test.
  T0p = beam_transmission(500, 200, [0 0 0],    10, pb);
  Tpp = beam_transmission(500, 200, [0 1e-2 0], 10, pb);
  Tmp = beam_transmission(500, 200, [0 -1e-2 0],10, pb);
  [npass,nfail] = chk(npass,nfail,'passivity: eta_s>0 attenuates, eta_s<0 amplifies (sign mutation)', ...
      Tpp < T0p*(1-1e-6) && Tmp > T0p*(1+1e-6), ...
      sprintf('T(eta=0)=%.5f  T(+eta_s)=%.5f  T(-eta_s)=%.5f', T0p, Tpp, Tmp));

  % (9d) EVANESCENT FLOOR.  Under heavy host loss the least-decaying Bloch wave
  % becomes the bare evanescent Timoshenko near-field: gamma_min -> kappa L =
  % acosh(C_2) at f_B1, a geometric quantity independent of eta (Eq. (6.9)).
  % PURE RESONATOR LOSS: the nodal Bragg state does not engage the lossy attachment,
  % so gamma(f_B1) = 0 exactly (Section 6.3).  beam_cmu_d enforces this nodal root in
  % double precision, so its value cannot fail; the test instead uses the unsnapped
  % variable-precision value (beam_gamma_node_vpa), which must fall with the working
  % precision (32 vs 50 digits) rather than saturate.  Host loss, in contrast, must give
  % gamma(f_B1) > 0.
  [~, k2c] = beam_wavenumbers(fB1d, pb);  kL = abs(imag(k2c))*pb.L;
  fw = linspace(636.4470, 636.4835, 60);
  [~,~,gsat] = beam_cmu_d(fw, fB1d, [0 0.1 0], pb);  gsat = min(gsat);
  g32 = beam_gamma_node_vpa(1e-4, pb, 32);  g50 = beam_gamma_node_vpa(1e-4, pb, 50);   % eta_r = 1e-4
  [~,~,gspin] = beam_cmu_d(fB1d, fB1d, [0 1e-4 0], pb);
  [npass,nfail] = chk(npass,nfail,'host-loss floor and pure-resonator nodal pinning', ...
      abs(gsat - kL)/kL < 0.05 && g32 < 1e-12 && g50 < 1e-6*g32 && gspin > 1e-3, ...
      sprintf('heavy eta_s min=%.6f vs kappaL=%.6f; gamma_B(eta_r) vpa 32/50 digits = %.1e/%.1e; gamma_B(eta_s)=%.2e', ...
              gsat, kL, g32, g50, gspin));

  % (9e) TRANSMISSION energy: N=0 gives T=1 exactly, and a lossless array conserves
  % power T+R=1 across pass and stop bands (beam_transmission, stable Bloch-eigenbasis solver).
  T0 = beam_transmission(500, 200, [0 0 0], 0, pb);  wtr = 0;
  for f = [150 300 500 950]
    for N = [3 9 15]
      [Tt, Rr] = beam_transmission(f, 200, [0 0 0], N, pb);  wtr = max(wtr, abs(Tt + Rr - 1));
    end
  end
  [npass,nfail] = chk(npass,nfail,'transmission N=0 gives T=1 and lossless energy T+R=1', ...
      abs(T0 - 1) < 1e-9 && wtr < 1e-5, sprintf('T(N=0) = %.10f ; worst |T+R-1| = %.2e', T0, wtr));

  % (9f) The residual-passband transmission PEAK (N=100 cells), located by transmission_peaks: lossless,
  % the band transmits fully (T = 1, the passivity bound, attained); with eta_s = 1e-3 the contrast
  % against the floor of Fig. 10 (median T, 100-250 mHz above f_B1) is below the 3-dB threshold.
  % The lossless resonances are only ~1e-9 to 1e-8 Hz wide, so a coarse frequency sample would miss
  % them; the peaks are therefore resolved by local refinement, and both bounds on T are asserted.
  [Rt, ~, rokt, rmsgt] = safe_refs(pb, false);  ffl = linspace(fB1d+0.10, fB1d+0.25, 24);
  bnd = Rt.coinc(1:2);  if ~rokt, bnd = [fB1d, fB1d + 0.0372]; end   % fallback window only (width, Section 4.4)
  [~, Tp0] = transmission_peaks(fB1d, [0 0 0], 100, bnd, pb);
  [~, Tp1] = transmission_peaks(fB1d, [0 1e-3 0], 100, bnd, pb);
  db0 = 10*log10(max(Tp0)/median(beam_transmission(ffl, fB1d, [0 0 0], 100, pb)));
  db1 = 10*log10(max(Tp1)/median(beam_transmission(ffl, fB1d, [0 1e-3 0], 100, pb)));
  [npass,nfail] = chk(npass,nfail,'passband transmission peak (N=100) collapses with damping', ...
      rokt && max(Tp0) > 1 - 1e-6 && max(Tp0) <= 1 + 1e-9 && db1 < 3, ...
      sprintf('lossless max T = %.8f (contrast %.1f dB); eta_s=1e-3 -> %.2f dB %s', max(Tp0), db0, db1, rmsgt));

  % =====================================================================
  hdr('SUMMARY');
  fprintf('  %d passed, %d failed\n', npass, nfail);
  if nfail == 0
    fprintf('  ALL CHECKS PASS\n');
  else
    fprintf('  SOME CHECKS FAILED -- see the [FAIL] lines above\n');
  end
end

% =========================================================================
function hdr(s)
  fprintf('\n%s\n%s\n', repmat('=',1,74), s);
end

function [np, nf] = chk(np, nf, name, cond, detail)
  if cond
    fprintf('  [PASS] %s\n', name); np = np + 1;
  else
    fprintf('  [FAIL] %s\n', name); nf = nf + 1;
  end
  if nargin >= 5 && ~isempty(detail)
    fprintf('         %s\n', detail);
  end
end

function s = fmt(v)
  s = '';
  for i = 1:numel(v), s = [s sprintf('%11.4f', v(i))]; end %#ok<AGROW>
end

function c = beam_cmu_noks2(f, fr, p)
%  beam_cmu with the 1/ks^2 factor deliberately removed, to demonstrate that it
%  is load-bearing (criterion 7).
  f = f(:).'; N = numel(f);
  [k1,k2,k1s,k2s,ks2] = beam_wavenumbers(f, p);
  den = (k2s - k1s);                    % <-- ks2 removed on purpose
  a1  =  1./(2i*p.GA).*k1.*(k2s - ks2)./den;
  a2  = -1./(2i*p.GA).*k2.*(k1s - ks2)./den;
  w = 2*pi*f; K = p.m*(2*pi*fr)^2;
  KT = K*p.m*w.^2./(p.m*w.^2 - K);
  C1 = cos(k1*p.L); C2 = cos(k2*p.L); S1 = sin(k1*p.L); S2 = sin(k2*p.L);
  b1 = real((C1+C2) + 1i*KT.*(a1.*S1 + a2.*S2));
  b2 = real(C1.*C2  + 1i*KT.*(a1.*S1.*C2 + a2.*S2.*C1));
  d = b1.^2 - 4*b2;
  c = NaN(2,N); ok = d >= 0; s = sqrt(max(d,0));
  sg = sign(b1); sg(sg==0) = 1;
  r = (b1 + sg.*s)/2; lo = r; hi = b2./r;
  c(1,ok) = min(lo(ok),hi(ok)); c(2,ok) = max(lo(ok),hi(ok));
end

function fB = bragg_frequency(p)
%  root of k1(f) L = pi  -- the geometric (config (a)) Bragg frequency f_B1, Eqs. (4.7), (5.1)(a)
  g = @(f) real(beam_wavenumbers(f, p))*p.L - pi;   % k1 L - pi
  lo = 1; hi = 5000;
  while (hi - lo) > 1e-11
    mid = 0.5*(lo+hi);
    if g(mid) < 0, lo = mid; else, hi = mid; end
  end
  fB = 0.5*(lo+hi);
end

function f = branch_cross(cmufun, target, branch, frange, tol)
%  locate where a chosen Bloch branch crosses cos(mu) = target (above the shear cutoff the gaps
%  are branch-wise, Section 4.5)
  n = 4000; fs = linspace(frange(1), frange(2), n); f = [];
  prev = NaN;
  for i = 1:n
    c = cmufun(fs(i));
    if numel(c) < branch, prev = NaN; continue; end
    v = c(branch) - target;
    if isfinite(prev) && prev*v < 0
      lo = fs(i-1); hi = fs(i);
      while (hi-lo) > tol
        mid = 0.5*(lo+hi); cc = cmufun(mid);
        if (cc(branch)-target)*prev > 0, lo = mid; else, hi = mid; end
      end
      f = 0.5*(lo+hi); return
    end
    prev = v;
  end
end

function [lo, hi, W] = narrowest_band(cmufun, edges, frange)
%  narrowest propagating interval bounded by consecutive edges
  cuts = [frange(1), edges, frange(2)];
  lo = NaN; hi = NaN; W = Inf;
  for i = 1:numel(cuts)-1
    mid = 0.5*(cuts(i)+cuts(i+1));
    c = cmufun(mid);
    if ~isempty(c) && any(c >= -1 & c <= 1)
      if (cuts(i+1)-cuts(i)) < W && i > 1 && i < numel(cuts)-1
        W = cuts(i+1)-cuts(i); lo = cuts(i); hi = cuts(i+1);
      end
    end
  end
  if ~isfinite(W), W = 0; end
end

% =========================================================================
function A = beam_stateA(w, p)
%  4x4 Timoshenko field matrix, state y=[v, phi, M, V], dy/dx = A y.
%  (time convention e^{+i w t}; field equations of Section 4.1, with M = EI phi' and V = GA(v' + phi).)
  A = [ 0,           -1,        0,   1/p.GA ;
        0,            0,   1/p.EI,        0  ;
        0, -p.rhoI*w^2,        0,        1  ;
       -p.rhoA*w^2,   0,        0,        0  ];
end

function d = configA_det(f, m, fr, p)
%  config (a): SS cell x in [0,L] (v=0, M=0 at both ends) with a PHYSICAL resonator
%  (m,K) attached at each end.  Root of d=0 is the config-(a) eigenfrequency, built
%  without assuming the resonators sit at nodes.
  w  = 2*pi*f;  T = expm(beam_stateA(w, p)*p.L);
  wr = 2*pi*fr; KTv = m*w^2*wr^2/(w^2 - wr^2);
  R  = eye(4); R(4,1) = KTv;                     % resonator jump V+ = V- + KT v
  Tc = R*T*R;                                    % resonator at x=0 and at x=L
  M2 = [Tc(1,2), Tc(1,4); Tc(3,2), Tc(3,4)];     % unknowns phi0,V0 ; v(L)=0, M(L)=0
  d  = det(M2);
end

function r = bisect_real(g, lo, hi, tol)
%  bracketed bisection for a sign-changing real function; NaN if no sign change
  flo = g(lo); fhi = g(hi);
  if sign(flo) == sign(fhi), r = NaN; return; end
  while (hi - lo) > tol
    mid = 0.5*(lo+hi); fm = g(mid);
    if fm == 0, r = mid; return; end
    if sign(fm) == sign(flo), lo = mid; flo = fm; else, hi = mid; end
  end
  r = 0.5*(lo+hi);
end

function [R, V, ok, msg] = safe_refs(pb, wantV)
%  Fail soft: if a reference routine errors, the checks that use it FAIL with the
%  reason, and the remaining checks still run.
  ok = true;  msg = '';  V = struct('psiB_err', NaN, 'W', NaN, 'Wstar', NaN, 'frstar', NaN);
  try
    R = beam_reference_values(pb);
    if wantV, V = beam_vpa_reference(R, pb, 40); end
  catch e
    ok = false;  msg = ['[reference computation failed: ' e.message ']'];
    R = struct('fB1', NaN, 'frstar', NaN, 'coinc', NaN(1,3), 'mexch', NaN(1,3), 'deg200', NaN, 'slide200', NaN);
  end
end

function w = chain_eig(pc, mu)
%  2x2 generalized eigenproblem of the mass-in-mass cell at Bloch phase mu (independent of the
%  closed-form quadratic in chain_dispersion).  Returns sorted frequencies.
  k = pc.m*pc.wr^2;
  K = [2*pc.Kb*(1 - cos(mu)) + k, -k; -k, k];  Mm = diag([pc.M, pc.m]);
  w = sqrt(max(sort(real(eig(K, Mm))), 0)).';
end

function v = bfun(f, fr, pb, which)
  [~, b1, b2] = beam_cmu(f, fr, pb);
  switch which, case 'p', v = 1 - b1 + b2; case 'm', v = 1 + b1 + b2; otherwise, v = b1^2 - 4*b2; end
end

function worst = edge_root_err(fp, fm, fd, edges, fr)
%  Each edge must be a root of one of the edge functions (sign change within +-1e-6 Hz), solved by
%  fzero; an edge on the K_T pole f = fr is checked to sit on the pole.  Returns the worst |edge-root|.
  worst = 0;  fns = {fp, fm, fd};  fns = fns(~cellfun(@isempty, fns));
  for i = 1:numel(edges)
    e = edges(i);
    if abs(e - fr) <= 1e-6*max(fr, 1), worst = max(worst, abs(e - fr)); continue, end
    best = Inf;
    for j = 1:numel(fns)
      a = e - 1e-6;  b = e + 1e-6;  va = fns{j}(a);  vb = fns{j}(b);
      if isfinite(va) && isfinite(vb) && sign(va) ~= sign(vb)
        best = min(best, abs(e - fzero(fns{j}, [a b], optimset('TolX', eps))));
      end
    end
    worst = max(worst, best);
  end
end
