function transmission_convergence()
%TRANSMISSION_CONVERGENCE  Numerical convergence of the finite-specimen transmission results
%  (manuscript Fig. 10, fig10_finite_transmission.m, Section 6.3): the resolved in-band maximum
%  vs scan density, the profile grid, the stop-band floor sampling, and the 3-dB crossing.
%  The lossless in-band resonances are narrower than a plotting grid resolves, so the band
%  maximum is taken from the resolved peaks (transmission_peaks), not from the grid.
%  Writes only transmission_convergence_output.txt here.
%  Run from matlab/verification/:  transmission_convergence, then transmission_profile_check
%  (which appends part 5 to the same transcript).
  here = fileparts(mfilename('fullpath'));  root = fullfile(here, '..', '..');
  addpath(fullfile(root, 'matlab'));
  fid = fopen(fullfile(here, 'transmission_convergence_output.txt'), 'w');
  pr = @(varargin) cellfun(@(f) fprintf(f, varargin{:}), {1, fid}, 'UniformOutput', false);
  p = params_beam();  R = beam_reference_values(p);  fB = R.fB1;  band = R.coinc(1:2);
  ff = linspace(fB + 0.10, fB + 0.25, 24);  Ns = [20 50 100 200];
  pr('transmission_convergence  %s\nband [%.10f, %.10f] Hz\n', char(datetime('now')), band);

  % 1. resolved band maximum: default scan vs doubled scan
  pr('\n1. Resolved band maximum, default scan vs 2x scan\n');
  pr('   %4s %-7s %6s %6s %12s %12s %10s %9s\n','N','eta_s','peaks','peaks2','maxT dB','maxT2 dB','diff dB','C_T dB');
  worst = 0;  C0 = zeros(1, 4);
  for j = 1:4
    N = Ns(j);  fl = median(beam_transmission(ff, fB, [0 0 0], N, p));
    for e = [0 1e-6 1e-5 1e-4 1e-3]
      ns = max(2000, 20*(N + 1));
      [~, T1, i1] = transmission_peaks(fB, [0 e 0], N, band, p, ns);
      [~, T2, i2] = transmission_peaks(fB, [0 e 0], N, band, p, 2*ns + 1);
      d = 10*log10(max(T2)/max(T1));  worst = max(worst, abs(d));
      fle = median(beam_transmission(ff, fB, [0 e 0], N, p));
      pr('   %4d %-7.0e %6d %6d %12.4f %12.4f %+10.2e %9.2f\n', N, e, i1.npeaks, i2.npeaks, ...
         10*log10(max(T1)), 10*log10(max(T2)), d, 10*log10(max(T1)/fle));
      if e == 0
        C0(j) = 10*log10(max(T1)/fl);
        pr('        lossless: %d resonances, T in [%.6f, %.6f]\n', numel(T1), min(T1), max(T1));
      end
    end
  end
  pr('   worst |diff| = %.2e dB  [%s]\n', worst, pf(worst < 1e-3));
  pr('   lossless contrast C_T (Fig. 10 floor): %s dB for N_c = 20, 50, 100, 200\n', mat2str(round(C0, 1)));

  % 2. Fig. 10(a) curves (N = 100): 200-point grid vs 20000-point grid
  pr('\n2. Fig. 10(a) profiles, N_c = 100: 200-point linear interpolant vs 20000 points (dB)\n');
  f  = linspace(fB - 0.020, fB + 0.095, 200);  fF = linspace(f(1), f(end), 20000);
  out = fF < band(1) | fF > band(2);
  for e = [0 1e-5 1e-4 1e-3]
    tc = 10*log10(beam_transmission(f, fB, [0 e 0], 100, p));
    tf = 10*log10(beam_transmission(fF, fB, [0 e 0], 100, p));
    di = abs(interp1(f, tc, fF) - tf);
    if e == 0
      pr('   lossless, stop bands only: max %.3f dB\n', max(di(out)));
    else
      pr('   eta_s = %.0e: max %.3f dB (in band %.3f, stop bands %.3f)\n', e, max(di), max(di(~out)), max(di(out)));
    end
  end

  % 3. floor convention: 24 vs 240 samples
  pr('\n3. Floor: median of 24 vs 240 samples in 100-250 mHz above f_B1 (dB)\n');
  ff2 = linspace(fB + 0.10, fB + 0.25, 240);
  for N = Ns
    pr('   N = %3d: %.4f vs %.4f\n', N, 10*log10(median(beam_transmission(ff, fB, [0 1e-4 0], N, p))), ...
       10*log10(median(beam_transmission(ff2, fB, [0 1e-4 0], N, p))));
  end

  % 4. 3-dB crossing of the resolved contrast, eta solved directly (log scale)
  pr('\n4. eta_s at which the resolved C_T falls to 3 dB\n');
  for N = Ns
    CT = @(le) 10*log10(max(pk(fB, 10^le, N, band, p))/median(beam_transmission(ff, fB, [0 10^le 0], N, p))) - 3;
    le = fzero(CT, [-5 -2], optimset('TolX', 1e-6));
    pr('   N_c = %3d: eta_s = %.3e\n', N, 10^le);
  end
  fclose(fid);
end
function T = pk(fB, e, N, band, p), [~, T] = transmission_peaks(fB, [0 e 0], N, band, p); end
function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end
