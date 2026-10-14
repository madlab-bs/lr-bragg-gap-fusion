function transmission_profile_check()
%TRANSMISSION_PROFILE_CHECK  Fig. 10(a) damped profiles (N_c = 100): 2000-point linear
%  interpolant vs a 40000-point grid, in dB.  Read-only; appends to transmission_convergence_output.txt.
%  Run from matlab/verification/, after transmission_convergence:  transmission_profile_check
  here = fileparts(mfilename('fullpath'));  addpath(fullfile(here, '..', '..', 'matlab'));
  p = params_beam();  fB = beam_bragg_frequency(p);
  fid = fopen(fullfile(here, 'transmission_convergence_output.txt'), 'a');
  f = linspace(fB - 0.020, fB + 0.095, 2000);  fF = linspace(f(1), f(end), 40000);
  fprintf(fid, '\n5. Fig. 10(a) damped profiles, N_c = 100: 2000-point interpolant vs 40000 points\n');
  for e = [1e-5 1e-4 1e-3]
    tc = 10*log10(beam_transmission(f, fB, [0 e 0], 100, p));
    tf = 10*log10(beam_transmission(fF, fB, [0 e 0], 100, p));
    d = max(abs(interp1(f, tc, fF) - tf));
    fprintf(fid, '   eta_s = %.0e: max %.4f dB, range [%.1f, %.1f] dB\n', e, d, min(tf), max(tf));
    fprintf('eta_s = %.0e: max %.4f dB, range [%.1f, %.1f] dB\n', e, d, min(tf), max(tf));
  end
  fclose(fid);
end
