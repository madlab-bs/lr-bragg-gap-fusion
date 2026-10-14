function reference_values_check()
%REFERENCE_VALUES_CHECK  beam_reference_values (double) against the 50-digit vpa route
%  (headline_vpa_results_50.mat, produced independently: exact rational inputs, vpa bisection).
%  Also measures the effect of the narrow_band tolerance on the fig06 widths.  Writes only
%  reference_values_check_output.txt here.
%  Run from matlab/verification/, after headline_vpa(50):  reference_values_check
  here = fileparts(mfilename('fullpath'));  root = fullfile(here, '..', '..');
  addpath(fullfile(root, 'matlab'));
  out = fullfile(here, 'reference_values_check_output.txt');  fid = fopen(out, 'w');
  pr = @(varargin) cellfun(@(f) fprintf(f, varargin{:}), {1, fid}, 'UniformOutput', false);
  p = params_beam();  tic;  R = beam_reference_values(p);  t = toc;
  V = load(fullfile(here, 'headline_vpa_results_50.mat'));
  rows = {'f_B1', R.fB1, V.fB1; 'f_r*', R.frstar, V.frstar; ...
          'coincidence hi', R.coinc(2), V.hi; 'W', R.coinc(3), V.W; ...
          'mode-exchange lo', R.mexch(1), V.lo_star; 'W*', R.mexch(3), V.Wstar; ...
          'degeneracy (200)', R.deg200, V.fdeg; 'config (b) (200)', R.slide200, V.fslide};
  pr('beam_reference_values (%.1f s) vs 50-digit vpa\n', t);
  worst = 0;
  for i = 1:size(rows, 1)
    d = rows{i,2} - str2double(rows{i,3});  worst = max(worst, abs(d));
    pr('  %-18s %.12f   vpa %s   diff %+.2e Hz\n', rows{i,1}, rows{i,2}, rows{i,3}, d);
  end
  pr('  edges: coinc {%s | %s}; mexch {%s | %s}\n', R.edge_id.coinc{:}, R.edge_id.mexch{:});
  pr('  worst |diff| = %.2e Hz  [%s]\n', worst, pf(worst < 1e-9));
  lit = [636.260501167 0.0372190 0.0362342 364.255994574 364.256921853];
  cmp = [R.frstar R.coinc(3) R.mexch(3) R.deg200 R.slide200];
  pr('\nliterals replaced (literal - computed):\n');
  nm = {'frstar','W (fig10 hi - fB)','W*','fig08 degeneracy line','fig08 sliding line'};
  for i = 1:5, pr('  %-22s %+.3e Hz\n', nm{i}, lit(i) - cmp(i)); end
  % narrow_band tolerance in compute_tuning_tolerance (fig06 panel c): W at two tolerances
  fB = R.fB1;  pr('\nfig06 widths vs narrow_band tol (fr = f_B1, f_r*):\n');
  for fr = [fB R.frstar]
    for tol = [1e-8 1e-10 1e-12]
      [lo, hi, W] = narrow_band(@(f) beam_cmu(f, fr, p), fB + [-0.3 0.3], ...
                                struct('n', 18000, 'tol', tol, 'pole', fr));
      pr('  fr=%.9f tol=%.0e  W = %.12f  (W - reference %+.2e)\n', fr, tol, W, ...
         W - (fr == fB)*R.coinc(3) - (fr ~= fB)*R.mexch(3));
    end
  end
  fclose(fid);
end
function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end
