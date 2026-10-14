# Fusion of local-resonance and Bragg band gaps in periodic resonator arrays

**Maria Carrillo-Munoz**<sup>a</sup>, **Bhisham Sharma**<sup>a,\*</sup>, **Jung-San Chen**<sup>b</sup>

<sup>a</sup> Department of Mechanical and Aerospace Engineering, Michigan Technological University, Houghton, MI, USA
<sup>b</sup> Department of Engineering Science, National Cheng Kung University, Tainan, Taiwan
<sup>\*</sup> Corresponding author: bnsharma@mtu.edu

This repository has the code and data that reproduce the figures, tables, and computed values of the
manuscript. The citation will be added when the article is published.

## What is reproduced

| Manuscript item | Command | Output |
|---|---|---|
| Figs. 2–10 | `run_all_figures` (in `matlab/`) | `figures/figNN_*.pdf`, `.png`, and `figures/data/figNN_*.dat` |
| Fig. 1 (schematic) | `rsvg-convert` (see below) | `figures/fig01_systems_unit_cells.pdf`, `.png` |
| Table 1 (rod stop bands) | `validate_all`, criterion 1 | printed |
| Table 2 (beam band edges, f_r = 850 Hz) | `validate_all`, criterion 4 | printed |
| Section 6.3: damping thresholds, transmission contrast, 3-dB loss factors | `damping_study` | printed; recorded in `matlab/verification/damping_study_output.txt` |
| Headline frequencies and widths at 50 digits | `headline_vpa(50)`, then `reference_values_check` | `matlab/verification/*_output*.txt` |
| Fig. 9(b), resonator-only curve | `gamma_fB1_vpa` | `matlab/verification/gamma_fB1_vpa_output.txt` |
| Fig. 10, convergence of the transmission results | `transmission_convergence`, then `transmission_profile_check` | `matlab/verification/transmission_convergence_output.txt` |
| Section 6.2, random-host tests of the general-N result | `python nchan.py`, `python nchan_transversality.py` | `python/nchannel/*_output.txt` |
| Section 7.3, exploratory two-dimensional sweep | `python coincidence_sweep.py` | `python/exploratory_2d/coincidence_sweep_output.txt` |

`validate_figure_program` (10 checks) and `validate_all` (35 checks) are the validation suites of the
figure program and the models. Both report `0 failed` when the code is intact.

## Requirements

- MATLAB R2026a with the Symbolic Math Toolbox. The toolbox is needed for Fig. 9 (100-digit arithmetic),
  both validation suites, and the scripts in `matlab/verification/`. No other toolbox is used.
- Python 3.11 with NumPy (`python/nchannel/`) and SciPy (`python/exploratory_2d/`).
- `rsvg-convert` (librsvg) only to re-export Fig. 1 from its SVG source.

The results were produced with MATLAB R2026a (Update 2), Python 3.11.16, NumPy 2.4.6, and SciPy 1.17.1.

## Running

From `matlab/`:

```matlab
run_all_figures          % Figs. 2-10: PDF, PNG and plotted data (about one minute)
validate_figure_program  % 10 checks
validate_all             % 35 checks; prints the values of Tables 1 and 2
damping_study            % damping and transmission results of Section 6.3
```

From `matlab/verification/` (each script rewrites its own transcript in that folder):

```matlab
gamma_fB1_vpa               % Fig. 9(b) resonator-only curve, 32-100 digits
headline_vpa(50)            % headline values at 50 digits (about two minutes)
reference_values_check      % double-precision values against the 50-digit results (run after headline_vpa(50))
transmission_convergence    % Fig. 10 convergence
transmission_profile_check  % Fig. 10(a) profile grid (run after transmission_convergence; appends to its transcript)
```

From `python/nchannel/` and `python/exploratory_2d/`:

```bash
python nchan.py
python nchan_transversality.py
python coincidence_sweep.py
```

Fig. 1, from `figures/` (rsvg-convert 2.62):

```bash
rsvg-convert -f pdf -o fig01_systems_unit_cells.pdf fig01_systems_unit_cells.svg
rsvg-convert -f png -d 300 -p 300 -b white -o fig01_systems_unit_cells.png fig01_systems_unit_cells.svg
```

## Checking a run against this repository

`run_all_figures` rewrites `figures/data/*.dat`; with an unchanged code base they are identical byte for byte
to the committed files, so `git status` shows no change to them. The regenerated PDF and PNG files are
pixel-identical to the committed ones but carry new embedded creation times, so git lists them as modified. The transcripts in `matlab/verification/` and `python/` are
rewritten by their scripts; apart from run stamps and timings, `git diff` shows no change. The `.mat` files
that the verification scripts write are not tracked.

## File map

| Path | Contents |
|---|---|
| `matlab/figNN_*.m` | One script per figure; `figNN` is manuscript Fig. NN. |
| `matlab/params_beam.m`, `params_rod.m`, `params_chain.m` | All model constants. Beam constants are derived from the sandwich geometry. |
| `matlab/rod_cmu.m`, `chain_dispersion.m` | Rod and chain dispersion relations (Sections 2 and 3). |
| `matlab/beam_*.m` | Timoshenko beam: wavenumbers, Green's-function coefficients, Bloch relation, transfer matrix, damping, transmission, reference values, and the extended-precision routines (Sections 4–6). |
| `matlab/gap_edges.m`, `narrow_band.m`, `unitcell_modes.m`, `compute_tuning_tolerance.m`, `damping_metrics.m`, `transmission_peaks.m` | Band edges, pass-band widths, unit-cell spectra, damping metrics, and resolved transmission maxima. |
| `matlab/journal_style.m`, `export_figure.m`, `export_fig_data.m` | Figure style and export of the PDF, PNG, and `.dat` files. |
| `matlab/run_all_figures.m`, `validate_all.m`, `validate_figure_program.m`, `damping_study.m` | Entry points (see above). |
| `matlab/verification/` | Extended-precision and convergence checks, with their output transcripts. |
| `python/nchannel/` | Random-host tests of the general-N result (Section 6.2), with transcripts. |
| `python/exploratory_2d/` | The exploratory two-dimensional sweep of Section 7.3, with its transcript. It is not used as evidence for the criterion. |
| `figures/` | The ten manuscript figures (vector PDF and 300-dpi PNG) and the SVG source of Fig. 1. |
| `figures/data/` | Every plotted number, one `.dat` file per computed figure. |

## Data file format (`figures/data/*.dat`)

Plain text. A `# figure data: <name>` line is followed by blocks, each starting with `@<key> <rows> <cols>`
and followed by that many rows, written with `%.10g`. A complex array is written as two blocks, `<key>_re`
and `<key>_im`. Text fields are written as `@<key> <text>`.

## Numerical conventions

- All model constants are set in `params_beam.m`, `params_rod.m`, and `params_chain.m`; no other routine
  contains a physical constant.
- The mode-exchange tuning, band edges, and guide lines are computed by `beam_reference_values.m` and checked
  against an independent 40-digit calculation in the validation suites.
- Lossless transmission resonances are narrower than double precision resolves, so only their maximum
  (full transmission) is reported; the plotted transmission profiles are damped.

## License

The code is released under the MIT License (`LICENSE`). The figure data files in `figures/data/` are
released under the Creative Commons Attribution 4.0 International License (`figures/data/LICENSE`). The figure
images in `figures/` (PDF, PNG, and SVG) are included for comparison with regenerated figures and are not
covered by either license.
