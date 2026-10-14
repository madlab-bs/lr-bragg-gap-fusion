function run_all_figures()
%RUN_ALL_FIGURES  Generate the MATLAB figures of the manuscript (Figs. 2-10).
%  figNN_*.m produces manuscript Fig. NN.  Each writes figures/figNN_*.pdf and .png and its plotted
%  data figures/data/figNN_*.dat (the .dat files are reproduced byte for byte).
%  Run from matlab/:  run_all_figures   (about one minute; Fig. 9 needs the Symbolic Math Toolbox).
%  Figure 1 is a schematic authored in SVG: its source is
%  figures/fig01_systems_unit_cells.svg, exported with
%    rsvg-convert -f pdf -o fig01_systems_unit_cells.pdf fig01_systems_unit_cells.svg
  figs={@fig02_chain_dispersion_effective_mass, ...
        @fig03_rod_bloch_tunings, ...
        @fig04_rod_gap_map_flip,@fig05_beam_complex_bloch, ...
        @fig06_rod_beam_coincidence,@fig07_unitcell_mode_exchange, ...
        @fig08_transfer_matrix_reduction,@fig09_damping_attenuation, ...
        @fig10_finite_transmission};
  for j=1:numel(figs)
    fprintf('\n[%02d/%02d] %s\n',j,numel(figs),func2str(figs{j}));
    figs{j}();
  end
end
