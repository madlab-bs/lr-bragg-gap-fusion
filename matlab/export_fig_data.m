function export_fig_data(name, S, outdir)
%EXPORT_FIG_DATA  Write a figure's plotted arrays to figures/data/<name>.dat
%
%  The figure scripts call this BEFORE plotting, so that (i) the figures are
%  auditable -- every plotted number is on disk next to the image -- and (ii)
%  they can be rendered in environments where Octave's text renderer is
%  unavailable.  Format: a "# figure data" header line, then for each field an
%  "@key nrows ncols" line followed by the rows (or "@key text" for strings).
%  Complex arrays are written as <key>_re and <key>_im.

  if nargin < 3 || isempty(outdir)
    here = fileparts(mfilename('fullpath'));
    outdir = fullfile(here, '..', 'figures', 'data');
  end
  if ~exist(outdir, 'dir'), mkdir(outdir); end
  fn = fullfile(outdir, [name '.dat']);
  fid = fopen(fn, 'w');
  fprintf(fid, '# figure data: %s\n', name);
  ks = fieldnames(S);
  for i = 1:numel(ks)
    v = S.(ks{i});
    if ischar(v)
      fprintf(fid, '@%s %s\n', ks{i}, v);
    else
      v = squeeze(v);
      if isvector(v), v = v(:).'; end
      if isnumeric(v) && ~isreal(v)
        % fprintf prints only the real part of a complex array: write both parts.
        put(fid, [ks{i} '_re'], real(v));  put(fid, [ks{i} '_im'], imag(v));
      else
        put(fid, ks{i}, v);
      end
    end
  end
  fclose(fid);
  fprintf('  wrote %s\n', fn);
end

function put(fid, key, v)
  fprintf(fid, '@%s %d %d\n', key, size(v,1), size(v,2));
  for r = 1:size(v,1)
    fprintf(fid, '%.10g ', v(r,:)); fprintf(fid, '\n');
  end
end
