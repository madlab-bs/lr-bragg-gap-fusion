function files = export_figure(fig, basename, outdir)
%EXPORT_FIGURE  Export a manuscript figure as vector PDF and 300-dpi PNG.
%
%  MATLAB uses exportgraphics with vector content.  The print fallback keeps
%  the scripts compatible with Octave and older MATLAB releases.

  if nargin < 3 || isempty(outdir)
    here = fileparts(mfilename('fullpath'));
    outdir = fullfile(here, '..', 'figures');
  end
  if ~exist(outdir,'dir'), mkdir(outdir); end
  s = journal_style('colors');
  lg = findall(fig,'Type','Legend');
  if ~isempty(lg), set(lg,'Box','off','FontName','Helvetica','FontSize',s.legendFont,'FontWeight','normal'); end
  pdfname = fullfile(outdir, [basename '.pdf']);
  pngname = fullfile(outdir, [basename '.png']);
  if exist('exportgraphics','file') == 2
    exportgraphics(fig, pdfname, 'ContentType','vector', 'BackgroundColor','white');
    exportgraphics(fig, pngname, 'Resolution',300, 'BackgroundColor','white');
  else
    print(fig, pdfname, '-dpdf', '-painters');
    print(fig, pngname, '-dpng', '-r300');
  end
  files = {pdfname, pngname};
  fprintf('  wrote %s\n  wrote %s\n', pdfname, pngname);
end
