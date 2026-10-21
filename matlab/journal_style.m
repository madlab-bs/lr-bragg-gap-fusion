function out = journal_style(action, varargin)
%JOURNAL_STYLE  Shared JSV plotting style for every manuscript figure.
%
%  fig = journal_style('figure', width_in, height_in)
%  s   = journal_style('colors')
%        journal_style('axes', ax)
%        journal_style('dispersion', ax)
%        journal_style('panel', ax, '(a)')
%
%  Uses only base MATLAB/Octave graphics.  The palette is a restrained subset
%  of Sasha Trubetskoy's 20 simple, distinct colors.  Figures use no more than
%  four categorical hues and supplement color with line style or markers.
%  Typography and axes follow a compact mechanics-journal style.

  s.blue     = [ 67  99 216]/255;  % #4363d8
  s.orange   = [245 130  49]/255;  % #f58231
  s.teal     = [ 70 153 144]/255;  % #469990
  s.maroon   = [128   0   0]/255;  % #800000
  s.navy     = [  0   0 117]/255;  % #000075
  s.lavender = [220 190 255]/255;  % #dcbeff
  s.mint     = [170 255 195]/255;  % #aaffc3, used only as a translucent field
  s.apricot  = [255 216 177]/255;  % #ffd8b1, used only as a translucent field
  s.gray     = [ 80  80  80]/255;  % #505050, strong reference-line gray
  s.white    = [255 255 255]/255;  % #ffffff
  s.ink      = [  0   0   0]/255;  % #000000
  s.gap      = s.lavender;
  s.light    = s.white;
  s.curveWidth     = 1.0;
  s.secondaryWidth = 0.9;
  s.referenceWidth = 1.1;
  s.bandAlpha      = 0.16;
  s.fieldAlpha     = 0.10;
  % These source sizes are calibrated for the manuscript's approximately
  % 0.77 figure reduction: labels render at the 12-pt body-text scale.
  s.tickFont       = 13;
  s.labelFont      = 15.5;
  s.legendFont     = 13;
  s.panelFont      = 16;
  s.annotationFont = 14;
  % Aliases for older color names (green -> teal, purple -> maroon).
  s.green    = s.teal;
  s.purple   = s.maroon;

  switch lower(action)
    case 'colors'
      out = s;
    case 'figure'
      w = varargin{1}; h = varargin{2};
      out = figure('Color','w','Units','inches','Position',[0.5 0.5 w h], ...
                   'PaperUnits','inches','PaperPosition',[0 0 w h], ...
                   'PaperSize',[w h], 'Renderer','painters');
      % Fix the light theme so that the output does not depend on the system appearance
      % (with a dark system theme, MATLAB R2025a and later draw legend text in light gray).
      if isprop(out,'Theme'), out.Theme = 'light'; end
      set(out,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',s.tickFont, ...
              'DefaultTextFontName','Helvetica','DefaultTextFontSize',s.annotationFont, ...
              'DefaultLineLineWidth',s.curveWidth);
    case 'axes'
      ax = varargin{1};
      set(ax, 'FontName','Helvetica', 'FontSize',s.tickFont, 'FontWeight','normal', ...
              'LineWidth',0.5, 'Box','off', 'TickDir','in', ...
              'TickLength',[0.007 0.007], 'XColor',s.ink, 'YColor',s.ink, ...
              'Color',s.white, 'Layer','top');
      grid(ax,'off');
      xl=get(ax,'XLabel'); yl=get(ax,'YLabel'); zl=get(ax,'ZLabel');
      set([xl yl zl],'FontName','Helvetica','FontSize',s.labelFont,'FontWeight','normal','Color',s.ink);
      out = ax;
    case 'dispersion'
      ax = varargin{1};
      pbaspect(ax,[1.25 1 1]);
      out = ax;
    case 'panel'
      ax = varargin{1}; lab = varargin{2};
      lab = strrep(strrep(lab,'(',''),')','');
      out = text(ax, -0.040, 1.020, lab, 'Units','normalized', ...
                 'Clipping','off', 'VerticalAlignment','bottom', 'FontName','Helvetica', ...
                 'FontSize',s.panelFont, 'FontWeight','bold', 'Color',s.ink);
    otherwise
      error('journal_style:unknownAction','Unknown action: %s', action);
  end
end
