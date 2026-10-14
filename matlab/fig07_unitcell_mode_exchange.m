function fig07_unitcell_mode_exchange()
%FIG07_UNITCELL_MODE_EXCHANGE  Manuscript Fig. 7: unit-cell spectra and mode exchange.
%  Writes figures/fig07_unitcell_mode_exchange.pdf and .png and the plotted data figures/data/fig07_unitcell_mode_exchange.dat.
%  Run from matlab/:  fig07_unitcell_mode_exchange   (run_all_figures makes Figs. 2-10).
  p=params_beam(); s=journal_style('colors');
  fB=beam_bragg_frequency(p); R=beam_reference_values(p); frstar=R.frstar;  % computed, not hard-coded
  frw=linspace(100,1000,61); fb=NaN(size(frw)); fc=NaN(size(frw));
  for j=1:numel(frw)
    u=unitcell_modes('beam',frw(j),[80 900],p);
    if ~isempty(u.config_b), fb(j)=u.config_b(1); end
    if ~isempty(u.config_c), fc(j)=u.config_c(1); end
  end
  frz=sort(unique([linspace(fB-0.75,fB+0.75,51) frstar fB]));
  bz=NaN(size(frz));
  for j=1:numel(frz)
    u=unitcell_modes('beam',frz(j),fB+[-0.55 0.55],p);
    if ~isempty(u.config_b), bz(j)=u.config_b(1); end
  end

  fig=journal_style('figure',7.0,5.20);
  ax1=subplot(2,2,1); set(ax1,'Position',[0.10 0.60 0.39 0.30]); hold(ax1,'on');
  plot(ax1,frw/fB,ones(size(frw)),'Color',s.ink, ...
    'LineStyle','-','LineWidth',1.0,'Marker','s','MarkerIndices',1:12:numel(frw), ...
    'MarkerSize',4,'MarkerFaceColor',s.white);
  plot(ax1,frw/fB,fb/fB,'Color',s.blue, ...
    'LineStyle','--','LineWidth',1.0,'Marker','o','MarkerIndices',5:12:numel(frw), ...
    'MarkerSize',4,'MarkerFaceColor',s.white);
  plot(ax1,frw/fB,fc/fB,'Color',s.orange, ...
    'LineStyle','-.','LineWidth',1.0,'Marker','d','MarkerIndices',9:12:numel(frw), ...
    'MarkerSize',4,'MarkerFaceColor',s.white);
  plot(ax1,frw/fB,frw/fB,'Color',s.gray, ...
    'LineStyle',':','LineWidth',0.9);
  xlabel(ax1,'\it\omega\rm_r/\it\omega\rm_{B,1}'); ylabel(ax1,'\it\omega\rm/\it\omega\rm_{B,1}');
  xlim(ax1,[100 1000]/fB); ylim(ax1,[80 900]/fB);
  journal_style('axes',ax1); journal_style('dispersion',ax1);
  journal_style('panel',ax1,'(a)');

  xz=1e3*(frz/fB-1); yz=1e3*(bz/fB-1);
  xstar=1e3*(frstar/fB-1);
  ax2=subplot(2,2,2); set(ax2,'Position',[0.61 0.60 0.34 0.30]); hold(ax2,'on');
  plot(ax2,xz,zeros(size(xz)),'Color',s.ink,'LineStyle','-','LineWidth',1.0);
  plot(ax2,xz,yz,'Color',s.blue,'LineStyle','--','LineWidth',1.0);
  xline(ax2,0,':','Color',[0.55 0.55 0.55],'LineWidth',0.8);
  plot(ax2,xstar,0,'o','Color',s.blue,'MarkerFaceColor',s.blue, ...
    'LineWidth',1.0,'MarkerSize',5);
  xlabel(ax2,'10^3(\it\omega\rm_r/\it\omega\rm_{B,1}-1)');
  ylabel(ax2,'10^3(\it\omega\rm/\it\omega\rm_{B,1}-1)');
  xlim(ax2,1e3*([-0.75 0.75]/fB)); ylim(ax2,1e3*([-0.35 0.35]/fB));
  journal_style('axes',ax2); journal_style('dispersion',ax2);
  journal_style('panel',ax2,'(b)');

  ax3=subplot(2,3,4); set(ax3,'Position',[0.05 0.05 0.27 0.36]);
  draw_cell(ax3,'a',s); journal_style('panel',ax3,'(c)');
  ax4=subplot(2,3,5); set(ax4,'Position',[0.365 0.05 0.27 0.36]);
  draw_cell(ax4,'b',s); journal_style('panel',ax4,'(d)');
  ax5=subplot(2,3,6); set(ax5,'Position',[0.68 0.05 0.27 0.36]);
  draw_cell(ax5,'c',s); journal_style('panel',ax5,'(e)');

  S.fr_wide=frw; S.end_resonators=fB;
  S.center_sliding=fb; S.center_supported=fc;
  S.fr_zoom=frz; S.center_sliding_zoom=bz; S.fr_star=frstar;
  S.note=['Unit-cell families identified by physical resonator placement ' ...
    'and boundary conditions; the zoom shows the displaced mode exchange.'];
  export_fig_data('fig07_unitcell_mode_exchange',S);
  export_figure(fig,'fig07_unitcell_mode_exchange'); close(fig);
end

function draw_cell(ax,kind,s)
  axes(ax); cla(ax); hold(ax,'on'); axis(ax,[0 1 -0.18 1.00]); axis(ax,'off');
  if kind=='a'
    col=s.ink;
    plot(ax,[0.12 0.88],[0.30 0.30],'Color',col,'LineWidth',2.0);
    support(ax,0.12,'ss'); support(ax,0.88,'ss');
    resonator(ax,0.12,s,col); resonator(ax,0.88,s,col);
    text(ax,0.5,0.93,'end resonators','Color',col, ...
      'HorizontalAlignment','center','FontSize',s.annotationFont);
    text(ax,0.5,-0.06,'supported: \itv\rm=\itM_b\rm=0','HorizontalAlignment','center', ...
      'FontSize',s.annotationFont);
  elseif kind=='b'
    col=s.blue;
    plot(ax,[0.12 0.88],[0.30 0.30],'Color',col,'LineWidth',2.0);
    support(ax,0.12,'slide'); support(ax,0.88,'slide'); resonator(ax,0.50,s,col);
    text(ax,0.5,0.93,'center resonator','Color',col, ...
      'HorizontalAlignment','center','FontSize',s.annotationFont);
    text(ax,0.5,-0.06,'sliding: \it\phi\rm=\itV\rm=0','HorizontalAlignment','center', ...
      'FontSize',s.annotationFont);
  else
    col=s.orange;
    plot(ax,[0.12 0.88],[0.30 0.30],'Color',col,'LineWidth',2.0);
    support(ax,0.12,'ss'); support(ax,0.88,'ss'); resonator(ax,0.50,s,col);
    text(ax,0.5,0.93,'center resonator','Color',col, ...
      'HorizontalAlignment','center','FontSize',s.annotationFont);
    text(ax,0.5,-0.06,'supported: \itv\rm=\itM_b\rm=0','HorizontalAlignment','center', ...
      'FontSize',s.annotationFont);
  end
end

function resonator(ax,x,s,col)
  yy=linspace(0.33,0.68,15); xx=x+0.025*[0 0 (-1).^(1:11) 0 0];
  plot(ax,xx,yy,'Color',s.gray,'LineWidth',0.9);
  rectangle('Parent',ax,'Position',[x-0.07 0.68 0.14 0.14], ...
    'FaceColor',s.white,'EdgeColor',col,'LineWidth',1.0);
end

function support(ax,x,kind)
  if strcmp(kind,'ss')
    plot(ax,[x-0.07 x x+0.07 x-0.07],[0.18 0.30 0.18 0.18], ...
      'k-','LineWidth',0.8);
    plot(ax,[x-0.09 x+0.09],[0.14 0.14],'k-','LineWidth',0.8);
  else
    plot(ax,[x x],[0.13 0.30],'k-','LineWidth',1.0);
    plot(ax,x+[-0.035 0.035],[0.13 0.13],'ko','MarkerSize',3);
  end
end
