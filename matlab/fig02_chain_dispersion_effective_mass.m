function fig02_chain_dispersion_effective_mass()
%FIG02_CHAIN_DISPERSION_EFFECTIVE_MASS  Manuscript Fig. 2: mass-in-mass chain, Bloch branches and effective mass.
%  Writes figures/fig02_chain_dispersion_effective_mass.pdf and .png and the plotted data figures/data/fig02_chain_dispersion_effective_mass.dat.
%  Run from matlab/:  fig02_chain_dispersion_effective_mass   (run_all_figures makes Figs. 2-10).
  p=params_chain(); q=linspace(0,pi,500); d=chain_dispersion(p,q);
  w=[linspace(0.02,p.wr-0.004,600), linspace(p.wr+0.004,2.35,900)];
  me=d.Meff(w)/p.M; s=journal_style('colors');
  gapGray=[0.90 0.90 0.90];
  guideGray=[0.70 0.70 0.70];
  fig=journal_style('figure',7.0,3.45);
  ax1=subplot(1,2,1); set(ax1,'Position',[0.10 0.18 0.37 0.70]); hold(ax1,'on');
  patch(ax1,[0 1 1 0],[d.gap(1) d.gap(1) d.gap(2) d.gap(2)], ...
    gapGray,'EdgeColor','none','FaceAlpha',1);
  % Black/teal color hierarchy; branch identity is encoded twice:
  % acoustic = black solid; optical = teal dash-dot.  Reference lines use
  % different patterns so no two functional classes share the same style.
  plot(ax1,q/pi,d.w_lower,'Color',s.ink,'LineStyle','-', ...
    'LineWidth',1.1);
  plot(ax1,q/pi,d.w_upper,'Color',s.teal,'LineStyle','-.', ...
    'LineWidth',1.1);
  yline(ax1,p.wr,':','Color',s.gray,'LineWidth',s.referenceWidth);
  xlabel(ax1,'\it\mu\rm/\it\pi'); ylabel(ax1,'\it\omega\rm/\it\omega\rm_r');
  xlim(ax1,[0 1]); ylim(ax1,[0 2.3]);
  journal_style('axes',ax1); journal_style('dispersion',ax1); journal_style('panel',ax1,'(a)');
  text(ax1,0.60,0.80,'acoustic','Color',s.ink,'FontSize',s.annotationFont, ...
    'HorizontalAlignment','center','VerticalAlignment','top');
  text(ax1,0.62,1.96,'optical','Color',s.teal,'FontSize',s.annotationFont, ...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
  % Sit the label in the part of the gap field above the dotted omega_r guide.
  text(ax1,0.50,p.wr+0.035,'LR gap','Color',s.gray, ...
    'FontSize',s.annotationFont,'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom','Tag','lrgap');

  ax2=subplot(1,2,2); set(ax2,'Position',[0.59 0.18 0.37 0.70]); hold(ax2,'on');
  patch(ax2,[-12 0 0 -12],[p.wr p.wr d.negmass(2) d.negmass(2)], ...
    gapGray,'EdgeColor','none','FaceAlpha',1);
  % w is sampled on either side of the pole; a NaN separator stops MATLAB
  % drawing a connector along omega=omega_r. S.w and S.Meff are unchanged.
  ip=find(w<p.wr,1,'last');
  plot(ax2,[me(1:ip) NaN me(ip+1:end)],[w(1:ip) NaN w(ip+1:end)], ...
    'Color',s.blue,'LineStyle','-','LineWidth',1.1);
  xline(ax2,0,'-','Color',guideGray,'LineWidth',0.6);
  yline(ax2,p.wr,':','Color',guideGray,'LineWidth',0.6);
  yline(ax2,d.negmass(2),'--','Color',s.orange,'LineWidth',s.referenceWidth);
  xlabel(ax2,'\itM\rm_{eff}/\itM');
  xlim(ax2,[-12 12]); ylim(ax2,[0 2.3]);
  journal_style('axes',ax2); journal_style('dispersion',ax2); journal_style('panel',ax2,'(b)');
  set([ax1 ax2],'YTick',[0 0.5 1 1.5 2]);
  S.qa=q/pi; S.w_lower=d.w_lower; S.w_upper=d.w_upper; S.w=w; S.Meff=me;
  S.gap=d.gap; S.negmass=d.negmass; S.note='Chain dispersion and the exact negative-effective-mass interval.';
  export_fig_data('fig02_chain_dispersion_effective_mass',S);
  export_figure(fig,'fig02_chain_dispersion_effective_mass'); close(fig);
end
