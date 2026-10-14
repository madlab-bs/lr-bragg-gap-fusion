function fig09_damping_attenuation()
%FIG09_DAMPING_ATTENUATION  Manuscript Fig. 9: damping and attenuation at coincidence.
%  Writes figures/fig09_damping_attenuation.pdf and .png and the plotted data figures/data/fig09_damping_attenuation.dat.
%  Run from matlab/:  fig09_damping_attenuation   (run_all_figures makes Figs. 2-10).
  p=params_beam(); s=journal_style('colors'); fB=beam_bragg_frequency(p);
  f=sort(unique([linspace(fB-0.035,fB+0.085,900) fB])); etaS=[0 1e-5 1e-4]; etaR=1e-3;
  G=zeros(4,numel(f));
  for j=1:3, [~,~,G(j,:)]=beam_cmu_d(f,fB,[0 etaS(j) 0],p); end
  [~,~,G(4,:)]=beam_cmu_d(f,fB,[0 0 etaR],p);
  ee=logspace(-6,-1,51); M=damping_metrics(fB,p,ee); gs=M.gmin_s;
  % Resonator-only curve: gamma(f_B1) computed in 100-digit vpa without enforcing the
  % nodal root (no snap), so the plotted values are computed, not set to the exact zero
  % of Section 6.3 (as M.gmin_r is); they are the vpa round-off floor, below 1e-50.
  gr=beam_gamma_node_vpa(ee,p,100);
  [~,k2]=beam_wavenumbers(fB,p); floorv=abs(imag(k2))*p.L;
  fig=journal_style('figure',7.0,3.60);
  ax1=subplot(1,2,1); set(ax1,'Position',[0.10 0.18 0.36 0.65]); hold(ax1,'on');
  h1=plot(ax1,(f-fB)*1e3,G(1,:),'Color',s.ink, ...
    'LineStyle','-','LineWidth',1.0);
  h2=plot(ax1,(f-fB)*1e3,G(2,:),'Color',s.blue, ...
    'LineStyle','--','LineWidth',1.0);
  h3=plot(ax1,(f-fB)*1e3,G(3,:),'Color',s.orange, ...
    'LineStyle','-.','LineWidth',1.0);
  h4=plot(ax1,(f-fB)*1e3,G(4,:),'Color',s.ink, ...
    'LineStyle','--','LineWidth',0.9,'Marker','o', ...
    'MarkerIndices',1:90:numel(f),'MarkerSize',3.2,'MarkerFaceColor',s.white);
  xline(ax1,0,':','Color',[0.58 0.58 0.58],'LineWidth',0.8);
  yline(ax1,floorv,':','Color',s.gray,'LineWidth',0.9);
  % Mark the exact nodal state gamma(f_B1)=0 of the resonator-only curve (a single grid point,
  % otherwise hidden by the f_B1 guide line).
  iB=find(f==fB,1);
  plot(ax1,0,G(4,iB),'o','Color',s.ink,'MarkerFaceColor',s.ink,'MarkerSize',3.6,'LineWidth',0.9);
  xlabel(ax1,'\itf\rm-\itf\rm_{B,1} (mHz)'); ylabel(ax1,'\it\gamma');
  xlim(ax1,1e3*[f(1)-fB f(end)-fB]); ylim(ax1,[0 0.09]);
  lg1=legend(ax1,[h1 h2 h3 h4],{'lossless','\it\eta\rm_s=10^{-5}', ...
    '\it\eta\rm_s=10^{-4}','\it\eta\rm_r=10^{-3}'},'Orientation','horizontal', ...
    'NumColumns',2,'FontSize',s.legendFont,'Box','off');
  set(lg1,'Units','normalized','Position',[0.09 0.866 0.39 0.080]);
  journal_style('axes',ax1); journal_style('dispersion',ax1);
  journal_style('panel',ax1,'(a)');

  ax2=subplot(1,2,2); set(ax2,'Position',[0.60 0.18 0.36 0.65]); hold(ax2,'on');
  hp=semilogx(ax2,ee,gs,'Color',s.blue,'LineStyle','-', ...
    'LineWidth',1.0);
  hr=semilogx(ax2,ee,gr,'Color',s.ink,'LineStyle','--', ...
    'LineWidth',0.9,'Marker','o','MarkerIndices',1:5:numel(ee), ...
    'MarkerSize',3.2,'MarkerFaceColor',s.white);
  set(ax2,'XScale','log');
  yline(ax2,floorv,':','Color',s.gray,'LineWidth',0.9);
  text(ax2,0.95,0.80,'\it\kappa L','Units','normalized','Color',s.gray, ...
    'HorizontalAlignment','right','VerticalAlignment','bottom', ...
    'FontSize',s.annotationFont);
  xlabel(ax2,'\it\eta'); ylabel(ax2,'\it\gamma\rm_{min}');
  xlim(ax2,[ee(1) ee(end)]); ylim(ax2,[0 1.25*floorv]);
  set(ax2,'XTick',[1e-6 1e-4 1e-2]);
  lg2=legend(ax2,[hp hr],{'host shear','resonator only'}, ...
    'Orientation','horizontal','NumColumns',2,'FontSize',s.legendFont,'Box','off');
  set(lg2,'Units','normalized','Position',[0.585 0.895 0.31 0.045]);
  journal_style('axes',ax2); journal_style('dispersion',ax2);
  journal_style('panel',ax2,'(b)');
  S.f=f; S.gamma_profiles=G; S.eta_s=etaS; S.eta_r=etaR; S.eta=ee; S.gamma_host=gs; S.gamma_resonator=gr;
  S.kappaL=floorv; S.note='Host loss fills the residual-pass-band dip; pure resonator loss preserves the exact Bragg node.';
  export_fig_data('fig09_damping_attenuation',S);
  export_figure(fig,'fig09_damping_attenuation'); close(fig);
end
