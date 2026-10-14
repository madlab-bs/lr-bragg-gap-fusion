function fig03_rod_bloch_tunings()
%FIG03_ROD_BLOCH_TUNINGS  Manuscript Fig. 3: rod Bloch dispersion at three resonator tunings.
%  Writes figures/fig03_rod_bloch_tunings.pdf and .png and the plotted data figures/data/fig03_rod_bloch_tunings.dat.
%  Run from matlab/:  fig03_rod_bloch_tunings   (run_all_figures makes Figs. 2-10).
  p=params_rod(); s=journal_style('colors'); fr=[700 p.fB1 1400];
  gapGray=[0.90 0.90 0.90];
  guideGray=[0.55 0.55 0.55];
  f=linspace(20,2200,16000); fig=journal_style('figure',7.0,2.85);
  pos=[0.075 0.17 0.255 0.68; 0.385 0.17 0.255 0.68; 0.695 0.17 0.255 0.68];
  S.f=f; S.fr=fr; S.cmu=NaN(3,numel(f));
  for j=1:3
    ax=subplot(1,3,j); set(ax,'Position',pos(j,:)); hold(ax,'on');
    [~,g]=gap_edges(@(x) rod_cmu(x,fr(j),p),[20 2200], ...
      struct('n',40000,'tol',1e-8,'pole',fr(j)));
    for k=1:size(g,1)
      patch(ax,[0 1 1 0],[g(k,1) g(k,1) g(k,2) g(k,2)]/p.fB1, ...
        gapGray,'EdgeColor','none','FaceAlpha',1);
    end
    c=rod_cmu(f,fr(j),p); S.cmu(j,:)=c;
    % Plot on the full grid with NaN in the stop bands, so separate pass
    % bands are never joined across a gap. The exact gap_edges
    % frequencies are inserted for display only, so each branch
    % ends on mu=0 or mu=1 at the shaded boundary; S.cmu is unchanged.
    fe=[f g(:).']; ce=rod_cmu(fe,fr(j),p);
    ie=numel(f)+1:numel(fe); ce(ie)=max(-1,min(1,ce(ie)));
    [fe,ix]=sort(fe); ce=ce(ix);
    mu=NaN(size(ce)); ok=isfinite(ce)&abs(ce)<=1; mu(ok)=acos(ce(ok))/pi;
    plot(ax,mu,fe/p.fB1,'Color',s.ink,'LineStyle','-', ...
      'LineWidth',1.1);
    mu0=acos(cos(pi*f/p.fB1))/pi;
    plot(ax,mu0,f/p.fB1,'Color',s.gray,'LineStyle','-.', ...
      'LineWidth',0.9);
    yline(ax,1,':','Color',guideGray,'LineWidth',0.7);
    xlabel(ax,'\it\mu\rm/\it\pi');
    if j==1, ylabel(ax,'\it\omega\rm/\it\omega\rm_{B,1}'); end
    xlim(ax,[0 1]); ylim(ax,[0 2.2]);
    text(ax,0.50,1.025,sprintf('\\omega_r/\\omega_{B,1}=%.1f',fr(j)/p.fB1), ...
      'Units','normalized','Clipping','off','HorizontalAlignment','center', ...
      'VerticalAlignment','bottom','FontSize',s.annotationFont);
    journal_style('axes',ax); journal_style('dispersion',ax); journal_style('panel',ax,sprintf('(%c)',96+j));
    if j==1
      % The caption identifies the black solid and gray dash-dot curves;
      % in-axis curve labels would collide with the folded branches here.
      text(ax,0.035,1.015,'\it\omega\rm_{B,1}','Color',guideGray, ...
        'FontSize',s.annotationFont,'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom');
    end
  end
  S.note='Angular-frequency versus Bloch-phase rod dispersion at classical, coincident, and flipped tunings.';
  export_fig_data('fig03_rod_bloch_tunings',S);
  export_figure(fig,'fig03_rod_bloch_tunings'); close(fig);
end
