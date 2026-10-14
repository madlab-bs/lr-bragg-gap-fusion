function fig05_beam_complex_bloch()
%FIG05_BEAM_COMPLEX_BLOCH  Manuscript Fig. 5: complex Bloch dispersion of the beam away from coincidence.
%  Writes figures/fig05_beam_complex_bloch.pdf and .png and the plotted data figures/data/fig05_beam_complex_bloch.dat.
%  Run from matlab/:  fig05_beam_complex_bloch   (run_all_figures makes Figs. 2-10).
  p=params_beam(); s=journal_style('colors'); cases=[200 850];
  f=linspace(30,1350,12000); fB=beam_bragg_frequency(p);
  fig=journal_style('figure',7.0,5.35); S.f=f; S.fr=cases;
  pos=[0.10 0.58 0.36 0.30; 0.60 0.58 0.36 0.30; ...
       0.10 0.12 0.36 0.30; 0.60 0.12 0.36 0.30];
  branchStyle={{s.ink,'-'},{s.maroon,'--'}};
  guideGray=[0.48 0.48 0.48];
  for j=1:2
    [~,mu]=beam_cmu_d(f,cases(j),[0 0 0],p);
    S.(['mu' num2str(j) '_re'])=real(mu)/pi;
    S.(['mu' num2str(j) '_im'])=abs(imag(mu));
    ax=subplot(2,2,j); set(ax,'Position',pos(j,:)); hold(ax,'on');
    xr1=break_root_jumps(real(mu(1,:))/pi,0.35);
    xr2=break_root_jumps(real(mu(2,:))/pi,0.35);
    h1=plot(ax,xr1,f/fB,'Color',branchStyle{1}{1}, ...
      'LineStyle',branchStyle{1}{2},'LineWidth',1.0);
    h2=plot(ax,xr2,f/fB,'Color',branchStyle{2}{1}, ...
      'LineStyle',branchStyle{2}{2},'LineWidth',1.0);
    yline(ax,1,':','Color',guideGray,'LineWidth',0.8);
    yline(ax,p.fc/fB,'-.','Color',s.gray,'LineWidth',0.9);
    xlabel(ax,'Re(\it\mu\rm)/\it\pi');
    if j==1, ylabel(ax,'\it\omega\rm/\it\omega\rm_{B,1}'); end
    xlim(ax,[0 1]); ylim(ax,[30 1350]/fB);
    text(ax,0.50,1.045,sprintf('\\it\\omega\\rm_r/\\it\\omega\\rm_{B,1}\\approx%.2f',cases(j)/fB), ...
      'Units','normalized','Clipping','off','HorizontalAlignment','center', ...
      'VerticalAlignment','bottom','FontSize',s.annotationFont);
    journal_style('axes',ax); journal_style('dispersion',ax);
    journal_style('panel',ax,sprintf('(%c)',96+j));
    if j==2, set(ax,'YTickLabel',[]); end
    if j==1, axfirst=ax; L=[h1 h2]; end

    bx=subplot(2,2,2+j); set(bx,'Position',pos(2+j,:)); hold(bx,'on');
    g1=break_root_jumps(abs(imag(mu(1,:))),0.35);
    g2=break_root_jumps(abs(imag(mu(2,:))),0.35);
    plot(bx,g1,f/fB,'Color',branchStyle{1}{1}, ...
      'LineStyle',branchStyle{1}{2},'LineWidth',1.0);
    plot(bx,g2,f/fB,'Color',branchStyle{2}{1}, ...
      'LineStyle',branchStyle{2}{2},'LineWidth',1.0);
    yline(bx,1,':','Color',guideGray,'LineWidth',0.8);
    yline(bx,p.fc/fB,'-.','Color',s.gray,'LineWidth',0.9);
    xlabel(bx,'|Im(\it\mu\rm)|');
    if j==1, ylabel(bx,'\it\omega\rm/\it\omega\rm_{B,1}'); end
    xlim(bx,[0 1.2]); ylim(bx,[30 1350]/fB);
    journal_style('axes',bx); journal_style('dispersion',bx);
    journal_style('panel',bx,sprintf('(%c)',98+j));
    if j==2, set(bx,'YTickLabel',[]); end
  end
  set(findall(fig,'Type','axes'),'YTick',[0.5 1.0 1.5 2.0]);
  lg=legend(axfirst,L,{'\it\mu\rm_1','\it\mu\rm_2'}, ...
    'Orientation','horizontal','NumColumns',2,'FontSize',s.legendFont,'Box','off');
  set(lg,'Units','normalized','Position',[0.39 0.915 0.22 0.040]);
  S.fB=fB; S.fc=p.fc; S.note='Complex Bloch branches for classical and flipped beam tunings away from coincidence.';
  export_fig_data('fig05_beam_complex_bloch',S);
  export_figure(fig,'fig05_beam_complex_bloch'); close(fig);
end

function x=break_root_jumps(x,tol)
% Do not draw artificial connectors when the attenuation-ordered roots swap.
  jump=find(abs(diff(x))>tol)+1;
  x(jump)=NaN;
end
