function fig04_rod_gap_map_flip()
%FIG04_ROD_GAP_MAP_FLIP  Manuscript Fig. 4: evolution of the rod stop bands with resonator tuning.
%  Writes figures/fig04_rod_gap_map_flip.pdf and .png and the plotted data figures/data/fig04_rod_gap_map_flip.dat.
%  Run from matlab/:  fig04_rod_gap_map_flip   (run_all_figures makes Figs. 2-10).
  p=params_rod(); s=journal_style('colors'); rr=linspace(0.40,1.60,181);
  gapGray=[0.90 0.90 0.90]; guideGray=[0.55 0.55 0.55];
  fig=journal_style('figure',4.8,3.45);
  ax=axes(fig,'Position',[0.15 0.17 0.80 0.74]); hold(ax,'on');
  lo1=NaN(size(rr)); hi1=lo1; lo2=lo1; hi2=lo1;
  for j=1:numel(rr)
    fr=rr(j)*p.fB1;
    [~,g]=gap_edges(@(f) rod_cmu(f,fr,p),[0.20 2.10]*p.fB1, ...
      struct('n',12000,'tol',2e-6,'pole',fr));
    g=g/p.fB1;
    if abs(rr(j)-1)<10*eps
      % At coincidence the two intervals meet at the fixed Bragg edge.
      lo1(j)=g(1,1); hi1(j)=1; lo2(j)=1; hi2(j)=g(1,2);
    elseif size(g,1)>=2
      lo1(j)=g(1,1); hi1(j)=g(1,2);
      lo2(j)=g(2,1); hi2(j)=g(2,2);
    end
  end
  fill(ax,[rr fliplr(rr)],[lo1 fliplr(hi1)],gapGray, ...
    'EdgeColor','none','FaceAlpha',1);
  fill(ax,[rr fliplr(rr)],[lo2 fliplr(hi2)],gapGray, ...
    'EdgeColor','none','FaceAlpha',1);
  left=rr<=1; right=rr>=1;
  plot(ax,rr,lo1,'Color',s.ink,'LineWidth',0.9);
  plot(ax,rr,hi2,'Color',s.ink,'LineWidth',0.9);
  plot(ax,rr(left),hi1(left),'Color',s.ink,'LineWidth',0.9);
  plot(ax,rr(right),lo2(right),'Color',s.ink,'LineWidth',0.9);
  yline(ax,1,'-','Color',s.ink,'LineWidth',1.0);
  plot(ax,rr,rr,'-.','Color',s.gray,'LineWidth',0.9);
  xline(ax,1,':','Color',guideGray,'LineWidth',0.7);
  xlabel(ax,'\it\omega\rm_r/\it\omega\rm_{B,1}'); ylabel(ax,'\it\omega\rm/\it\omega\rm_{B,1}');
  xlim(ax,[rr(1) rr(end)]); ylim(ax,[0.2 2.05]);
  set(ax,'XTick',[0.4 0.7 1.0 1.3 1.6],'YTick',[0.5 1.0 1.5 2.0]);
  journal_style('axes',ax); journal_style('dispersion',ax);
  S.ratio=rr; S.lower_gap=[lo1;hi1]; S.upper_gap=[lo2;hi2];
  S.note='Continuous rod stop-band evolution through exact LR-Bragg coincidence.';
  export_fig_data('fig04_rod_gap_map_flip',S);
  export_figure(fig,'fig04_rod_gap_map_flip'); close(fig);
end
