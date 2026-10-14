function fig08_transfer_matrix_reduction()
%FIG08_TRANSFER_MATRIX_REDUCTION  Manuscript Fig. 8: 4x4 transfer matrix vs the exact 2x2 reduction.
%  Writes figures/fig08_transfer_matrix_reduction.pdf and .png and the plotted data figures/data/fig08_transfer_matrix_reduction.dat.
%  Run from matlab/:  fig08_transfer_matrix_reduction   (run_all_figures makes Figs. 2-10).
  p=params_beam(); s=journal_style('colors'); R=beam_reference_values(p); frstar=R.frstar; fB=R.fB1;  % computed reference values
  f1=linspace(363.85,364.65,601); [ct1,mt1]=beam_bloch_tm(f1,200,[0 0 0],p);
  [cq1,mq1]=beam_cmu_d(f1,200,[0 0 0],p); e1=pair_error(ct1,cq1);
  f2=linspace(fB-0.07,fB+0.09,601); [ct2,mt2]=beam_bloch_tm(f2,frstar,[0 0 0],p);
  [cq2,mq2]=beam_cmu_d(f2,frstar,[0 0 0],p); e2=pair_error(ct2,cq2);
  fig=journal_style('figure',7.0,3.45);
  ax1=subplot(1,2,1); set(ax1,'Position',[0.12 0.18 0.35 0.68]); hold(ax1,'on');
  h1=plot(ax1,break_root_jumps(real(mq1(1,:))/pi,0.35),f1,'Color',s.blue, ...
    'LineStyle','-','LineWidth',1.0);
  h2=plot(ax1,break_root_jumps(real(mq1(2,:))/pi,0.35),f1,'Color',s.orange, ...
    'LineStyle','-.','LineWidth',1.0);
  id=1:25:numel(f1);
  hm=plot(ax1,real(mt1(1,id))/pi,f1(id),'o','Color',s.ink, ...
    'MarkerFaceColor',s.white,'LineWidth',0.8,'MarkerSize',3.2);
  plot(ax1,real(mt1(2,id))/pi,f1(id),'o','Color',s.ink, ...
    'MarkerFaceColor',s.white,'LineWidth',0.8,'MarkerSize',3.2, ...
    'HandleVisibility','off');
  yline(ax1,R.deg200,':','Color',s.gray,'LineWidth',0.9);
  yline(ax1,R.slide200,'--','Color',s.gray,'LineWidth',0.9);
  xlabel(ax1,'Re(\it\mu\rm)/\it\pi'); ylabel(ax1,'\itf\rm (Hz)');
  xlim(ax1,[0 0.035]); ylim(ax1,[f1(1) f1(end)]);
  journal_style('axes',ax1); journal_style('dispersion',ax1); journal_style('panel',ax1,'(a)');

  ax2=subplot(1,2,2); set(ax2,'Position',[0.61 0.18 0.35 0.68]); hold(ax2,'on');
  plot(ax2,break_root_jumps(real(mq2(1,:))/pi,0.35),(f2-fB)*1e3,'Color',s.blue, ...
    'LineStyle','-','LineWidth',1.0);
  plot(ax2,break_root_jumps(real(mq2(2,:))/pi,0.35),(f2-fB)*1e3,'Color',s.orange, ...
    'LineStyle','-.','LineWidth',1.0);
  id=1:24:numel(f2);
  plot(ax2,real(mt2(1,id))/pi,(f2(id)-fB)*1e3,'o','Color',s.ink, ...
    'MarkerFaceColor',s.white,'LineWidth',0.8,'MarkerSize',3.2);
  plot(ax2,real(mt2(2,id))/pi,(f2(id)-fB)*1e3,'o','Color',s.ink, ...
    'MarkerFaceColor',s.white,'LineWidth',0.8,'MarkerSize',3.2);
  yline(ax2,0,':','Color',s.gray,'LineWidth',0.9);
  xlabel(ax2,'Re(\it\mu\rm)/\it\pi'); ylabel(ax2,'\itf\rm-\itf\rm_{B,1} (mHz)');
  xlim(ax2,[0 1]); ylim(ax2,1e3*[f2(1)-fB f2(end)-fB]);
  journal_style('axes',ax2); journal_style('dispersion',ax2);
  journal_style('panel',ax2,'(b)');
  lg=legend(ax1,[h1 h2 hm],{'2\times2 branch 1','2\times2 branch 2', ...
    '4\times4 roots'},'Orientation','horizontal','NumColumns',3, ...
    'FontSize',s.legendFont,'Box','off');
  set(lg,'Units','normalized','Position',[0.18 0.895 0.64 0.045]);
  S.f_degenerate=f1; S.cmu_tm_degenerate=ct1; S.cmu_2x2_degenerate=cq1; S.error_degenerate=e1;
  S.f_exchange=f2; S.cmu_tm_exchange=ct2; S.cmu_2x2_exchange=cq2; S.error_exchange=e2;
  S.guide_degeneracy=R.deg200; S.guide_sliding=R.slide200; S.fr_star=frstar;  % plotted guides
  S.remu_over_pi_2x2_degenerate=real(mq1)/pi; S.remu_over_pi_tm_degenerate=real(mt1)/pi;   % plotted quantities
  S.remu_over_pi_2x2_exchange=real(mq2)/pi;   S.remu_over_pi_tm_exchange=real(mt2)/pi;
  S.note='Independent 4x4 transfer-matrix roots against the exact 2x2 reduction at the degeneracy and mode exchange.';
  export_fig_data('fig08_transfer_matrix_reduction',S);
  export_figure(fig,'fig08_transfer_matrix_reduction'); close(fig);
end
function x=break_root_jumps(x,tol)
% Do not draw artificial connectors when the ordered roots swap rows (the
% same display-only helper as fig05). Stored arrays are unchanged.
  jump=find(abs(diff(x))>tol)+1;
  x(jump)=NaN;
end
function e=pair_error(a,b)
  e=zeros(1,size(a,2));
  for j=1:size(a,2)
    e(j)=min(max(abs(a(:,j)-b(:,j))),max(abs(a(:,j)-b([2 1],j))));
  end
end
