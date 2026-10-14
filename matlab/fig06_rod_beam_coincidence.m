function fig06_rod_beam_coincidence()
%FIG06_ROD_BEAM_COINCIDENCE  Manuscript Fig. 6: rod-beam comparison near coincidence; width vs detuning.
%  Writes figures/fig06_rod_beam_coincidence.pdf and .png and the plotted data figures/data/fig06_rod_beam_coincidence.dat.
%  Run from matlab/:  fig06_rod_beam_coincidence   (run_all_figures makes Figs. 2-10).
  pr=params_rod(); pb=params_beam(); s=journal_style('colors');
  fBb=beam_bragg_frequency(pb); gapGray=[0.90 0.90 0.90];
  guideGray=[0.55 0.55 0.55]; fig=journal_style('figure',7.0,5.65);

  % Rod: the coincidence gap is finite and contains no real Bloch branch.
  fr=linspace(250,1900,22000);
  cr=rod_cmu(fr,pr.fB1,pr); ok=isfinite(cr)&abs(cr)<=1;
  [~,gr]=gap_edges(@(f) rod_cmu(f,pr.fB1,pr),[250 1900], ...
    struct('n',60000,'tol',1e-8,'pole',pr.fB1));
  ax1=subplot(2,2,1); set(ax1,'Position',[0.10 0.59 0.36 0.29]); hold(ax1,'on');
  for k=1:size(gr,1)
    patch(ax1,[0 1 1 0],([gr(k,1) gr(k,1) gr(k,2) gr(k,2)]/pr.fB1-1), ...
      gapGray,'EdgeColor','none','FaceAlpha',1);
  end
  mur=NaN(size(fr)); mur(ok)=acos(cr(ok))/pi;
  plot(ax1,mur,fr/pr.fB1-1,'Color',s.ink,'LineStyle','-','LineWidth',1.0);
  yline(ax1,0,':','Color',guideGray,'LineWidth',0.8);
  xlabel(ax1,'\it\mu\rm/\it\pi'); ylabel(ax1,'(\it\omega\rm-\it\omega\rm_{B,1})/\it\omega\rm_{B,1}');
  xlim(ax1,[0 1]); ylim(ax1,[250 1900]/pr.fB1-1);
  journal_style('axes',ax1); journal_style('dispersion',ax1);
  journal_style('panel',ax1,'(a)');

  % Beam: gray stop bands leave the narrow residual pass band white.
  fb=linspace(fBb-0.08,fBb+0.10,18000);
  [~,mub]=beam_cmu_d(fb,fBb,[0 0 0],pb);
  [flo,fhi,Wcoin]=narrow_band(@(f) beam_cmu(f,fBb,pb), ...
    fBb+[-0.08 0.10],struct('n',60000,'tol',1e-9,'pole',fBb));
  yb=1e6*(fb/fBb-1); ylo=1e6*(flo/fBb-1); yhi=1e6*(fhi/fBb-1);
  ax2=subplot(2,2,2); set(ax2,'Position',[0.60 0.59 0.36 0.29]); hold(ax2,'on');
  patch(ax2,[0 1 1 0],[-130 -130 ylo ylo],gapGray,'EdgeColor','none','FaceAlpha',1);
  patch(ax2,[0 1 1 0],[yhi yhi 160 160],gapGray,'EdgeColor','none','FaceAlpha',1);
  x1=break_root_jumps(real(mub(1,:))/pi,0.35);
  x2=break_root_jumps(real(mub(2,:))/pi,0.35);
  plot(ax2,x1,yb,'Color',s.ink,'LineStyle','-','LineWidth',1.0);
  plot(ax2,x2,yb,'Color',s.maroon,'LineStyle','--','LineWidth',1.0);
  yline(ax2,0,':','Color',guideGray,'LineWidth',0.8);
  xlabel(ax2,'Re(\it\mu\rm)/\it\pi');
  ylabel(ax2,'10^6(\it\omega\rm-\it\omega\rm_{B,1})/\it\omega\rm_{B,1}');
  xlim(ax2,[0 1]); ylim(ax2,[-130 160]);
  journal_style('axes',ax2); journal_style('dispersion',ax2);
  journal_style('panel',ax2,'(b)');

  % Exact geometry comparison at the reference mass ratio.
  R=beam_reference_values(pb); frstar=R.frstar;   % computed mode-exchange tuning, Eq. (5.2)
  delta=sort(unique([linspace(-0.0012,0.0012,81) frstar/fBb-1 0]));
  T=compute_tuning_tolerance(pb.theta_c,delta,18000);
  WrRel=T.Wrod_normalized; WbRel=T.Wbeam_normalized;
  [~,ic]=min(abs(delta)); [~,is]=min(abs(fBb*(1+delta)-frstar));
  ax3=subplot(2,2,[3 4]); set(ax3,'Position',[0.13 0.115 0.82 0.31]); hold(ax3,'on');
  hr=plot(ax3,1e3*delta,1e6*WrRel,'Color',s.ink, ...
    'LineStyle','-','LineWidth',1.0);
  hb=plot(ax3,1e3*delta,1e6*WbRel,'Color',s.navy, ...
    'LineStyle','--','LineWidth',1.0);
  xline(ax3,0,':','Color',guideGray,'LineWidth',0.8);
  plot(ax3,0,1e6*WrRel(ic),'o','Color',s.ink, ...
    'MarkerFaceColor',s.white,'LineWidth',1.0,'MarkerSize',5);
  plot(ax3,0,1e6*WbRel(ic),'o','Color',s.navy, ...
    'MarkerFaceColor',s.white,'LineWidth',1.0,'MarkerSize',5);
  plot(ax3,1e3*delta(is),1e6*WbRel(is),'d','Color',s.navy, ...
    'MarkerFaceColor',s.white,'LineWidth',1.0,'MarkerSize',6);
  xlabel(ax3,'10^3(\it\omega\rm_r/\it\omega\rm_{B,1}-1)');
  ylabel(ax3,'10^6\itW_n');
  xlim(ax3,1e3*[delta(1) delta(end)]); ylim(ax3,[0 330]);
  legend(ax3,[hr hb],{'rod','Timoshenko beam'},'Location','north', ...
    'Orientation','horizontal','NumColumns',2,'FontSize',s.legendFont,'Box','off');
  journal_style('axes',ax3); journal_style('panel',ax3,'(c)');
  set(ax3,'Position',[0.13 0.115 0.82 0.31]);

  S.f_rod=fr; S.cmu_rod=cr; S.f_beam=fb; S.mu_beam=real(mub)/pi;
  S.beam_passband=[flo fhi Wcoin]; S.delta=delta; S.theta_c=pb.theta_c;
  S.Wrod_normalized=WrRel; S.Wbeam_normalized=WbRel;
  S.fr_star=frstar; S.Wbeam_at_coincidence=WbRel(ic);
  S.Wbeam_minimum=WbRel(is); S.fB=[pr.fB1 fBb];
  S.note='Direct rod-beam coincidence comparison and normalized pass-band width versus detuning.';
  export_fig_data('fig06_rod_beam_coincidence',S);
  export_figure(fig,'fig06_rod_beam_coincidence'); close(fig);
end

function x=break_root_jumps(x,tol)
% Do not draw artificial connectors when the attenuation-ordered roots swap.
  jump=find(abs(diff(x))>tol)+1;
  x(jump)=NaN;
end
