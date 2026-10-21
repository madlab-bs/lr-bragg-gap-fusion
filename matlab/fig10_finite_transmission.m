function fig10_finite_transmission()
%FIG10_FINITE_TRANSMISSION  Manuscript Fig. 10: finite-specimen transmission and observability.
%  Writes figures/fig10_finite_transmission.pdf and .png and the plotted data figures/data/fig10_finite_transmission.dat.
%  Run from matlab/:  fig10_finite_transmission   (run_all_figures makes Figs. 2-10).
  p=params_beam(); s=journal_style('colors'); fB=beam_bragg_frequency(p);
  R=beam_reference_values(p); lo=R.coinc(1); hi=R.coinc(2);   % computed residual pass band
  % Lossless resonances are 1e-9 to 1e-8 Hz wide, and near the band edge narrower than double
  % precision resolves: no lossless curve is drawn.  Only damped profiles, on a grid
  % converged to < 0.01 dB (transmission_profile_check.m), are plotted.  The lossless MAXIMUM is reliable (T = 1
  % attains the passivity bound) and is exported below as the source of the text.
  Nspec=100; etas=[1e-5 1e-4 1e-3];
  f=linspace(fB-0.020,fB+0.095,2000); T=zeros(numel(etas),numel(f));
  for j=1:numel(etas), T(j,:)=10*log10(beam_transmission(f,fB,[0 etas(j) 0],Nspec,p)); end
  Ns=[20 50 100 200]; ee=logspace(-6,-2.3,16); P=zeros(numel(Ns),numel(ee));
  ff=linspace(fB+0.10,fB+0.25,24);
  T0max=zeros(1,numel(Ns)); C0=T0max;               % lossless maximum and contrast, per N_c
  for j=1:numel(Ns)
    [~,Tp]=transmission_peaks(fB,[0 0 0],Ns(j),[lo hi],p); T0max(j)=max(Tp);
    C0(j)=10*log10(T0max(j)/median(beam_transmission(ff,fB,[0 0 0],Ns(j),p)));
  end
  for j=1:numel(Ns)
    for k=1:numel(ee)
      [~,Tp]=transmission_peaks(fB,[0 ee(k) 0],Ns(j),[lo hi],p); pk=max(Tp);   % resolved band maximum
      fl=median(beam_transmission(ff,fB,[0 ee(k) 0],Ns(j),p));
      P(j,k)=10*log10(pk/fl);
    end
  end
  fig=journal_style('figure',7.0,3.60);
  ax1=subplot(1,2,1); set(ax1,'Position',[0.125 0.18 0.345 0.65]); hold(ax1,'on');
  patch(ax1,1e3*[lo-fB hi-fB hi-fB lo-fB],[-210 -210 -100 -100], ...
    [0.82 0.82 0.82],'EdgeColor','none','FaceAlpha',0.16);
  cols=[s.blue;s.orange;s.maroon]; styles={'--','-.',':'}; H=gobjects(1,numel(etas));
  for j=1:numel(etas)
    H(j)=plot(ax1,(f-fB)*1e3,T(j,:),'Color',cols(j,:), ...
      'LineStyle',styles{j},'LineWidth',1.0);
  end
  xlabel(ax1,'\itf\rm-\itf\rm_{B,1} (mHz)'); ylabel(ax1,'10 log_{10}\itT\rm (dB)');
  xlim(ax1,1e3*[f(1)-fB f(end)-fB]); ylim(ax1,[-210 -100]);
  lg1=legend(ax1,H,{'\it\eta\rm_s=10^{-5}','\it\eta\rm_s=10^{-4}','\it\eta\rm_s=10^{-3}'}, ...
    'Orientation','horizontal','NumColumns',2,'FontSize',s.legendFont,'Box','off');
  set(lg1,'Units','normalized','Position',[0.09 0.866 0.39 0.080]);
  journal_style('axes',ax1); journal_style('dispersion',ax1);
  journal_style('panel',ax1,'(a)');

  ax2=subplot(1,2,2); set(ax2,'Position',[0.60 0.18 0.36 0.65]); hold(ax2,'on');
  cols2=[s.ink;s.blue;s.orange;s.maroon]; styles2={'-','--','-.',':'};
  marks={'s','o','d','^'}; H2=gobjects(1,numel(Ns));
  for j=1:numel(Ns)
    H2(j)=semilogx(ax2,ee,P(j,:),'Color',cols2(j,:), ...
      'LineStyle',styles2{j},'LineWidth',0.95,'Marker',marks{j}, ...
      'MarkerSize',3.4,'MarkerFaceColor',s.white);
  end
  set(ax2,'XScale','log');
  yline(ax2,3,'--','Color',[0.48 0.48 0.48],'LineWidth',0.8);
  xlabel(ax2,'\it\eta\rm_s'); ylabel(ax2,'\itC_T\rm (dB)');
  xlim(ax2,[ee(1) ee(end)]); ylim(ax2,[-5 max(30,1.05*max(P(:)))]);
  yl=ylim(ax2); ythr=(3-yl(1))/diff(yl);
  text(ax2,0.96,ythr+0.02,'3 dB','Units','normalized', ...
    'Color',[0.42 0.42 0.42],'HorizontalAlignment','right', ...
    'VerticalAlignment','bottom','FontSize',s.annotationFont);
  lg2=legend(ax2,H2,arrayfun(@(n)sprintf('\\itN\\rm_c=%d',n),Ns,'UniformOutput',false), ...
    'Orientation','horizontal','NumColumns',2,'FontSize',s.legendFont,'Box','off');
  set(lg2,'Units','normalized','Position',[0.59 0.875 0.39 0.080]);
  journal_style('axes',ax2); journal_style('dispersion',ax2);
  journal_style('panel',ax2,'(b)');
  S.f=f; S.TdB=T; S.eta_profiles=etas; S.N_cells_profile=Nspec;
  S.eta=ee; S.N_cells=Ns; S.peak_floor=P;
  S.passband=[lo hi]; S.lossless_Tmax=T0max; S.lossless_contrast_dB=C0;
  S.note='Finite-cell transmission. Damped profiles only (lossless resonances unresolvable in double precision); peak_floor and lossless_* use resolved band maxima.';
  export_fig_data('fig10_finite_transmission',S);
  export_figure(fig,'fig10_finite_transmission'); close(fig);
end
