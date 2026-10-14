function nfail = validate_figure_program()
%VALIDATE_FIGURE_PROGRAM  Numerical checks on the quantities plotted in the manuscript figures.
%  NFAIL = VALIDATE_FIGURE_PROGRAM() runs 10 checks and returns the number that fail.
%  Run from matlab/:  validate_figure_program   (expected summary: 10 passed, 0 failed).
  pb=params_beam(); pr=params_rod(); pc=params_chain(); fB=beam_bragg_frequency(pb);
  nfail=0; npass=0;
  fprintf('\nJSV figure-program validation\n');
  [npass,nfail]=check(npass,nfail,'Reference systems are physically admissible', ...
    pc.M>0&&pc.m>0&&pr.EA>0&&pb.EI>0&&pb.GA>0);
  ch=chain_dispersion(pc);
  % Chain edges (Section 2.2, Eqs. (2.2), (2.5)) against the independent 2x2 generalized eigenproblem
  % of the cell at mu = pi and mu = 0, not against quoted reference values.
  kc=pc.m*pc.wr^2; Mm=diag([pc.M pc.m]);
  wpi=sqrt(sort(real(eig([4*pc.Kb+kc -kc;-kc kc],Mm)))); w0=sqrt(max(sort(real(eig([kc -kc;-kc kc],Mm))),0));
  [npass,nfail]=check(npass,nfail,'Chain edges/effective-mass zero', ...
    abs(ch.w_acoustic_top-wpi(1))<1e-12 && abs(ch.w_optical_bot-w0(2))<1e-12);
  % Rod states (Section 3): at f_B1, cos(mu) = -1 for the classical (f_r = 700) and flipped
  % (f_r = 1400) tunings, Eq. (3.3).  The coincident value -(1+Gamma_1) of Eq. (3.8)
  % (Gamma_1 = pr.gamma1) is tested as a LIMIT: the substituted value at f = f_r = f_B1 is
  % algebraically -(1+Gamma_1), so comparing it would be circular.  The unsubstituted relation at
  % f_B1(1 +- d) must approach -(1+Gamma_1) linearly in d from both sides.
  rl=@(d) rod_cmu(pr.fB1*(1+d),pr.fB1,pr)+(1+pr.gamma1);
  % d = 1e-6, 1e-7: for smaller d, round-off in K_T ~ 1/d dominates the linear approach.
  e6=max(abs([rl(1e-6) rl(-1e-6)])); e7=max(abs([rl(1e-7) rl(-1e-7)]));
  [npass,nfail]=check(npass,nfail,'Rod classical/coincident/flipped states', ...
    abs(rod_cmu(pr.fB1,700,pr)+1)<1e-12 && abs(rod_cmu(pr.fB1,1400,pr)+1)<1e-12 && ...
    e6<1e-4 && e7/e6>0.05 && e7/e6<0.2);
  % Flip criterion, Eq. (3.6): sign of K_T(f_B1), Eq. (3.1).
  KT=@(fr) pr.m*(2*pi*fr)^2*pr.m*(2*pi*pr.fB1)^2/(pr.m*(2*pi*pr.fB1)^2-pr.m*(2*pi*fr)^2);
  [npass,nfail]=check(npass,nfail,'Exact flip sign',KT(700)>0&&KT(1400)<0);
  fr5=pr.fB1+1;
  [~,~,W5]=narrow_band(@(f) rod_cmu(f,fr5,pr),pr.fB1+[-4 4], ...
    struct('n',30000,'tol',1e-9,'pole',fr5));
  W5_asym=2*abs(fr5-pr.fB1)/(pr.gamma1+2);   % small-detuning rod width, Eq. (3.9)
  [~,mu6]=beam_cmu_d([199 500 850 1100],200,[0 0 0],pb);   % f = f_r (K_T pole) excluded
  [npass,nfail]=check(npass,nfail,'Complex Bloch roots finite away from poles',all(isfinite(mu6(:))));
  % Widths (Section 4.4, Eq. (4.11); Fig. 6(c)) by three routes: closed forms and edge roots
  % (beam_reference_values), independent bisection (narrow_band), and a 40-digit evaluation
  % (beam_vpa_reference).  The rod width must vanish at coincidence; the beam width must be
  % smaller at the mode-exchange tuning f_r* than at coincidence.
  [R,V,refok,refmsg]=safe_refs(pb);  % fail soft: a broken reference fails its checks, not the run
  [~,~,Wb]=narrow_band(@(f) beam_cmu(f,fB,pb),fB+[-0.3 0.3],struct('n',60000,'tol',1e-12,'pole',fB));
  [~,~,Wr]=narrow_band(@(f) rod_cmu(f,pr.fB1,pr),pr.fB1+[-8 8],struct('n',30000,'tol',1e-8,'pole',pr.fB1));
  frstar=R.frstar;  Wbstar=NaN;
  if refok
    [~,~,Wbstar]=narrow_band(@(f) beam_cmu(f,frstar,pb),fB+[-0.3 0.3], ...
      struct('n',60000,'tol',1e-12,'pole',frstar));
  end
  [npass,nfail]=check(npass,nfail,'Exact width-versus-detuning comparison', ...
    refok && Wr==0 && Wbstar<Wb && abs(W5-W5_asym)/W5_asym<5e-3 && ...
    max(abs([Wb-V.W, R.coinc(3)-V.W, Wbstar-V.Wstar, R.mexch(3)-V.Wstar]))<1e-10, refmsg);
  % f_r*, the mode-exchange tuning (Section 5.2; Eq. (5.3) is its first-order estimate), against
  % its 40-digit value and against the unit-cell crossing it defines (config (b) reaches f_B1):
  % d(config_b)/d(f_r) = 0.20, so 5*tol_bisect resolves an f_r* error of ~2.5e-8 Hz.
  okx=false;
  if refok
    u=unitcell_modes('beam',frstar,fB+[-0.3 0.3],pb);
    okx=~isempty(u.config_b) && abs(u.config_b(1)-fB)<5*pb.tol_bisect && abs(frstar-V.frstar)<1e-9;
  end
  [npass,nfail]=check(npass,nfail,'Mode exchange at f_r*',okx,refmsg);
  fdg=R.deg200; if ~isfinite(fdg), fdg=364.2559945; end   % fallback test point only (degeneracy edge, Section 5.4)
  ft=[364.0 fdg 364.5 fB]; ct=beam_bloch_tm(ft,200,[0 0 0],pb); cq=beam_cmu_d(ft,200,[0 0 0],pb);
  err=pair_error_local(ct,cq);
  [npass,nfail]=check(npass,nfail,'4x4 roots equal exact 2x2 roots', ...
    max(err([1 3 4]))<1e-9 && err(2)<1e-6);
  % Pure resonator loss leaves gamma(f_B1) = 0 (Section 6.3).  beam_cmu_d enforces this nodal root
  % in double precision, so that value cannot fail; the unsnapped variable-precision value must
  % fall with the working precision (exact zero), not saturate.  Host loss must give gamma > 0.
  g32=beam_gamma_node_vpa(1e-3,pb,32); g50=beam_gamma_node_vpa(1e-3,pb,50);
  [~,~,gs]=beam_cmu_d(fB,fB,[0 1e-4 0],pb);
  [npass,nfail]=check(npass,nfail,'Resonator-loss pinning and host-loss attenuation', ...
    g32<1e-12 && g50<1e-6*g32 && gs>1e-3);
  [T0,R0]=beam_transmission(500,200,[0 0 0],12,pb);
  [npass,nfail]=check(npass,nfail,'Lossless finite specimen conserves energy',abs(T0+R0-1)<1e-5);
  fprintf('  %d passed, %d failed\n',npass,nfail);
end
function [np,nf]=check(np,nf,name,ok,msg)
  if ok, np=np+1; fprintf('  [PASS] %s\n',name);
  else, nf=nf+1; fprintf('  [FAIL] %s\n',name); if nargin>4 && ~isempty(msg), fprintf('         %s\n',msg); end
  end
end
function [R,V,ok,msg]=safe_refs(pb)
  ok=true; msg='';
  try
    R=beam_reference_values(pb); V=beam_vpa_reference(R,pb,40);
  catch e
    ok=false; msg=['reference computation failed: ' e.message];
    R=struct('fB1',NaN,'frstar',NaN,'coinc',NaN(1,3),'mexch',NaN(1,3),'deg200',NaN,'slide200',NaN);
    V=struct('fB1',NaN,'frstar',NaN,'W',NaN,'Wstar',NaN,'psiB_err',NaN);
  end
end
function e=pair_error_local(a,b)
  e=zeros(1,size(a,2));
  for j=1:size(a,2), e(j)=min(max(abs(a(:,j)-b(:,j))),max(abs(a(:,j)-b([2 1],j)))); end
end
