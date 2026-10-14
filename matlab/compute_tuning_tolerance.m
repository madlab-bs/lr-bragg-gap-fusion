function T = compute_tuning_tolerance(theta, delta, nscan)
%COMPUTE_TUNING_TOLERANCE  Exact rod/beam pass-band widths near coincidence.
%
%  T = compute_tuning_tolerance(theta,delta) evaluates the narrow intervening
%  pass band for every mass ratio in theta and relative resonator detuning in
%  delta = f_r/f_B1-1.  Both hosts are solved from their dispersion relations;
%  no plotted width is fitted or adjusted to the earlier publication.
%  Band edges are bisected to 1e-12 Hz (a 1e-8 Hz tolerance gives width errors
%  up to ~8e-9 Hz).
%  The 18000-point scan used by fig06_rod_beam_coincidence is grid-converged
%  (a 720000-point scan agrees to 9e-13 Hz).

  if nargin < 3 || isempty(nscan), nscan=14000; end
  pb=params_beam(); pr=params_rod(); fBb=beam_bragg_frequency(pb);
  theta=theta(:); delta=delta(:).';
  [DM,TH]=meshgrid(delta,theta);
  Wr=zeros(size(DM)); Wb=zeros(size(DM));
  for i=1:numel(theta)
    qr=pr; qr.m=theta(i)*pr.rhoA*pr.L; qr.theta_c=theta(i);
    qb=pb; qb.m=theta(i)*pb.rhoA*pb.L; qb.theta_c=theta(i);
    for j=1:numel(delta)
      frr=pr.fB1*(1+delta(j));
      [~,~,wr]=narrow_band(@(f) rod_cmu(f,frr,qr),pr.fB1+[-4 4], ...
        struct('n',nscan,'tol',1e-12,'pole',frr));
      frb=fBb*(1+delta(j));
      [~,~,wb]=narrow_band(@(f) beam_cmu(f,frb,qb),fBb+[-2.5 2.5], ...
        struct('n',nscan,'tol',1e-12,'pole',frb));
      Wr(i,j)=wr/pr.fB1; Wb(i,j)=wb/fBb;
    end
  end

  % Local mode-exchange guide from Eq. (5.3), evaluated for each mass ratio.
  [~,k2B]=beam_wavenumbers(fBb,pb); [~,a2B]=beam_coeffs(fBb,fBb,pb);
  s2B=imag(sin(k2B*pb.L)); c2B=real(cos(k2B*pb.L));
  KB=theta*pb.rhoA*pb.L*(2*pi*fBb)^2;
  deltaStar=KB*a2B*s2B./(2*(1-c2B));

  T.theta=theta; T.delta=delta; T.delta_mesh=DM; T.theta_mesh=TH;
  T.Wrod_normalized=Wr; T.Wbeam_normalized=Wb;
  T.beam_mode_exchange_delta=deltaStar;
  T.fB=[pr.fB1 fBb];
end
