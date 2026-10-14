function [T, Rf, t1] = beam_transmission(f, fr, eta, N, p)
%BEAM_TRANSMISSION  Power transmission through N damped unit cells, LOSSLESS leads.
%
%  Finite-specimen transmission of Section 6.3 (fig10_finite_transmission).  N cells
%  (each = field transfer T=expm(A L) then resonator jump V(+)=V(-)+K_T v, damped
%  via eta) between semi-infinite LOSSLESS Timoshenko leads.  Incident propagating
%  +k1 wave -> transmitted (+k1,+k2) and reflected (-k1,-k2), the k2 waves being
%  evanescent below the shear cutoff; T=|t1|^2.
%
%  STABLE Bloch-eigenbasis solve: the cell transfer Tc=R*T is
%  diagonalized, Tc = V diag(lam) V^{-1}, and the N-cell propagation is applied as
%  lam^N to the DECAYING Bloch modes (|lam|<=1, referenced to the left face) and
%  lam^{-N} to the GROWING modes (|lam|>1, referenced to the right face).  All
%  coefficients are bounded (no exp(+kappa N L) overflow), so it is accurate for
%  arbitrary N and through the coincidence resonance.  Band-structure gamma is the
%  PRIMARY evidence (Section 6.3); this transmission is a corollary.
%
%  eta = [eta_b,eta_s,eta_r].  f may be a vector; T,Rf,t1 same size.

  if nargin < 5 || isempty(p), p = params_beam(); end
  f = f(:).';  Mf = numel(f);
  T = zeros(1,Mf);  Rf = zeros(1,Mf);  t1 = zeros(1,Mf);

  for i = 1:Mf
    w = 2*pi*f(i);
    [k1,k2,~,~,ks2] = beam_wavenumbers_d(f(i), [0 0 0], p);   % lossless leads
    wR1 = wstate(k1,ks2,p.EI,p.GA);  wL1 = wstate(-k1,ks2,p.EI,p.GA);
    wR2 = wstate(k2,ks2,p.EI,p.GA);  wL2 = wstate(-k2,ks2,p.EI,p.GA);

    Tf  = expm(beam_stateA_d(w, eta, p)*p.L);
    Kc  = p.m*(2*pi*fr)^2*(1 + 1i*eta(3));
    KTv = Kc*p.m*w^2/(p.m*w^2 - Kc);
    Rj  = eye(4); Rj(4,1) = KTv;
    Tc  = Rj*Tf;

    [V,Lam] = eig(Tc);  lam = diag(Lam);  Vi = inv(V);
    % state at left face:  z0 = wR1 + r1 wL1 + r2 wL2   (unknowns r1,r2,t1,t2)
    % state at right face: zN = t1 wR1 + t2 wR2 = Tc^N z0
    C0 = Vi*[wL1, wL2, zeros(4,1), zeros(4,1)];  c0c = Vi*wR1;   % c0 = C0 x + c0c
    CN = Vi*[zeros(4,1), zeros(4,1), wR1, wR2];                   % cN = CN x
    A = zeros(4);  b = zeros(4,1);
    for m = 1:4
      if abs(lam(m)) <= 1
        A(m,:) = CN(m,:) - lam(m)^N * C0(m,:);
        b(m)   = lam(m)^N * c0c(m);
      else
        A(m,:) = C0(m,:) - lam(m)^(-N) * CN(m,:);
        b(m)   = -c0c(m);
      end
    end
    x = A\b;                          % [r1; r2; t1; t2]
    t1(i) = x(3);
    T(i)  = abs(x(3))^2;
    Rf(i) = abs(x(1))^2;
  end
end

function s = wstate(k, ks2, EI, GA)
  s = [ 1; 1i*(k^2-ks2)/k; EI*(k^2-ks2); -1i*GA*ks2/k ];
end
