function [cmu, mu, gam, b1, b2] = beam_cmu_d(f, fr, eta, p)
%BEAM_CMU_D  Damped complex Bloch structure of the Timoshenko beam (Section 6.3).
%
%  The quadratic Eq. (4.5) with COMPLEX beta1, beta2 -- the real() of
%  beam_cmu.m is dropped: with damping beta1, beta2 are genuinely complex (the
%  realness statement of Section 4.2 no longer holds) and cos(mu) is complex.
%
%     cos^2(mu) - beta1 cos(mu) + beta2 = 0                        Eq. (4.5)
%     beta1 = (C1+C2) + i K_T (a1 S1 + a2 S2)                       Eq. (4.6)
%     beta2 =  C1 C2  + i K_T (a1 S1 C2 + a2 S2 C1)                 Eq. (4.6)
%
%  eta = [eta_b, eta_s, eta_r].  Returns, for f of length N:
%     cmu : 2xN complex, the two Bloch branches cos(mu)      (row1 = smaller |Im mu|)
%     mu  : 2xN complex propagation constant, mu = acos(cmu) = delta + i*gamma (Section 6.3)
%     gam : 1xN, min over the two branches of |Im mu| = decay per cell of the
%           LEAST-decaying Bloch wave (Section 6.3)
%     b1,b2 : the complex coefficients beta1, beta2 of Eq. (4.6)
%
%  The lossless propagation test (real cos(mu) in [-1,1], see gap_edges) no
%  longer applies with damping; gamma = |Im mu| replaces it (Section 6.3).  Roots
%  are taken in a cancellation-safe form that reduces to the Vieta extraction of
%  beam_cmu when beta1, beta2 are real.

  if nargin < 4 || isempty(p),   p   = params_beam(); end
  if nargin < 3 || isempty(eta), eta = [0 0 0];       end
  f = f(:).';  N = numel(f);

  [k1, k2]  = beam_wavenumbers_d(f, eta, p);
  [a1, a2]  = beam_coeffs_d(f, fr, eta, p);
  C1 = cos(k1*p.L);  C2 = cos(k2*p.L);
  S1 = sin(k1*p.L);  S2 = sin(k2*p.L);

  w  = 2*pi*f;
  Kc = p.m*(2*pi*fr)^2*(1 + 1i*eta(3));
  KT = Kc*p.m*w.^2./(p.m*w.^2 - Kc);

  b1 = (C1 + C2) + 1i*KT.*(a1.*S1 + a2.*S2);          % complex (NO real())
  b2 = C1.*C2    + 1i*KT.*(a1.*S1.*C2 + a2.*S2.*C1);

  % ---- cancellation-safe complex quadratic  c^2 - b1 c + b2 = 0 -------------
  d   = b1.^2 - 4*b2;
  sq  = sqrt(d);
  sgn = ones(1, N);
  sgn(real(conj(b1).*sq) < 0) = -1;      % add the sqrt in the non-canceling sense
  rbig = (b1 + sgn.*sq)/2;               % larger-magnitude root
  z    = (rbig == 0);
  rsmall = zeros(1, N);
  rsmall(~z) = b2(~z)./rbig(~z);         % product of roots = b2
  rsmall(z)  = (b1(z) - sgn(z).*sq(z))/2;

  % ---- pole at eta_r = 0, f = fr : one branch diverges, other -> Psi, Eq. (4.9)
  if eta(3) == 0
    at = abs(f - fr) <= p.tol_pole*max(fr, 1);
    if any(at)
      num  = a1(at).*S1(at).*C2(at) + a2(at).*S2(at).*C1(at);
      den2 = a1(at).*S1(at) + a2(at).*S2(at);
      rsmall(at) = num./den2;            % Psi (finite branch)
      rbig(at)   = Inf;                  % divergent branch (fully attenuating)
    end
  end

  roots2 = [rsmall; rbig];               % 2xN

  % Exact nodal Bragg state under pure resonator loss (Section 6.3).  When the
  % host is lossless and sin(k1 L)=0, c=C1=+/-1 is an algebraic root for every
  % complex K_T.  Enforce that identity at numerical nodes before acos; near
  % resonance, multiplying roundoff in sin(k1 L) by a large K_T otherwise
  % creates a fictitious attenuation of the pinned state.
  if eta(1) == 0 && eta(2) == 0
    node = abs(S1) < 1e-12;
    jj = find(node);
    for q = 1:numel(jj)
      j = jj(q); [~,ir] = min(abs(roots2(:,j) - C1(j)));
      roots2(ir,j) = C1(j);
    end
  end

  mu2    = acos(roots2);                 % 2xN complex
  mu2(~isfinite(roots2)) = 1i*1e6;       % divergent branch -> gamma = +Inf marker

  % order each column by |Im mu| ascending: row 1 = least-decaying branch
  [~, ord] = sort(abs(imag(mu2)), 1);
  cmu = zeros(2, N);  mu = zeros(2, N);
  lin = ord + repmat((0:N-1)*2, 2, 1);   % column-wise linear indices
  cmu(:) = roots2(lin);
  mu(:)  = mu2(lin);
  gam = abs(imag(mu(1, :)));             % least decay per cell
end
