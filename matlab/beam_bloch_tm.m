function [cmu, mu, gam] = beam_bloch_tm(f, fr, eta, p)
%BEAM_BLOCH_TM  Damped Bloch cos(mu) by the INDEPENDENT 4x4 transfer matrix.
%
%  No Green's function, no a1/a2, no lattice sum: one-cell field transfer
%  T = expm(A L) from beam_stateA_d, point resonator jump V(+) = V(-) + K_T v,
%  cell transfer Tc = R T, and cos(mu) = (lam + 1/lam)/2 from eig(Tc).  This is
%  the damped continuation of validate_all criterion 8, and the independent
%  cross-check of the complex quadratic beam_cmu_d with damping (Section 6.3).
%
%  eta = [eta_b, eta_s, eta_r].  Returns 2xN cos(mu), mu, and 1xN gamma (least
%  |Im mu|).  The four eigenvalues form two reciprocal pairs (lam, 1/lam), one
%  per Bloch branch; each pair gives the same cos(mu), deduplicated below.

  if nargin < 4 || isempty(p),   p   = params_beam(); end
  if nargin < 3 || isempty(eta), eta = [0 0 0];       end
  f = f(:).';  N = numel(f);
  cmu = zeros(2, N);

  for i = 1:N
    w   = 2*pi*f(i);
    T   = expm(beam_stateA_d(w, eta, p)*p.L);
    Kc  = p.m*(2*pi*fr)^2*(1 + 1i*eta(3));
    KTv = Kc*p.m*w^2/(p.m*w^2 - Kc);
    R   = eye(4);  R(4,1) = KTv;              % resonator jump on V (state comp 4)
    lam = eig(R*T);
    c4  = 0.5*(lam + 1./lam);                 % 4 values, equal in reciprocal pairs
    cmu(:, i) = dedupe2(c4);
  end

  mu  = acos(cmu);
  gam = min(abs(imag(mu)), [], 1);
end

function c2 = dedupe2(c4)
%  Collapse 4 cos(mu) values (two equal pairs) to the 2 distinct branch values.
  c4 = c4(:);
  c2 = zeros(2,1);  n = 0;
  while ~isempty(c4)
    c0   = c4(1);
    near = abs(c4 - c0) < 1e-6*max(1, abs(c0));
    n = n + 1;
    if n > 2                                   % more than two clusters: stop rather than average
      error('beam_bloch_tm:clusters', 'more than two distinct cos(mu) clusters: %s', mat2str(c4.', 6));
    end
    c2(n) = mean(c4(near));
    c4(near) = [];
  end
  % At an exact branch collision all four transfer-matrix eigenvalues produce
  % the same cos(mu) and the loop finds one cluster.  The second output would
  % otherwise stay at its initial zero, a spurious root at the discriminant
  % edge; a double Bloch root is therefore returned twice.
  if n == 1, c2(2) = c2(1); end
  % sort by |Im(acos)| so row 1 is the least-decaying branch
  if abs(imag(acos(c2(2)))) < abs(imag(acos(c2(1)))), c2 = c2([2 1]); end
end
