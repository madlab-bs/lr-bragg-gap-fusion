function [k1, k2, k1sq, k2sq, ks2] = beam_wavenumbers_d(f, eta, p)
%BEAM_WAVENUMBERS_D  Damped Timoshenko wavenumbers (Section 6.3).
%
%  Hysteretic (structural) damping via complex moduli under e^{+iwt} (Section 6.3):
%  EI -> EI(1+i eta_b),  GA -> GA(1+i eta_s).  With complex moduli the groups
%  A,B,C, the discriminant, and both wavenumbers become complex.  The manuscript
%  states only the lossless branch rule Im k <= 0 (Section 4.1); here the same rule
%  is applied to BOTH complex square roots by ksqrt_decay below: take the principal
%  root and flip its sign if Im k > 0, so that each wave e^{-i k |x|} decays away
%  from its source.  At eta = 0 this reduces bit-for-bit to beam_wavenumbers.
%
%     eta = [eta_b, eta_s, eta_r]   (only eta_b, eta_s used here)
%
%  At eta = [0 0 0] this reproduces beam_wavenumbers.m exactly (continuity
%  regression, criterion 9 of validate_all).

  if nargin < 3 || isempty(p),   p   = params_beam(); end
  if nargin < 2 || isempty(eta), eta = [0 0 0];       end

  EIc = p.EI*(1 + 1i*eta(1));
  GAc = p.GA*(1 + 1i*eta(2));
  Ac  = p.rhoI/EIc + p.rhoA/GAc;                 % coefficients of Eq. (4.2), complex
  Bc  = p.rhoI*p.rhoA/(EIc*GAc);
  Cc  = p.rhoA/EIc;

  w    = 2*pi*f;
  disc = Ac.^2.*w.^4 - 4*(Bc.*w.^4 - Cc.*w.^2);  % discriminant of Eq. (4.2) in k^2, complex
  sq   = sqrt(disc);                             % principal branch (Re >= 0)
  k1sq = (Ac.*w.^2 + sq)/2;
  k2sq = (Ac.*w.^2 - sq)/2;
  ks2  = w.^2*p.rhoA/GAc;                         % k_s^2 of Eq. (4.4), complex

  k1 = ksqrt_decay(k1sq);
  k2 = ksqrt_decay(k2sq);
end

function k = ksqrt_decay(ksq)
%  Complex square root with the outgoing/decaying branch Im(k) <= 0.
%  For real ksq this is bit-identical to beam_wavenumbers: ksq>0 -> +sqrt (Im=0, kept);
%  ksq<0 -> principal sqrt is +i*sqrt|.| (Im>0, flipped) -> -i*sqrt|.|.
  k  = sqrt(ksq);
  fl = imag(k) > 0;
  k(fl) = -k(fl);
end
