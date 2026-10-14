function fB = beam_bragg_frequency(p, n)
%BEAM_BRAGG_FREQUENCY  Geometric Timoshenko Bragg frequency k1*L = n*pi, Eq. (4.7).
%
%  The single numerical definition of the resonator-independent beam Bragg
%  frequencies f_B,n (Section 4.3), solved by bisection.  Plotting and damping
%  routines call it instead of using decimal literals, which are not accurate
%  enough at the exactly pinned nodal state when the resonator stiffness is
%  complex (Section 6.3).

  if nargin < 1 || isempty(p), p = params_beam(); end
  if nargin < 2 || isempty(n), n = 1; end
  g = @(f) real(beam_wavenumbers(f, p))*p.L - n*pi;
  lo = max(1e-9, 0.25*n*600);  hi = max(5000, 4*n*1000);
  while g(lo) > 0, lo = 0.5*lo; end
  while g(hi) < 0, hi = 2*hi; end
  for it = 1:100
    mid = 0.5*(lo + hi);
    if g(mid) < 0, lo = mid; else, hi = mid; end
  end
  fB = 0.5*(lo + hi);
end
