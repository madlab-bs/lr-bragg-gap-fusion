function A = beam_stateA_d(w, eta, p)
%BEAM_STATEA_D  Damped 4x4 Timoshenko field matrix, state y=[v, phi, M, V].
%  dy/dx = A y for the Timoshenko equations of Section 4.1 (e^{+iwt}), with
%  M = EI phi', V = GA(v' + phi) and the complex moduli EI(1+i eta_b),
%  GA(1+i eta_s) of Section 6.3.  At eta = 0 it equals the lossless state
%  matrix beam_stateA used in validate_all (criterion 8).
  if nargin < 3 || isempty(p),   p   = params_beam(); end
  if nargin < 2 || isempty(eta), eta = [0 0 0];       end
  EIc = p.EI*(1 + 1i*eta(1));
  GAc = p.GA*(1 + 1i*eta(2));
  A = [ 0,            -1,         0,   1/GAc ;
        0,             0,   1/EIc,        0  ;
        0,  -p.rhoI*w^2,         0,        1  ;
       -p.rhoA*w^2,    0,         0,        0 ];
end
