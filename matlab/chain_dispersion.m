function out = chain_dispersion(p, qa)
%CHAIN_DISPERSION  Discrete mass-in-mass chain (Section 2, Eqs. (2.1)-(2.5)).
%
%   out = chain_dispersion(p)      band edges only
%   out = chain_dispersion(p, qa)  also the two branches w_-(qa), w_+(qa)
%
%  Here qa is the manuscript's mu (phase per cell) and Kb its K_b.
%  Dispersion relation, Eq. (2.1):  Meff(w)*w^2 = 2*Kb*(1 - cos(qa))
%  Effective mass, Eq. (2.2):       Meff(w) = M + m*wr^2/(wr^2 - w^2)
%
%  For fixed qa the relation is a quadratic in Lam = w^2 (Section 2.3):
%     M*Lam^2 - [(M+m)*wr^2 + S(qa)]*Lam + S(qa)*wr^2 = 0,   S = 2*Kb*(1-cos qa)

  if nargin < 1 || isempty(p), p = params_chain(); end

  M = p.M; Kb = p.Kb; m = p.m; wr = p.wr;

  % ---- zone-boundary band edges: the two roots of Eq. (2.4) ---------------
  co   = [M, -((M+m)*wr^2 + 4*Kb), 4*Kb*wr^2];
  r    = sort(roots(co));
  out.Delta        = ((M+m)*wr^2 + 4*Kb)^2 - 16*Kb*M*wr^2;   % discriminant of (2.4)
  out.w_acoustic_top = sqrt(r(1));                            % sqrt(Lambda_-)
  out.w_spectrum_top = sqrt(r(2));                            % sqrt(Lambda_+)

  % ---- q = 0 optical mode w_0, Eq. (2.2) and Section 2.2 ------------------
  out.w_optical_bot = wr*sqrt(1 + m/M);
  out.w_pole        = wr;                                     % pole of Meff

  % ---- the single gap, Eq. (2.5) ------------------------------------------
  out.gap        = [out.w_acoustic_top, out.w_optical_bot];
  out.gap_width  = diff(out.gap);
  out.negmass    = [wr, out.w_optical_bot];                   % Eq. (2.3)

  % ---- effective mass handle ----------------------------------------------
  out.Meff = @(w) M + m*wr^2./(wr^2 - w.^2);

  % ---- branches, if requested ---------------------------------------------
  if nargin >= 2 && ~isempty(qa)
    qa = qa(:).';
    S  = 2*Kb*(1 - cos(qa));
    nq = numel(qa);
    out.qa = qa;
    out.w_lower = zeros(1,nq);
    out.w_upper = zeros(1,nq);
    for i = 1:nq
      rr = sort(roots([M, -((M+m)*wr^2 + S(i)), S(i)*wr^2]));
      out.w_lower(i) = sqrt(max(rr(1),0));
      out.w_upper(i) = sqrt(max(rr(2),0));
    end
  end
end
