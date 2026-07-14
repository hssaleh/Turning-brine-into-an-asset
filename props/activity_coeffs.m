function [gamma_pm, phi] = activity_coeffs(m, I_s, par)
%ACTIVITY_COEFFS Mean activity coefficient and osmotic coefficient of the
%NaCl-equivalent brine electrolyte, by Davies (dilute) or single-
%electrolyte Pitzer (concentrated brine, I_s > 0.5) ion-interaction
%theory.
%
% Purpose: Implements Eq. (2.12) Davies and a single-dominant-electrolyte
%          reduction of the multi-ion Pitzer framework, Eqs. (4.1)-(4.3),
%          for the NaCl-equivalent lumped salt that represents brine at
%          high ionic strength (I_s 1.3-2.5 mol/kg per Sec. 4).
% Inputs : m    - NaCl-equivalent molality [mol/kg water]
%          I_s  - ionic strength [mol/kg]
%          par  - parameter struct (par.stream.A_DH, par.pitzer.*)
% Outputs: gamma_pm - mean activity coefficient (-)
%          phi      - osmotic coefficient (-)
% Assumptions: the full multi-ion Pitzer matrix (Eq. 4.1, mixing terms
%          Theta/Psi) is reduced to the dominant Na-Cl binary
%          interaction; this is the documented "Pitzer-lite" simplifi-
%          cation appropriate for NaCl-dominant seawater-brine (Na+/Cl-
%          account for >85% of total ionic strength, Table 1.3). A
%          smooth Davies-to-Pitzer blend over I_s in [0.3, 0.7] avoids a
%          discontinuity at the literature crossover point.
% Reference: Pitzer (1973) J. Phys. Chem. 77, 268; Davies equation;
%          model Sec. 2.4 (Eq. 2.12) and Sec. 4.1 (Eqs. 4.1-4.3).

A_DH = par.stream.A_DH;
sqI  = sqrt(I_s);

% --- Davies (Eq. 2.12), |z+ z-| = 1 for NaCl
log10_gamma_Davies = -A_DH .* (sqI./(1+sqI) - 0.3.*I_s);
gamma_Davies = 10.^log10_gamma_Davies;
phi_Davies   = 0.90 .* ones(size(m));  % dilute-regime fallback (Sec. 2.5 baseline phi(S=35)~0.90); blend weight -> 0 for brine I_s>0.7

% --- Single-electrolyte Pitzer (Na-Cl binary), Eq. (4.1)-(4.3) reduced
A_phi = par.pitzer.A_phi;
b     = par.pitzer.b;
beta0 = par.pitzer.NaCl.beta0;
beta1 = par.pitzer.NaCl.beta1;
Cphi  = par.pitzer.NaCl.Cphi;
alphaP = 2.0;

f_gamma = -A_phi .* ( sqI./(1+b.*sqI) + (2/b).*log(1+b.*sqI) );
expI    = exp(-alphaP.*sqI);
B_gamma = 2.*beta0 + (2.*beta1./(alphaP.^2.*max(I_s,eps))) .* (1 - (1+alphaP.*sqI-0.5.*alphaP.^2.*I_s).*expI);
C_gamma = 1.5.*Cphi;
ln_gamma_Pitzer = f_gamma + m.*B_gamma + m.^2.*C_gamma;
gamma_Pitzer = exp(ln_gamma_Pitzer);

phi_Pitzer = 1 + ( -A_phi.*sqI./(1+b.*sqI) + m.*(beta0 + beta1.*expI) + m.^2.*Cphi );
phi_Pitzer = max(min(phi_Pitzer, 1.5), 0.4);

% --- Smooth blend over the Davies-Pitzer crossover band
w = (I_s - 0.3) ./ (0.7-0.3);
w = max(min(w,1),0);
gamma_pm = (1-w).*gamma_Davies + w.*gamma_Pitzer;
phi      = (1-w).*phi_Davies   + w.*phi_Pitzer;
end
