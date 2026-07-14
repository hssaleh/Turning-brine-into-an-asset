function out = red_power(stream, par)
%RED_POWER Reverse-electrodialysis (RED) energy-recovery pathway.
%
% Purpose: Implements Eqs. (3.1)-(3.11) of the brine-valorization model:
%          cell-pair EMF, stack resistance, matched-load gross power,
%          pumping parasitic, net power and the Gibbs-mixing thermo-
%          dynamic ceiling / second-law efficiency.
% Inputs : stream - struct from M00_streams.m (uses c_NaCl_b/_f,
%             gamma_b/_f, T_f, Q_b)
%          par    - parameter struct (par.red.*, par.const.*)
% Outputs: out (struct), N x 1 vectors:
%             E_cp [V], N_cp [-], R_int [Ohm], P_RED_gross [W], P_d [W/m2],
%             P_pump [W], P_net [W], dG_mix_rate [W], eta_II [-]
% Assumptions: dilute (feed) and draw (brine) channel flows are matched
%          1:1 to the brine flow Q_b (lab/industrial RED stack practice);
%          channel-sizing flow per cell pair vdot_cp and membrane area
%          per cell pair A_m are fixed design constants (par.red).
% Reference: model Sec. 3.1, Eqs. (3.1)-(3.11).

R = par.const.R; F = par.const.F;
T = stream.T_f;
% Draw = brine; dilute = low-salinity treated-effluent counter-stream
% (par.sge.S_dilute), not full seawater, which sets the realisable gradient.
[c_dil, gam_dil] = dilute_state(par);
c_D = stream.c_NaCl_b; c_F = c_dil .* ones(size(c_D));
gamma_D = stream.gamma_b; gamma_F = gam_dil .* ones(size(c_D));
Q_b = stream.Q_b;
A_m_per_cp = 1.0;  % m^2 per cell pair (design assumption)

a_D = gamma_D.*c_D; a_F = gamma_F.*c_F;
a_F = max(a_F, 1e-6);

%% Eq. 3.2 cell-pair EMF
E_cp = 2*par.red.alpha .* (R.*T/F) .* log(a_D./a_F);
E_cp = max(E_cp, 0);

%% Eq. 3.10 sizing: number of cell pairs from brine flow / channel-sizing flow
N_cp = Q_b ./ par.red.vdot_cp;

%% Eq. 3.4 area resistance per cell pair
kappa_D = par.red.Lambda_NaCl .* c_D;   % S/m
kappa_F = par.red.Lambda_NaCl .* c_F;   % S/m
kappa_F = max(kappa_F, 1e-6);
r_cp = par.red.R_CEM + par.red.R_AEM + par.red.d_D./(kappa_D.*par.red.epsilon_sp) ...
       + par.red.d_F./(kappa_F.*par.red.epsilon_sp);

%% Eq. 3.3, 3.5, 3.6 stack EMF, internal resistance, gross power
E_stack = N_cp .* E_cp;
R_int   = N_cp .* r_cp ./ A_m_per_cp;
P_RED_gross = E_stack.^2 ./ (4.*R_int);

%% Eq. 3.7 areal power density (benchmark, 1-4 W/m^2 literature)
P_d = E_cp.^2 ./ (4.*r_cp);

%% Eq. 3.8 pumping parasitic (matched flow both compartments)
mu = 1.0e-3 .* ones(size(Q_b)); % representative brine/seawater dynamic viscosity, Pa s (engineering default)
dP = 12.*mu.*par.red.L_path.*par.red.v ./ (par.red.d_h.^2 .* par.red.epsilon_sp);
P_pump = (Q_b.*dP + Q_b.*dP) ./ par.red.eta_pump;

%% Eq. 3.11 thermodynamic ceiling (Gibbs free energy of mixing rate, mole-fraction proxy)
Vw_molar = par.const.M_w/par.const.rho_w;
x_D = min(max(c_D.*Vw_molar,eps), 0.5); x_F = min(max(c_F.*Vw_molar,eps), 0.5);
Q_F = par.sge.avail_dilute .* Q_b; Q_M = Q_b + Q_F; c_M = (Q_b.*c_D + Q_F.*c_F)./Q_M;
x_M = min(max(c_M.*Vw_molar,eps), 0.5);
% Available mixing work rate = -dG_mix (positive, spontaneous mixing releases work)
dG_mix_rate = R.*T.*( Q_b.*c_D.*log(x_D./x_M) + Q_F.*c_F.*log(x_F./x_M) );
dG_mix_rate = max(dG_mix_rate, eps);

%% Eq. 3.9 net power, capped at the thermodynamically admissible fraction of
%% the Gibbs-mixing ceiling (a finite stack cannot exceed eta_II_max*dG_mix)
P_gross_capped = min(P_RED_gross, par.red.eta_II_max .* dG_mix_rate);
P_net = par.red.eta_DC .* P_gross_capped - P_pump;
P_net = max(P_net, 0);

eta_II = min(P_net ./ dG_mix_rate, par.red.eta_II_max);  % Eq. 3.11, second-law efficiency

%% Economically-relevant membrane area: only enough stack to extract the
%% gradient-limited (capped) power at the achievable areal density P_d. The
%% channel-flow count N_cp over-sizes the membrane several-fold because the
%% Gibbs-mixing ceiling, not membrane area, limits the deliverable power.
A_econ = min(N_cp, P_gross_capped ./ max(P_d, eps));

out.E_cp = E_cp; out.N_cp = N_cp; out.A_econ = A_econ; out.R_int = R_int;
out.P_RED_gross = P_RED_gross; out.P_gross_capped = P_gross_capped;
out.P_d = P_d; out.P_pump = P_pump;
out.P_net = P_net; out.dG_mix_rate = dG_mix_rate; out.eta_II = eta_II;
end

% ----------------------------------------------------------------------
function [c_dil, gam_dil] = dilute_state(par)
%DILUTE_STATE Molar concentration (mol m^-3) and activity coefficient of the
%low-salinity treated-effluent counter-stream used by RED and PRO.
S = par.sge.S_dilute;                               % g/kg
m = (S./(1000-S))./par.const.M_NaCl;                % mol/kg water
rho = 1000;                                          % near pure water
c_dil = m .* rho ./ (1 + S/1000);                    % mol/m^3
I = m;                                               % 1:1 electrolyte
gam_dil = 10.^(-par.stream.A_DH .* (sqrt(I)./(1+sqrt(I)) - 0.3.*I));
end
