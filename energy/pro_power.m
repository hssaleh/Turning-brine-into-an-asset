function out = pro_power(stream, par)
%PRO_POWER Pressure-retarded-osmosis (PRO) energy-recovery pathway.
%
% Purpose: Implements Eqs. (3.12)-(3.18): implicit water-flux solution
%          with internal/external concentration polarisation, salt
%          reverse flux, areal/maximum power density, gross/net power.
% Inputs : stream - struct from M00_streams.m (pi_b, pi_f, c_NaCl_b/_f, Q_b)
%          par    - parameter struct (par.pro.*, par.red.eta_pump/eta_DC)
% Outputs: out (struct), N x 1 vectors:
%             J_w [m/s], J_s [mol/m2/s], W_d [W/m2], W_d_max [W/m2],
%             A_PRO [m2], P_gross [W], P_pump [W], P_net [W]
% Assumptions: membrane area sized proportionally to brine flow
%          (par.pro.A_PRO_per_Qb); the bulk osmotic pressures pi_b, pi_f
%          from M00_streams.m serve as pi_D, pi_F in the flux equation
%          (Eq. 3.12), with ICP/ECP captured through the exponential
%          polarisation factors only. Implicit equation solved by fixed-
%          point (damped) iteration, vectorised over plants.
% Reference: model Sec. 3.2, Eqs. (3.12)-(3.18).

% Draw = pressurised brine; dilute feed = low-salinity treated effluent
% (par.sge.S_dilute) -> larger effective osmotic driving force.
pi_D = stream.pi_b;
m_dil = (par.sge.S_dilute./(1000-par.sge.S_dilute))./par.const.M_NaCl;   % mol/kg
pi_F = 0.90 .* 2 .* m_dil .* par.const.R .* stream.T_f .* par.const.rho_w; % Eq.2.13, dilute (phi=0.9, nu=2)
Q_b = stream.Q_b;
A_w = par.pro.A_w; B_s = par.pro.B_s;
K_ICP = par.pro.S_struct ./ par.pro.D_s;
k_d   = par.pro.k_d;

dpi_eff = max(pi_D - pi_F, eps);
dP = dpi_eff/2;                                    % par.pro default operating point

%% Eq. 3.12 implicit water flux, damped fixed-point iteration
J_w = A_w .* dpi_eff * 0.3;                        % initial guess
for it = 1:60
    J_w_new = A_w .* ( pi_D.*exp(-J_w./k_d) - pi_F.*exp(J_w.*K_ICP) ) - A_w.*dP;
    J_w_new = max(J_w_new, 0);
    J_w = 0.5*J_w + 0.5*J_w_new;                   % damping for stability
end

%% Eq. 3.13 salt reverse flux
J_s = B_s .* ( stream.c_NaCl_b.*exp(-J_w./k_d) - stream.c_NaCl_f.*exp(J_w.*K_ICP) );

%% Eq. 3.14-3.15 areal / maximum power density
W_d     = J_w .* dP;
W_d_max = A_w .* dpi_eff.^2 / 4;

%% Membrane area sizing: enough membrane to permeate a dilute counter-flow of
%% order Q_b into the draw (A_PRO = Q_dilute / J_w). The fixed-area placeholder
%% (par.pro.A_PRO_per_Qb) grossly under-sized the stack; physical sizing scales
%% area with the available mixing flow, giving MW-class power as in PRO pilots.
Q_dil = par.sge.avail_dilute .* Q_b;
A_PRO = Q_dil ./ max(J_w, 1e-9);
P_gross = W_d .* A_PRO .* par.pro.eta_turb;
P_pump  = (Q_b .* dP) ./ par.red.eta_pump .* (1 - par.pro.eta_ERD);

%% Thermodynamic admissibility cap: PRO draws on the brine/effluent mixing
%% free energy; its net output cannot exceed eta_II_max*dG_mix. The dilute
%% stream is the low-salinity effluent (par.sge.S_dilute), consistent with the
%% osmotic driving force used above.
R = par.const.R; T = stream.T_f;
Vw = par.const.M_w/par.const.rho_w;
m_d = (par.sge.S_dilute./(1000-par.sge.S_dilute))./par.const.M_NaCl;        % mol/kg
cF  = m_d .* par.const.rho_w ./ (1 + par.sge.S_dilute/1000);                % mol/m^3
cD = stream.c_NaCl_b;
xD = min(max(cD.*Vw,eps),0.5); xF = min(max(cF.*Vw,eps),0.5);
QF = par.sge.avail_dilute .* Q_b; cM = (Q_b.*cD + QF.*cF)./(Q_b+QF); xM = min(max(cM.*Vw,eps),0.5);
dG_mix_rate = max(R.*T.*( Q_b.*cD.*log(xD./xM) + QF.*cF.*log(xF./xM) ), eps);
P_gross_capped = min(P_gross, par.pro.eta_II_max .* dG_mix_rate);
P_net   = par.red.eta_DC .* P_gross_capped - P_pump;
P_net   = max(P_net, 0);
eta_II  = min(P_net ./ dG_mix_rate, par.pro.eta_II_max);

% Economically-relevant membrane area: only enough to extract the capped
% power at the operating areal density W_d (the gradient, not area, limits).
A_econ = min(A_PRO, P_gross_capped ./ max(W_d, eps));

out.J_w = J_w; out.J_s = J_s; out.W_d = W_d; out.W_d_max = W_d_max;
out.A_PRO = A_PRO; out.A_econ = A_econ; out.P_gross = P_gross; out.P_gross_capped = P_gross_capped;
out.P_pump = P_pump; out.P_net = P_net; out.dpi_eff = dpi_eff;
out.dG_mix_rate = dG_mix_rate; out.eta_II = eta_II;
end
