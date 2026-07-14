function out = crystallizer(stream, par)
%CRYSTALLIZER (M2) Evaporative / eutectic crystallisation of halite (NaCl)
%and sylvite (KCl), with the associated water-removal energy.
%
% Purpose: Implements Eqs. (4.8)-(4.9) extended to the actual salt
%          products. Brine is concentrated to halite saturation and NaCl
%          crystallises; the residual liquor is concentrated further and
%          KCl crystallises. Because desalination brine is ~78% NaCl by
%          dissolved mass, halite is the single largest recoverable solid;
%          recovering it requires evaporative concentration, which is why
%          this stage carries a large thermal/electrical energy demand and
%          is the pathway whose economic selection is most sensitive to the
%          availability of low-grade heat.
% Inputs : stream - struct from M00_streams.m (Q_b, rho_b, S_b, m_s_b,
%             m_i_b, T_f)
%          par    - parameter struct (par.m2.*, par.ions.*, par.thermal.*,
%             par.const.*)
% Outputs: out.mdot_NaCl, out.mdot_KCl [kg/s]; out.V_evap [m3/s];
%          out.P_M2 [kW] (thermal-route process power); out.P_M2_MVC [kW]
%          (electrical MVC alternative); out.E_heat, out.E_elec split [kW]
% Assumptions: NaCl recovery is sodium-limited at par.m2.eta_NaCl of the
%          brine Na inventory delivered as NaCl; KCl recovery is potassium-
%          limited at par.m2.eta_KCl. A fraction par.m2.f_heat of the
%          evaporation duty is met by low-grade/solar heat (priced at
%          par.m2.p_heat), the remainder electrically.
% Reference: model Sec. 4.3, Eqs. (4.8)-(4.9); ionic composition Sec. 1.3.

names = par.ions.name;
idxNa = find(strcmpi(names,'Na')); idxK = find(strcmpi(names,'K'));
m_Na = stream.m_i_b(:,idxNa); m_K = stream.m_i_b(:,idxK);     % mol/kg water
conv = stream.rho_b ./ (1+stream.S_b/1000);                   % mol/kg-water -> mol/m^3

%% Eq. 4.8 water removed to reach halite saturation
V_evap = stream.Q_b .* max(1 - stream.m_s_b./par.m2.m_sat_NaCl, 0);   % m3/s

%% Sodium-limited halite and potassium-limited sylvite mass rates
mdot_NaCl = par.m2.eta_NaCl .* stream.Q_b .* (m_Na.*conv) .* par.m1.M_NaCl;  % kg/s
mdot_KCl  = par.m2.eta_KCl  .* stream.Q_b .* (m_K .*conv) .* par.m1.M_KCl;   % kg/s

%% Eq. 4.9 evaporation specific energy (thermal route) and process power
SEC_evap = (1/par.m2.GOR_cryst) .* par.const.rho_w .* par.thermal.h_fg / 3.6e6;  % kWh/m3
P_M2  = V_evap .* SEC_evap .* 3600;            % kW (thermal-equivalent)
P_M2_MVC = V_evap .* par.m2.SEC_MVC .* 3600;   % kW (electrical MVC alternative)

%% Heat / electricity split of the thermal route
E_heat = par.m2.f_heat .* P_M2;                % kW supplied as low-grade heat
E_elec = (1-par.m2.f_heat) .* P_M2;            % kW supplied electrically

out.mdot_NaCl = mdot_NaCl; out.mdot_KCl = mdot_KCl;
out.V_evap = V_evap; out.P_M2 = P_M2; out.P_M2_MVC = P_M2_MVC;
out.E_heat = E_heat; out.E_elec = E_elec; out.SEC_evap = SEC_evap;
end
