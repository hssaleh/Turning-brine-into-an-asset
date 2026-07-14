function out = lithium_recovery(stream, par)
%LITHIUM_RECOVERY (M3) Selective lithium recovery from end-brine.
%
% Purpose: Implements Eqs. (4.10)-(4.13): Langmuir sorbent uptake
%          (diagnostic), recovered lithium mass rate, battery-grade
%          Li2CO3 product, and stage specific energy.
% Inputs : stream - struct from M00_streams.m (Q_b, S_b)
%          par    - parameter struct (par.m3.*)
% Outputs: out.q_Li [mg/g], out.mdot_Li [kg/s], out.mdot_Li2CO3 [kg/s],
%          out.P_M3 [kW]
% Assumptions: brine Li concentration scales linearly with salinity from
%          the seawater mean (c_Li,in = S_b/35 * 0.18 mg/L); recovered Li
%          mass rate uses the overall efficiency directly on the inlet
%          concentration (Eq. 4.11 reduction, c_Li,out implicit in eta_Li).
% Reference: model Sec. 4.4, Eqs. (4.10)-(4.13).

c_Li_in_mgL = (stream.S_b/35) .* par.m3.c_Li35;     % mg/L
c_Li_in = c_Li_in_mgL * 1e-3;                        % kg/m^3 (mg/L == g/m3)

%% Eq. 4.10 Langmuir uptake (diagnostic sorbent loading at inlet concentration)
q_Li = par.m3.q_max .* par.m3.K_L .* c_Li_in_mgL ./ (1 + par.m3.K_L.*c_Li_in_mgL);

%% Eq. 4.11-4.12 recovered Li and Li2CO3 mass rates
mdot_Li = par.m3.eta_Li .* stream.Q_b .* c_Li_in;
mdot_Li2CO3 = mdot_Li .* (par.m3.M_Li2CO3 ./ (2*par.m3.M_Li));

%% Eq. 4.13 stage specific energy -> process power
P_M3 = mdot_Li2CO3 .* par.m3.SEC_Li .* 3600;   % kW

out.q_Li = q_Li; out.mdot_Li = mdot_Li; out.mdot_Li2CO3 = mdot_Li2CO3; out.P_M3 = P_M3;
end
