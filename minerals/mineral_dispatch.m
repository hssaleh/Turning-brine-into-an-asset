function mineral = mineral_dispatch(stream, par)
%MINERAL_DISPATCH Dispatcher aggregating all M1-M5 resource-recovery
%pathways into product-mass vectors, energy demand and revenue rate.
%
% Purpose: Runs precip_sequence (M1), crystallizer (M2), lithium_recovery
%          (M3), trace_recovery (M4) and zld_water (M5), then assembles
%          the annual product-mass vector mprod_yr, mineral-train energy
%          demand E_min_yr, reagent cost rate and gross (market-
%          uncapped) revenue rate, per the M20_minerals mapping in
%          Sec. 4 and the energy roll-up of Eq. (6.3).
% Inputs : stream - struct from M00_streams.m
%          par    - parameter struct
% Outputs: mineral (struct):
%             mprod_yr [N x 7 t/yr] columns [CaSO4,MgOH2,NaCl,KCl,
%               Li2CO3,Br,Rb]
%             E_min_yr [N x 1 kWh/yr], E_ZLD_yr [N x 1 kWh/yr]
%             reagent_cost_yr [N x 1 USD/yr]
%             revenue_yr [N x 1 USD/yr] (gross, pre market-absorption cap)
%             plus pass-through sub-structs m1,m2,m3,m4,m5
% Reference: model Sec. 4, Eq. (6.3).

f_av = par.stream.f_av;
ann = 8760*3600*f_av/1000;     % kg/s -> t/yr factor: *3600*8760*fav /1000 (kg->t)

m1 = precip_sequence(stream, par);
m2 = crystallizer(stream, par);
m3 = lithium_recovery(stream, par);
m4 = trace_recovery(stream, par);
m5 = zld_water(stream, par);

% Columns: CaSO4, MgOH2 (M1 gypsum/brucite) | NaCl, KCl (M2 crystalliser) |
% Li2CO3 (M3) | Br, Rb (M4). Halite/sylvite come from the evaporative
% crystalliser, not from raw-brine precipitation (which is far below halite
% saturation), so M1 contributes only the gypsum and brucite columns.
mprod_yr = [m1.mdot_solid(:,1:2), m2.mdot_NaCl, m2.mdot_KCl, ...
            m3.mdot_Li2CO3, m4.mdot_Br, m4.mdot_Rb] .* ann;  % t/yr

E_min_yr = (m1.P_M1 + m2.E_elec + m3.P_M3) .* 8760 .* f_av;  % kWh/yr (M2 electric share only)
E_ZLD_yr = m5.P_M5 .* 8760 .* f_av;                          % kWh/yr

reagent_cost_yr = (m1.mdot_reagent(:,1).*par.m1.p_lime) .* 8760.*3600.*f_av;   % lime for Mg(OH)2 stage

prices = [par.m1.p_CaSO4, par.m1.p_MgOH2, par.m1.p_NaCl, par.m1.p_KCl, par.m3.p_Li2CO3, par.m4.p_Br, par.m4.p_Rb]; % USD/kg
revenue_yr = (mprod_yr*1000) * prices(:);   % t/yr*1000 -> kg/yr, USD/yr

mineral.mprod_yr = mprod_yr;
mineral.product_names = {'CaSO4_2H2O','Mg_OH_2','NaCl','KCl','Li2CO3','Br','Rb'};
mineral.E_min_yr = E_min_yr; mineral.E_ZLD_yr = E_ZLD_yr;
mineral.reagent_cost_yr = reagent_cost_yr; mineral.revenue_yr = revenue_yr;
mineral.Q_water_add_yr = m5.Q_water_add .* 8760.*3600.*f_av;  % m3/yr
mineral.S_residual = m5.S_residual; mineral.Q_res = m5.Q_res;
mineral.m1 = m1; mineral.m2 = m2; mineral.m3 = m3; mineral.m4 = m4; mineral.m5 = m5;
end
