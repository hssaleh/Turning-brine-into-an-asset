function k = fleet_kpis(plant, R)
%FLEET_KPIS Collapse a per-plant pipeline result into scalar fleet-level
%key performance indicators.
%
% Purpose: Provides the scalar KPI vector tracked by the Monte-Carlo
%          propagation and reported as the deterministic headline
%          results. Capacity-weighted means are used for intensive
%          quantities (net-LCOW, carbon intensity, efficiencies) so that
%          large plants dominate appropriately; extensive quantities
%          (energy, CO2, NPV) use the global-scaled sums already formed
%          in R.agg.global.
% Inputs : plant - inventory table (C_des used as capacity weight)
%          R     - struct from run_pipeline.m
% Outputs: k (struct) of scalar KPIs (see field names below).
% Reference: model Sec. 6-9 reporting layer.

w = plant.C_des(:); w = w./sum(w);   % capacity weights

k.NPV_global_B      = R.agg.global.NPV_B_USD;                 % billion USD (merchant)
k.E_rec_TWh         = R.agg.global.E_rec_TWh_yr;              % TWh/yr
k.CO2_mit_Mt        = R.agg.global.CO2_mitigation_Mt_yr;      % Mt CO2e/yr
k.water_Mm3         = R.agg.global.water_Mm3_yr;              % Mm3/yr
k.MgOH2_Mt_constr   = R.agg.global.mineral_Mt_yr_constr(2);   % Mt/yr (market-absorbable)
k.Li2CO3_kt_constr  = R.agg.global.mineral_Mt_yr_constr(5)*1e3; % kt/yr (market-absorbable)

k.LCOW_net_mean     = sum(w.*R.econ.LCOW_net);                % USD/m3 (merchant)
k.dLCOW_mean        = sum(w.*R.econ.dLCOW);                   % USD/m3
k.CI_water_mean     = sum(w.*R.lca.CI_water);                 % kg CO2e/m3
k.EROI_mean         = sum(w.*R.lca.EROI);                     % -
k.VI_mean           = sum(w.*R.agg.VI);                       % -
k.eta_II_RED_mean   = sum(w.*R.exergy.eta_II_RED);           % -
k.integ_prem_B      = R.agg.GRF.*sum(R.econ.integration_premium)./1e9; % billion USD

k.frac_NPV_pos      = mean(R.econ.NPV > 0);                  % fraction of plants viable
k.frac_LCOW_down    = mean(R.econ.LCOW_net < plant.LCOW_base(:)); % fraction with net-LCOW below baseline
end
