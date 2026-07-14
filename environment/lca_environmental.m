function lca = lca_environmental(plant, stream, energy, mineral, massE, econ, par)
%LCA_ENVIRONMENTAL Life-cycle and avoided-discharge environmental
%assessment of the brine-valorization retrofit.
%
% Purpose: Implements Eqs. (8.1)-(8.13): operational GWP from the net
%          plant electricity balance, embodied GWP from membranes and
%          process equipment, avoided-virgin-production credits for
%          recovered minerals, the avoided-discharge salinity-load
%          indicator, energy-return-on-investment (EROI) and the net
%          carbon intensity of water produced.
% Inputs : plant   - inventory table (EF_grid per region)
%          stream  - struct from M00_streams.m
%          energy  - struct from energy_dispatch.m
%          mineral - struct from mineral_dispatch.m
%          massE   - struct from mass_energy_integration.m
%          econ    - struct from techno_economics.m (CAPEX for embodied
%                    structural-material proxy)
%          par     - parameter struct
% Outputs: lca (struct):
%             GWP_op_yr [kg CO2e/yr] operational (net electricity balance)
%             GWP_emb_yr [kg CO2e/yr] embodied (membranes + process energy)
%             GWP_avoided_min_yr [kg CO2e/yr] avoided virgin-production
%               credit (negative = credit)
%             GWP_net_yr [kg CO2e/yr] = GWP_op_yr + GWP_emb_yr + GWP_avoided_min_yr
%             CO2_mitigation_t_yr [t/yr] = -GWP_net_yr/1000 (positive = net mitigation)
%             CI_water [kg CO2e/m3] carbon intensity of water produced
%             salt_diverted_t_yr [t/yr], f_salinity_avoided [-]
%             EROI [-] energy return on (process + embodied) energy invested
% Assumptions: embodied GWP of process steelwork/structure is proxied as
%          a fixed emission factor per USD of CAPEX is NOT used (no such
%          factor in the source spec); instead embodied energy is taken
%          as par.lca.E_emb_frac of the primary desalination demand
%          (documented bookkeeping assumption) converted to GWP via the
%          regional grid factor, plus explicit membrane-material GWP
%          (par.lca.EF_membrane) annualized over the membrane lifetime.
%          Avoided-discharge salinity load is reported relative to the
%          total dissolved-solids mass rate in the brine stream.
% Reference: model Sec. 8, Eqs. (8.1)-(8.13).

EF_grid = plant.EF_grid(:);
f_av = par.stream.f_av;

%% Eq. 8.1-8.2 operational GWP (net SELECTED-portfolio electricity balance)
GWP_op_yr = -econ.E_net_sel_yr .* EF_grid;   % net consumer (E_net<0) -> positive emissions

%% Eq. 8.3-8.4 embodied GWP: membrane materials (annualized, selected) + embodied process energy
A_mem = econ.sel(:,1).*energy.red.N_cp + econ.sel(:,2).*energy.pro.A_PRO;  % m^2 (only built membranes)
GWP_membrane_yr = A_mem .* par.lca.EF_membrane ./ par.econ.t_mem_life;
E_emb_yr = par.lca.E_emb_frac .* massE.E_dem_yr;                  % kWh/yr, embodied-energy bookkeeping
GWP_emb_yr = GWP_membrane_yr + E_emb_yr.*EF_grid;

%% Eq. 8.5-8.7 avoided virgin-production credit for recovered minerals (selected)
mprod = econ.mprod_selected;   % t/yr [CaSO4,MgOH2,NaCl,KCl,Li2CO3,Br,Rb]
GWP_avoided_min_yr = -1000.*( mprod(:,3).*par.lca.EF_prim_NaCl + ...
                               mprod(:,2).*par.lca.EF_prim_MgOH2 + ...
                               mprod(:,5).*par.lca.EF_prim_Li2CO3 );  % t->kg, negative = credit
GWP_avoided_water_yr = -mineral.Q_water_add_yr .* par.lca.EF_water;   % avoided marginal water production

%% Eq. 8.8 net GWP and CO2-mitigation potential
GWP_net_yr = GWP_op_yr + GWP_emb_yr + GWP_avoided_min_yr + GWP_avoided_water_yr;
CO2_mitigation_t_yr = -GWP_net_yr ./ 1000;

%% Eq. 8.9 carbon intensity of water produced
CI_water = GWP_net_yr ./ max(econ.Q_p_m3yr, eps);

%% Eq. 8.10-8.11 avoided-discharge salinity-load indicator
TDS_discharge_yr = stream.Q_b .* stream.rho_b .* (stream.S_b/1000) .* 8760.*3600.*f_av ./ 1000;  % t/yr
salt_diverted_t_yr = sum(mprod(:,1:4),2);   % t/yr (CaSO4,MgOH2,NaCl,KCl - the bulk salt-load fraction)
f_salinity_avoided = salt_diverted_t_yr ./ max(TDS_discharge_yr, eps);

%% Eq. 8.12-8.13 energy return on investment (selected portfolio)
E_invested_yr = econ.E_proc_sel_yr + E_emb_yr;
EROI = econ.E_rec_sel_yr ./ max(E_invested_yr, eps);

lca.GWP_op_yr = GWP_op_yr; lca.GWP_emb_yr = GWP_emb_yr;
lca.GWP_avoided_min_yr = GWP_avoided_min_yr; lca.GWP_avoided_water_yr = GWP_avoided_water_yr;
lca.GWP_net_yr = GWP_net_yr; lca.CO2_mitigation_t_yr = CO2_mitigation_t_yr;
lca.CI_water = CI_water; lca.TDS_discharge_yr = TDS_discharge_yr;
lca.salt_diverted_t_yr = salt_diverted_t_yr; lca.f_salinity_avoided = f_salinity_avoided;
lca.EROI = EROI; lca.E_invested_yr = E_invested_yr;
end
