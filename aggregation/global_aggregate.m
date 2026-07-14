function agg = global_aggregate(plant, energy, mineral, econ, lca, par)
%GLOBAL_AGGREGATE Composite viability index, electricity-price cost-
%crossover, market-absorption-constrained global potential and the brine-
%resource-atlas backbone.
%
% Purpose: Implements Eqs. (9.1)-(9.7) using a deliberate two-tier logic
%   that is itself a central result of the paper:
%   (A) PLANT-LEVEL VIABILITY is assessed at merchant prices - any single
%       plant is a price-taker whose output the market absorbs, so its
%       NPV / net-LCOW / viability index come straight from the
%       merchant-price techno-economics (econ.*).
%   (B) The GLOBAL resource potential is reported in two versions: the
%       technical potential (selected pathways, merchant) and the
%       MARKET-ABSORPTION-CONSTRAINED potential, where each product's
%       global recovered mass is capped at its addressable global-market
%       share (design-doc "Moat 3"). The ratio of the two exposes which
%       products are demand-limited (e.g. bulk Mg(OH)2 / NaCl) versus
%       chemistry/energy-limited (e.g. Li, Br) - the headline message
%       that the binding global constraint is market demand, not
%       thermodynamics.
% Inputs : plant, energy, mineral, econ, lca - upstream per-plant
%          structs/tables; par - parameter struct.
% Outputs: agg (struct):
%             GRF [-] sample-to-global capacity scale factor
%             f_absorb [1 x 7 -] absorbable fraction per product
%             VI [N x 1, 0-1] composite viability index (+ sub-scores)
%             p_elec_crossover [USD/kWh] per-plant breakeven elec price
%             global.* totals: technical and market-constrained tonnage,
%               energy, water, CO2 mitigation, NPV.
% Assumptions: per-plant viability assumes merchant-price sales (standard
%          price-taker assumption in plant techno-economics); the market
%          cap binds only on the aggregate, signalling the deployable
%          share rather than re-pricing individual plants.
% Reference: model Sec. 9, Eqs. (9.1)-(9.7).

%% Eq. 9.1 sample-to-global capacity scale factor
GRF = par.market.global_C_des_total_m3day ./ sum(plant.C_des);

%% Eq. 9.2-9.3 market-absorption cap on aggregate selected production
product_names = mineral.product_names;   % {CaSO4,MgOH2,NaCl,KCl,Li2CO3,Br,Rb}
demand_map = struct('Mg_OH_2','MgOH2','NaCl','NaCl','KCl','KCl','Li2CO3','Li2CO3','Br','Br');
nP = numel(product_names);
mprod_global_tech = GRF .* sum(mineral.mprod_yr, 1);   % t/yr per product (technical chemistry, all pathways)
f_absorb = ones(1, nP);
for p = 1:nP
    pname = product_names{p};
    if isfield(demand_map, pname)
        key = demand_map.(pname);
        cap_mass = par.market.global_demand_t_yr.(key) .* par.market.addressable_share.(key);
        if mprod_global_tech(p) > 0
            f_absorb(p) = min(1, cap_mass ./ mprod_global_tech(p));
        end
    end
end
mprod_global_constr = mprod_global_tech .* f_absorb;   % t/yr per product (deployable)

%% Eq. 9.5 per-plant electricity-price cost-crossover (merchant, breakeven vs ACC)
Eproc = econ.E_proc_sel_yr;
Lf = par.econ.C_labor_frac_opex;
OPEX_fixed = econ.OPEX_maint + econ.OPEX_mem + econ.OPEX_reagent;
p_elec_crossover = ( (econ.Rev_total - econ.ACC_total).*(1-Lf) - OPEX_fixed ) ./ max(Eproc, eps);
p_elec_crossover(Eproc <= 0) = NaN;

%% Eq. 9.6 composite viability index (merchant per-plant, bounded sub-scores)
VI_fin    = 1 ./ (1 + exp(-econ.NPV ./ max(econ.CAPEX_total, eps)));   % logistic of NPV/CAPEX
VI_thermo = lca.EROI ./ (1 + lca.EROI);                                % Michaelis-Menten of EROI
VI_cost   = plant.LCOW_base(:) ./ max(econ.LCOW_net + plant.LCOW_base(:), eps);
VI_cost   = min(max(VI_cost, 0), 1);
VI_env    = 1 ./ (1 + exp( lca.GWP_net_yr ./ (0.5.*abs(lca.GWP_op_yr) + 1) ));
VI = par.viability.w_fin.*VI_fin + par.viability.w_thermo.*VI_thermo + ...
     par.viability.w_cost.*VI_cost + par.viability.w_env.*VI_env;

%% Eq. 9.7 global resource-potential aggregation
g.GRF = GRF;
g.E_rec_TWh_yr            = GRF .* sum(econ.E_rec_sel_yr) ./ 1e9;       % TWh/yr
g.mineral_Mt_yr_tech      = mprod_global_tech ./ 1e6;                  % Mt/yr (technical)
g.mineral_Mt_yr_constr    = mprod_global_constr ./ 1e6;               % Mt/yr (market-absorbable)
g.product_names           = product_names;
g.water_Mm3_yr            = GRF .* sum(mineral.Q_water_add_yr.*any(econ.sel(:,10),2)) ./ 1e6;
g.CO2_mitigation_Mt_yr    = GRF .* sum(lca.CO2_mitigation_t_yr) ./ 1e6; % Mt CO2e/yr
g.NPV_B_USD               = GRF .* sum(econ.NPV) ./ 1e9;               % billion USD (merchant)
g.revenue_B_USD_tech      = GRF .* sum(econ.Rev_total) ./ 1e9;
g.revenue_B_USD_constr    = (mprod_global_constr * mineralPriceVec(par).') ./ 1e9 ...
                            + GRF.*sum(econ.Rev_energy + econ.Rev_water + econ.Rev_EC)./1e9;

agg.GRF = GRF; agg.f_absorb = f_absorb; agg.product_names = product_names;
agg.mprod_global_tech = mprod_global_tech; agg.mprod_global_constr = mprod_global_constr;
agg.p_elec_crossover = p_elec_crossover;
agg.VI = VI; agg.VI_fin = VI_fin; agg.VI_thermo = VI_thermo; agg.VI_cost = VI_cost; agg.VI_env = VI_env;
agg.global = g;
end

% ----------------------------------------------------------------------
function v = mineralPriceVec(par)
%MINERALPRICEVEC USD/kg price vector aligned to the 7-product mass vector.
v = [par.m1.p_CaSO4, par.m1.p_MgOH2, par.m1.p_NaCl, par.m1.p_KCl, ...
     par.m3.p_Li2CO3, par.m4.p_Br, par.m4.p_Rb] .* 1000;   % USD/t
end
