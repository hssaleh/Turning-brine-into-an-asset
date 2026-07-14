function econ = techno_economics(plant, stream, energy, mineral, par)
%TECHNO_ECONOMICS Marginal, portfolio-optimised net-LCOW techno-economic
%assessment of the brine-valorization retrofit.
%
% Purpose: Implements Eqs. (7.1)-(7.18) with the paper's central
%          economic innovation: instead of capitalising every pathway on
%          every plant (which buries viable high-value recovery under
%          uneconomic bulk-evaporation cost), each of the ten pathways
%          (RED, PRO, ERD, ORC, EC, M1-M5) is evaluated on its own
%          marginal cash flow and only included in a plant's retrofit if
%          its standalone marginal NPV is positive. The reported plant
%          net-LCOW, NPV, IRR, payback and BCR therefore describe the
%          OPTIMAL valorization portfolio, and the difference between the
%          portfolio NPV and the best single-pathway NPV defines the
%          "integration premium" (Brine_Valorization_Paper_Design.docx).
% Inputs : plant   - inventory table (LCOW_base, p_elec per region)
%          stream  - struct from M00_streams.m
%          energy  - struct from energy_dispatch.m
%          mineral - struct from mineral_dispatch.m
%          par     - parameter struct
% Outputs: econ (struct), N x 1 unless noted:
%             Per-pathway [N x 10]: CAPEX_k, Rev_k, OPEX_k, ACC_k, NPV_k,
%               sel (logical selection mask). pathway_names {1x10}.
%             Selected-portfolio aggregates: CAPEX_total, ACC_total,
%               OPEX_total, Rev_total, CF_net, dLCOW, LCOW_net, NPV, IRR,
%               payback, payback_disc, BCR, integration_premium.
%             Component revenue/opex used downstream: Rev_energy,
%               Rev_minerals, Rev_water, Rev_EC, OPEX_maint, OPEX_mem,
%               OPEX_reagent.
%             Selected energy balance: E_rec_sel_yr, E_proc_sel_yr,
%               E_net_sel_yr [kWh/yr]; mprod_selected [N x 7 t/yr].
%             LCOE_energy [USD/kWh]; Q_p_m3yr [m3/yr].
% Assumptions: pathway selection uses uncapped (pre-market-absorption)
%          product value - the investment decision is taken on expected
%          merchant prices; the market-absorption haircut is applied to
%          realised global revenue downstream in global_aggregate.m.
%          Maintenance is charged per selected subsystem CAPEX; thermal
%          crystalliser (M2) electricity is charged at the grid price as
%          a conservative proxy (no separate low-grade-heat tariff).
% Reference: model Sec. 7, Eqs. (7.1)-(7.18).

i = par.econ.i; Nlife = par.econ.N; f_av = par.stream.f_av;
n_u = par.econ.n_u; CEPCI = par.econ.CEPCI_ratio; f_inst = par.econ.f_install;
pe = plant.p_elec(:); N = height(plant);
yr = 8760.*f_av;                          % operating hours per year
cap6 = @(Size, C_ref, Size_ref) C_ref .* (max(Size,0)./Size_ref).^n_u .* CEPCI .* f_inst;

%% ---- Eq. 7.1 subsystem six-tenths CAPEX --------------------------------
A_RED  = energy.red.A_econ;               % m^2 membrane (right-sized to capped power)
A_PRO  = energy.pro.A_econ;               % m^2 membrane (right-sized to capped power)
A_EC   = par.ec.A_cell_per_Qb .* par.ec.frac_brine_to_EC .* stream.Q_b;

CAPEX_RED = cap6(A_RED,                          par.econ.cap.RED.C_ref, par.econ.cap.RED.Size_ref);
CAPEX_PRO = cap6(A_PRO,                          par.econ.cap.PRO.C_ref, par.econ.cap.PRO.Size_ref);
CAPEX_ERD = cap6(stream.Q_b,                     par.econ.cap.ERD.C_ref, par.econ.cap.ERD.Size_ref);
CAPEX_ORC = cap6(max(energy.thermal.P_ORC_net,0),par.econ.cap.ORC.C_ref, par.econ.cap.ORC.Size_ref);
CAPEX_EC  = cap6(A_EC,                           par.econ.cap.EC.C_ref,  par.econ.cap.EC.Size_ref);
CAPEX_M1  = cap6(sum(mineral.m1.mdot_solid(:,1:2),2), par.econ.cap.M1.C_ref, par.econ.cap.M1.Size_ref);
CAPEX_M2  = cap6(mineral.m2.V_evap,              par.econ.cap.M2.C_ref,  par.econ.cap.M2.Size_ref);
CAPEX_M3  = cap6(mineral.m3.mdot_Li2CO3,         par.econ.cap.M3.C_ref,  par.econ.cap.M3.Size_ref);
CAPEX_M4  = cap6(mineral.m4.mdot_Br+mineral.m4.mdot_Rb, par.econ.cap.M4.C_ref, par.econ.cap.M4.Size_ref);
CAPEX_M5  = cap6(mineral.m5.Q_water_add,         par.econ.cap.M5.C_ref,  par.econ.cap.M5.Size_ref);

%% ---- per-pathway annual electricity (kWh/yr) --------------------------
E_RED = max(energy.red.P_net,0).*yr/1000;
E_PRO = max(energy.pro.P_net,0).*yr/1000;
E_ERD = max(energy.erd.P_net,0).*yr/1000;
E_ORC = max(energy.thermal.P_ORC_net,0).*yr/1000;
E_M1  = mineral.m1.P_M1 .* yr;            % P in kW -> kWh/yr
E_M2  = mineral.m2.E_elec .* yr;          % crystalliser electric share (kWh/yr)
E_M2heat = mineral.m2.E_heat .* yr;       % crystalliser low-grade-heat share (kWh/yr)
E_M3  = mineral.m3.P_M3 .* yr;
E_M5  = mineral.m5.P_M5 .* yr;
E_EC  = energy.ec.P_demand .* yr/1000;

%% ---- per-pathway attributable revenue (USD/yr) ------------------------
pc = par.econ.p_elec_credit;
mp = mineral.mprod_yr;                    % t/yr [CaSO4,MgOH2,NaCl,KCl,Li2CO3,Br,Rb]
Rev_RED = E_RED.*pc; Rev_PRO = E_PRO.*pc; Rev_ERD = E_ERD.*pc; Rev_ORC = E_ORC.*pc;
Rev_M1  = (mp(:,1).*par.m1.p_CaSO4 + mp(:,2).*par.m1.p_MgOH2).*1000;          % gypsum + brucite
Rev_M2  = (mp(:,3).*par.m1.p_NaCl  + mp(:,4).*par.m1.p_KCl).*1000;            % halite + sylvite
Rev_M3  = mp(:,5).*par.m3.p_Li2CO3.*1000;
Rev_M4  = (mp(:,6).*par.m4.p_Br + mp(:,7).*par.m4.p_Rb).*1000;
Rev_M5  = mineral.Q_water_add_yr .* par.m5.p_water;
Rev_EC  = max( (par.ec.p_H2*energy.ec.mdot_H2 + par.ec.p_Cl2*energy.ec.mdot_Cl2 + ...
                par.ec.p_NaOH*energy.ec.mdot_NaOH).*3600.*yr , 0);   % gross product value (USD/yr)

%% ---- per-pathway attributable OPEX (USD/yr) ---------------------------
mf = par.econ.C_maint_frac; cm = par.econ.c_mem; tl = par.econ.t_mem_life;
OPEX_RED = A_RED.*cm./tl + mf.*CAPEX_RED;
OPEX_PRO = A_PRO.*cm./tl + mf.*CAPEX_PRO;
OPEX_ERD = mf.*CAPEX_ERD;
OPEX_ORC = mf.*CAPEX_ORC;
OPEX_EC  = E_EC.*pe + mf.*CAPEX_EC;
OPEX_M1  = mineral.reagent_cost_yr + E_M1.*pe + mf.*CAPEX_M1;
OPEX_M2  = E_M2.*pe + E_M2heat.*par.m2.p_heat + mf.*CAPEX_M2;   % electric + low-grade-heat duty
OPEX_M3  = E_M3.*pe + par.m3.c_proc_Li.*mp(:,5).*1000 + mf.*CAPEX_M3;   % +DLE processing
OPEX_M4  = par.m4.c_proc_Br.*mp(:,6).*1000 + mf.*CAPEX_M4;              % +Br stripping processing
OPEX_M5  = E_M5.*pe + mf.*CAPEX_M5;

%% ---- assemble [N x 10] pathway matrices -------------------------------
pathway_names = {'RED','PRO','ERD','ORC','EC','M1','M2','M3','M4','M5'};
CAPEX_k = [CAPEX_RED CAPEX_PRO CAPEX_ERD CAPEX_ORC CAPEX_EC CAPEX_M1 CAPEX_M2 CAPEX_M3 CAPEX_M4 CAPEX_M5];
Rev_k   = [Rev_RED Rev_PRO Rev_ERD Rev_ORC Rev_EC Rev_M1 Rev_M2 Rev_M3 Rev_M4 Rev_M5];
OPEX_k  = [OPEX_RED OPEX_PRO OPEX_ERD OPEX_ORC OPEX_EC OPEX_M1 OPEX_M2 OPEX_M3 OPEX_M4 OPEX_M5];

CRF = i.*(1+i).^Nlife ./ ((1+i).^Nlife - 1);
annuity = (1-(1+i).^-Nlife)./i;
ACC_k = CAPEX_k .* CRF;
NPV_k = -CAPEX_k + (Rev_k - OPEX_k).*annuity;

%% ---- Eq. 7.11-7.13 portfolio selection (marginal NPV > 0) -------------
sel = NPV_k > 0;

CAPEX_total = sum(sel.*CAPEX_k, 2);
ACC_total   = CAPEX_total .* CRF;
Rev_total   = sum(sel.*Rev_k, 2);
OPEX_base   = sum(sel.*OPEX_k, 2);
OPEX_labor  = par.econ.C_labor_frac_opex .* OPEX_base ./ max(1-par.econ.C_labor_frac_opex, eps);
OPEX_total  = OPEX_base + OPEX_labor;
CF_net      = Rev_total - OPEX_total;

%% ---- Eq. 7.15-7.16 net-LCOW ------------------------------------------
Q_p_m3yr = stream.Q_p .* 86400 .* 365 .* f_av;
dLCOW = (ACC_total - CF_net) ./ max(Q_p_m3yr, eps);
LCOW_net = plant.LCOW_base(:) + dLCOW;

%% ---- Eq. 7.9-7.10 NPV, Eq. 7.14 BCR, integration premium -------------
NPV = -CAPEX_total + CF_net.*annuity;
BCR = (Rev_total.*annuity) ./ max(CAPEX_total + OPEX_total.*annuity, eps);
NPV_best_single = max(NPV_k.*sel, [], 2);          % best individually-selected pathway
integration_premium = NPV - NPV_best_single;        % value of co-deploying the portfolio

%% ---- Eq. 7.11 IRR (vectorized Newton-Raphson on the selected cash flow)
IRR = nan(N,1);
valid = CF_net > 0 & CAPEX_total > 0;
x = i*ones(N,1);
for it = 1:60
    ann_x = (1-(1+x).^-Nlife)./x;
    f  = -CAPEX_total + CF_net.*ann_x;
    dann = ( Nlife.*(1+x).^(-Nlife-1).*x - (1-(1+x).^-Nlife) ) ./ x.^2;
    step = f./(CF_net.*dann);
    step(~isfinite(step)) = 0;
    x = max(x - step, -0.99);
end
IRR(valid) = x(valid);

%% ---- Eq. 7.12-7.13 simple and discounted payback ---------------------
payback = CAPEX_total ./ max(CF_net, eps); payback(CF_net<=0) = Inf;
k_disc = CAPEX_total.*i./max(CF_net,eps);
payback_disc = log(1./max(1-k_disc,eps)) ./ log(1+i);
payback_disc(CF_net<=0 | k_disc>=1) = Inf;

%% ---- selected energy balance & selected mineral production ------------
selE_rec  = sel(:,1).*E_RED + sel(:,2).*E_PRO + sel(:,3).*E_ERD + sel(:,4).*E_ORC;
selE_proc = sel(:,5).*E_EC + sel(:,6).*E_M1 + sel(:,7).*E_M2 + sel(:,8).*E_M3 + sel(:,10).*E_M5;
E_rec_sel_yr = selE_rec; E_proc_sel_yr = selE_proc; E_net_sel_yr = selE_rec - selE_proc;

mprod_selected = mp;
mprod_selected(:,1:2) = mp(:,1:2).*sel(:,6);     % M1 products (gypsum, brucite)
mprod_selected(:,3:4) = mp(:,3:4).*sel(:,7);     % M2 products (halite, sylvite)
mprod_selected(:,5)   = mp(:,5).*sel(:,8);        % M3 (Li2CO3)
mprod_selected(:,6:7) = mp(:,6:7).*sel(:,9);     % M4 (Br,Rb)

%% ---- component fields used downstream (selected) ----------------------
Rev_energy   = sum(sel(:,1:4).*[Rev_RED Rev_PRO Rev_ERD Rev_ORC],2);
Rev_minerals = sel(:,6).*Rev_M1 + sel(:,7).*Rev_M2 + sel(:,8).*Rev_M3 + sel(:,9).*Rev_M4;
Rev_water    = sel(:,10).*Rev_M5;
Rev_EC_sel   = sel(:,5).*Rev_EC;
OPEX_maint   = sum(sel.*(mf.*CAPEX_k),2);
OPEX_mem     = sel(:,1).*(A_RED.*cm./tl) + sel(:,2).*(A_PRO.*cm./tl);
OPEX_reagent = sel(:,6).*mineral.reagent_cost_yr;

%% ---- Eq. 7.17-7.18 LCOE of the selected energy-recovery train --------
ACC_energy = sum(sel(:,1:4).*ACC_k(:,1:4),2);
LCOE_energy = (ACC_energy + OPEX_mem) ./ max(E_rec_sel_yr, eps);

%% ---- pack ------------------------------------------------------------
econ.pathway_names = pathway_names;
econ.CAPEX_k = CAPEX_k; econ.Rev_k = Rev_k; econ.OPEX_k = OPEX_k;
econ.ACC_k = ACC_k; econ.NPV_k = NPV_k; econ.sel = sel;
econ.CAPEX_RED=CAPEX_RED; econ.CAPEX_PRO=CAPEX_PRO; econ.CAPEX_ERD=CAPEX_ERD;
econ.CAPEX_ORC=CAPEX_ORC; econ.CAPEX_EC=CAPEX_EC; econ.CAPEX_M1=CAPEX_M1;
econ.CAPEX_M2=CAPEX_M2; econ.CAPEX_M3=CAPEX_M3; econ.CAPEX_M4=CAPEX_M4; econ.CAPEX_M5=CAPEX_M5;
econ.CAPEX_total=CAPEX_total; econ.CRF=CRF; econ.ACC_total=ACC_total;
econ.OPEX_total=OPEX_total; econ.OPEX_labor=OPEX_labor;
econ.OPEX_maint=OPEX_maint; econ.OPEX_mem=OPEX_mem; econ.OPEX_reagent=OPEX_reagent;
econ.Rev_total=Rev_total; econ.Rev_energy=Rev_energy; econ.Rev_minerals=Rev_minerals;
econ.Rev_water=Rev_water; econ.Rev_EC=Rev_EC_sel;
econ.CF_net=CF_net; econ.dLCOW=dLCOW; econ.LCOW_net=LCOW_net;
econ.NPV=NPV; econ.IRR=IRR; econ.payback=payback; econ.payback_disc=payback_disc;
econ.BCR=BCR; econ.integration_premium=integration_premium;
econ.NPV_best_single=NPV_best_single;
econ.E_rec_sel_yr=E_rec_sel_yr; econ.E_proc_sel_yr=E_proc_sel_yr; econ.E_net_sel_yr=E_net_sel_yr;
econ.mprod_selected=mprod_selected;
econ.LCOE_energy=LCOE_energy; econ.Q_p_m3yr=Q_p_m3yr;
econ.n_selected = sum(sel,2);
end
