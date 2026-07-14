function par = config_parameters()
%CONFIG_PARAMETERS Master parameter struct for the brine-valorization model.
%
% Project   : Turning Brine into an Asset - Global Techno-Economic and
%             Environmental Assessment of Energy and Mineral Recovery
%             from Desalination Brine
% Purpose   : Single source of truth for every constant, default value and
%             uncertainty range used by the brine-valorization model
%             (Sections 1-10 of Brine_Valorization_Mathematical_Model.md).
% Author    : <AUTHOR>
% Date      : <DATE>
% Version   : 1.0
%
% Description:
%   Returns a nested struct PAR with one field group per model section.
%   No other file in the project hard-codes a numeric constant; every
%   downstream module receives PAR (or a Monte-Carlo-perturbed copy of
%   it) and reads its inputs from it. This keeps the model auditable and
%   reproducible: changing a single default here propagates everywhere.
%
% Inputs  : none.
% Outputs : par (struct) - see field groups below.
% Assumptions:
%   - Default values are literature-centred baselines per the model
%     specification; uncertainty ranges (par.unc) are P5-P95 envelopes
%     used by the Monte-Carlo / Sobol modules, not point predictions.
% References:
%   Sharqawy, Lienhard & Zubair (2010) Desalin. Water Treat. 16, 354-380.
%   Pitzer (1973) J. Phys. Chem. 77, 268.
%   Brine_Valorization_Mathematical_Model.md, Sections 1-10.

par = struct();

%% ---------------------------------------------------------- 1.1 Universal constants
par.const.R        = 8.314462618;      % J mol^-1 K^-1, universal gas constant
par.const.F        = 96485.332;        % C mol^-1, Faraday constant
par.const.N_A      = 6.02214076e23;    % mol^-1, Avogadro constant
par.const.g        = 9.80665;          % m s^-2
par.const.T_ref    = 298.15;           % K
par.const.P_ref    = 101325;           % Pa
par.const.rho_w    = 997.0;            % kg m^-3
par.const.M_w      = 0.0180153;        % kg mol^-1
par.const.M_NaCl   = 0.0584428;        % kg mol^-1
par.const.kWh_to_J = 3.6e6;            % J per kWh

%% ---------------------------------------------------------- 1.3 Seawater ionic composition (35 g/kg reference)
% Order: Cl, Na, SO4, Mg, Ca, K, HCO3, Br, Sr, B, Li, Rb
par.ions.name = {'Cl','Na','SO4','Mg','Ca','K','HCO3','Br','Sr','B','Li','Rb'};
par.ions.w35  = [19.345, 10.784, 2.712, 1.284, 0.4119, 0.3991, 0.1424, 0.0672, 0.00794, 0.00450, 1.8e-4, 1.2e-4]; % g/kg seawater @ S=35
par.ions.M    = [0.035453, 0.022990, 0.096063, 0.024305, 0.040078, 0.039098, 0.061017, 0.079904, 0.087620, 0.010811, 0.006941, 0.085468]; % kg/mol
par.ions.z    = [-1, 1, -2, 2, 2, 1, -1, -1, 2, 0, 1, 1];

%% ---------------------------------------------------------- 2.x Stream-state defaults
par.stream.r       = 0.45;     % water recovery ratio (-)
par.stream.S_f      = 35;       % feed salinity g/kg
par.stream.T_f      = 298.15;   % feed/brine temperature K
par.stream.T0       = 298.15;   % dead-state temperature K
par.stream.P_b      = 5.5e6;    % brine gauge pressure, SWRO, Pa
par.stream.f_av     = 0.90;     % plant availability factor (-)
par.stream.S_p      = 0.30;     % permeate salinity g/kg
par.stream.R_salt   = 0.995;    % salt rejection (-)
par.stream.A_DH     = 0.509;    % Davies/Debye-Hueckel coefficient, kg^0.5 mol^-0.5 @ 25C
par.stream.i_vh     = 1.8;      % van't Hoff factor, NaCl-equivalent
par.stream.nu_NaCl  = 2;        % ions per formula unit

%% ---------------------------------------------------------- 3.1 RED parameters
par.red.alpha        = 0.90;     % membrane permselectivity (-)
par.red.R_CEM        = 2.0e-4;   % Ohm m^2
par.red.R_AEM        = 2.0e-4;   % Ohm m^2
par.red.d_D          = 1.5e-4;   % m
par.red.d_F          = 1.5e-4;   % m
par.red.epsilon_sp   = 0.825;    % spacer porosity (-)
par.red.Lambda_NaCl  = 1.26e-2;  % S m^2 mol^-1
par.red.eta_pump     = 0.75;     % (-)
par.red.eta_DC       = 0.95;     % (-)
par.red.v            = 0.01;     % m s^-1 cross-flow velocity
par.red.L_path       = 0.20;     % m
par.red.d_h          = 3.0e-4;   % m
par.red.vdot_cp       = 1.0e-6;   % m^3 s^-1 per cell pair, channel sizing flow (assumed)
par.red.eta_II_max    = 0.40;     % max realisable 2nd-law efficiency vs Gibbs ceiling (thermodynamic admissibility cap)
% Low-salinity dilute counter-stream for RED/PRO. Coastal salinity-gradient
% power couples the brine (draw) to a low-salinity stream - typically treated
% secondary wastewater effluent - not full seawater; this sets the realisable
% gradient (model Sec. 3 "low-salinity counter-stream"). g/kg.
par.sge.S_dilute      = 8.0;      % brackish/diluted counter-stream salinity, g/kg
                                  % (sets RED P_d in 1-4 W/m2 and PRO W_d,max in 3-10 W/m2 lit. envelopes)
par.sge.avail_dilute  = 1.0;      % dilute-stream availability as a multiple of Q_b (-)

%% ---------------------------------------------------------- 3.2 PRO parameters
par.pro.A_w      = 1.0e-12;   % m s^-1 Pa^-1
par.pro.B_s      = 1.0e-7;    % m s^-1
par.pro.S_struct = 3.0e-4;    % m
par.pro.D_s      = 1.5e-9;    % m^2 s^-1
par.pro.k_d      = 4.0e-5;    % m s^-1
par.pro.eta_turb = 0.85;      % (-)
par.pro.eta_ERD  = 0.96;      % (-)
par.pro.A_PRO_per_Qb = 6.0;   % m^2 membrane per (m^3 s^-1) of brine treated (sizing assumption)
par.pro.eta_II_max   = 0.40;  % max realisable 2nd-law efficiency vs Gibbs ceiling (thermodynamic admissibility cap)

%% ---------------------------------------------------------- 3.3 ERD parameters
par.erd.eta_ERD       = 0.96;  % isobaric exchanger efficiency (-)
par.erd.eta_turb_ERD  = 0.88;  % Pelton turbine efficiency (-)

%% ---------------------------------------------------------- 3.4 Thermal recovery parameters
par.thermal.eta_II_ORC = 0.55;     % second-law ORC efficiency fraction (-)
par.thermal.h_fg       = 2.40e6;   % J/kg at 40 C
par.thermal.GOR_md     = 3.0;      % membrane-distillation gained-output ratio (-)
par.thermal.Tb_by_tech = struct('SWRO', 298.15, 'MED', 313.15, 'MSF', 313.15); % K (25C, 40C, 40C)

%% ---------------------------------------------------------- 3.5 Electrochemical (chlor-alkali) parameters
par.ec.E0       = 2.19;        % V
par.ec.eta_act  = 0.40;        % V
par.ec.R_area   = 1.5e-4;      % Ohm m^2
par.ec.i_cd     = 4000;        % A m^-2
par.ec.eta_F    = 0.97;        % Faradaic efficiency (-)
par.ec.M_H2     = 2.016e-3;    % kg/mol
par.ec.M_Cl2    = 0.070906;    % kg/mol
par.ec.M_NaOH   = 0.039997;    % kg/mol
par.ec.z_H2     = 2;
par.ec.p_H2     = 4.0;         % USD/kg
par.ec.p_Cl2    = 0.20;        % USD/kg
par.ec.p_NaOH   = 0.40;        % USD/kg
par.ec.A_cell_per_Qb = 50;     % m^2 cell area per (m^3 s^-1) of brine routed to EC (sizing assumption)
par.ec.frac_brine_to_EC = 0.05; % fraction of brine flow assigned to the electrochemical pathway (option, not primary)

%% ---------------------------------------------------------- 4.x Mineral-recovery parameters
% Solubility products (mol/kg)^n
par.m1.Ksp_CaSO4   = 3.14e-5;
par.m1.Ksp_CaCO3   = 3.36e-9;
par.m1.Ksp_MgOH2   = 5.61e-12;
par.m1.Ksp_NaCl    = 37.3;
par.m1.Ksp_KCl     = 8.7;
par.m1.eta_rec     = 0.85;    % stage recovery efficiency (-)
par.m1.f_ex        = 1.10;    % reagent excess factor (-)
par.m1.e_mix       = 0.05;    % kWh/kg
par.m1.e_sep       = 0.10;    % kWh/kg
par.m1.e_dry       = 0.20;    % kWh/kg
par.m1.p_CaSO4     = 0.03;    % USD/kg
par.m1.p_MgOH2     = 0.50;    % USD/kg
par.m1.p_NaCl      = 0.06;    % USD/kg
par.m1.p_KCl       = 0.35;    % USD/kg
par.m1.p_Na2CO3    = 0.25;    % USD/kg reagent (carbonate, retained for reference)
par.m1.p_NaOH_r    = 0.40;    % USD/kg reagent (caustic, retained for reference)
par.m1.p_lime      = 0.12;    % USD/kg Ca(OH)2 reagent (industrial Mg(OH)2 precipitant)
par.m1.M_CaSO4     = 0.172172; % kg/mol (CaSO4.2H2O)
par.m1.M_MgOH2     = 0.058320; % kg/mol
par.m1.M_CaOH2     = 0.074093; % kg/mol (lime)
par.m1.M_NaCl      = 0.0584428;% kg/mol
par.m1.M_KCl       = 0.074551; % kg/mol

par.m2.GOR_cryst   = 8;       % (-)
par.m2.SEC_MVC     = 11;      % kWh/m^3 electrical (MVC alternative)
par.m2.m_sat_NaCl  = 6.15;    % mol/kg
par.m2.m_sat_KCl   = 4.8;     % mol/kg
par.m2.eta_NaCl    = 0.80;    % halite crystallisation recovery (Na-limited) (-)
par.m2.eta_KCl     = 0.55;    % sylvite fractional-crystallisation recovery (K-limited) (-)
par.m2.f_heat      = 0.70;    % fraction of evaporation duty met by low-grade/solar heat
par.m2.p_heat      = 0.02;    % USD/kWh thermal (low-grade heat tariff, << electricity)

par.m3.q_max       = 35;      % mg/g sorbent
par.m3.K_L         = 1.5;     % L/mg
par.m3.eta_Li      = 0.70;    % overall recovery efficiency (-)
par.m3.M_Li        = 0.006941; % kg/mol
par.m3.M_Li2CO3    = 0.073891; % kg/mol
par.m3.SEC_Li      = 40;      % kWh per kg Li2CO3
par.m3.p_Li2CO3    = 15;      % USD/kg
par.m3.c_proc_Li   = 4.0;     % USD/kg Li2CO3, DLE sorbent/reagent/processing OPEX
par.m3.c_Li35      = 0.18;    % mg/L Li in seawater at S=35 (reference)

par.m4.eta_Br      = 0.60;    % (-)
par.m4.eta_Rb      = 0.15;    % (-) conservative (dilute, speculative)
par.m4.p_Br        = 3.5;     % USD/kg
par.m4.p_Rb        = 25;      % USD/kg (bulk-equivalent, not specialty-grade)
par.m4.c_proc_Br   = 0.6;     % USD/kg Br, stripping/oxidation processing OPEX

par.m5.r_ZLD       = 0.50;    % incremental ZLD recovery (-)
par.m5.SEC_ZLD     = 12;      % kWh/m^3
par.m5.p_water     = 0.8;     % USD/m^3

%% ---------------------------------------------------------- 4.1 Pitzer-lite binary parameters (Na-Cl dominant brine)
par.pitzer.A_phi = 0.391;    % kg^0.5 mol^-0.5 @ 25C
par.pitzer.b     = 1.2;      % kg^0.5 mol^-0.5
par.pitzer.NaCl.beta0 = 0.0765;  par.pitzer.NaCl.beta1 = 0.2664;  par.pitzer.NaCl.Cphi = 0.00127;
par.pitzer.KCl.beta0  = 0.04835; par.pitzer.KCl.beta1  = 0.2122;  par.pitzer.KCl.Cphi  = -0.00084;
par.pitzer.MgCl2.beta0= 0.35235; par.pitzer.MgCl2.beta1= 1.6815;  par.pitzer.MgCl2.Cphi= 0.00519;
par.pitzer.CaCl2.beta0= 0.31590; par.pitzer.CaCl2.beta1= 1.6140;  par.pitzer.CaCl2.Cphi= -0.00034;
par.pitzer.Na2SO4.beta0=0.01958; par.pitzer.Na2SO4.beta1=1.1130;  par.pitzer.Na2SO4.Cphi=0.00497;

%% ---------------------------------------------------------- 6.x Technology specific energy of desalination
par.tech.SEC_des = struct('SWRO', 3.0, 'MED', 10.0, 'MSF', 15.0); % kWh/m^3 (thermal-equivalent for MED/MSF)

%% ---------------------------------------------------------- 7.x Techno-economic parameters
par.econ.i              = 0.08;   % discount rate (-)
par.econ.N              = 25;     % project lifetime, yr
par.econ.CEPCI_ratio    = 1.15;   % CEPCI_now/CEPCI_ref (-)
par.econ.f_install      = 1.45;   % installation/indirect factor (-)
par.econ.n_u            = 0.60;   % six-tenths scaling exponent (-)
par.econ.p_elec         = 0.08;   % USD/kWh operating electricity price
par.econ.p_elec_credit  = 0.08;   % USD/kWh credit price
par.econ.c_mem          = 15;     % USD/m^2, RED/PRO membrane cost
par.econ.t_mem_life     = 7;      % yr
par.econ.C_maint_frac   = 0.02;   % fraction of CAPEX per year
par.econ.LCOW_base      = 0.80;   % USD/m^3, baseline SWRO LCOW
par.econ.C_labor_frac_opex = 0.05; % fraction of OPEX attributed to labor (sizing assumption)
% Reference unit costs for the six-tenths capital scaling, Eq. 7.1 (USD, reference size)
% RED/PRO stack reference cost recalibrated to ~250/300 USD/m2 installed
% (literature stack+membrane+BoP), replacing earlier ~830-1000 USD/m2 values.
par.econ.cap.RED.C_ref  = 7.5e5;  par.econ.cap.RED.Size_ref  = 5000;  % ~250 USD/m2 installed
par.econ.cap.PRO.C_ref  = 9.0e5;  par.econ.cap.PRO.Size_ref  = 5000;  % ~300 USD/m2 installed
par.econ.cap.ERD.C_ref  = 0.8e6;  par.econ.cap.ERD.Size_ref  = 1.0;   % USD per (m^3/s) brine
par.econ.cap.ORC.C_ref  = 1.2e6;  par.econ.cap.ORC.Size_ref  = 1.0e6; % USD per W net
par.econ.cap.EC.C_ref   = 2.0e6;  par.econ.cap.EC.Size_ref   = 1000;  % USD per m^2 cell area
par.econ.cap.M1.C_ref   = 3.0e7;  par.econ.cap.M1.Size_ref   = 1.0;   % USD per (kg/s) solids (precip+dewater+dry train)
par.econ.cap.M2.C_ref   = 4.0e7;  par.econ.cap.M2.Size_ref   = 1.0;   % USD per (m^3/s) water removed (evap crystalliser)
par.econ.cap.M3.C_ref   = 3.0e7;  par.econ.cap.M3.Size_ref   = 1.0;   % USD per (kg/s) Li2CO3 (capital-intensive DLE)
par.econ.cap.M4.C_ref   = 1.5e7;  par.econ.cap.M4.Size_ref   = 1.0;   % USD per (kg/s) trace product
par.econ.cap.M5.C_ref   = 1.0e7;  par.econ.cap.M5.Size_ref   = 1.0;   % USD per (m^3/s) added water

%% ---------------------------------------------------------- 8.x LCA / environmental parameters
par.lca.EF_grid       = 0.45;  % kg CO2e/kWh
par.lca.EF_prim_NaCl  = 0.12;  % kg CO2e/kg
par.lca.EF_prim_MgOH2 = 1.20;  % kg CO2e/kg
par.lca.EF_prim_Li2CO3= 15;    % kg CO2e/kg
par.lca.EF_water      = 0.4;   % kg CO2e/m^3
par.lca.EF_membrane   = 8;     % kg CO2e/m^2
par.lca.EF_steel      = 1.9;   % kg CO2e/kg
par.lca.E_emb_frac    = 0.05;  % embodied energy as fraction of E_dem_yr (EROI bookkeeping, assumption)

%% ---------------------------------------------------------- 9.x Market-absorption constraint (design-document Moat 3)
% Caps each product's global recovery at a share of its addressable market
% so that bulk low-value salts cannot be reported at theoretical (market-
% crashing) scale. This is the paper's "market-absorption-constrained
% valorization" concept (Brine_Valorization_Paper_Design.docx, Sec. 3.4).
par.market.global_demand_t_yr = struct('NaCl', 3.0e8, 'MgOH2', 1.0e7, 'KCl', 6.5e7, ...
                                        'Li2CO3', 1.2e6, 'Br', 4.0e5);   % t/yr, indicative global demand
par.market.addressable_share  = struct('NaCl', 0.02, 'MgOH2', 0.05, 'KCl', 0.03, ...
                                        'Li2CO3', 0.15, 'Br', 0.10);     % fraction of demand desal brine may serve
par.market.global_C_des_total_m3day = 95e6;   % m3/day, indicative installed global desalination capacity
                                               % (IDA/GWI desalination inventory order-of-magnitude, ~2023;
                                               % used solely to scale the N_plant sample to a global total)

%% ---------------------------------------------------------- 10.x Uncertainty / sensitivity parameters
par.unc.N_MC  = 600;     % Monte-Carlo draws (reduced from literature default 10,000 for tractable runtime; documented deviation)
par.unc.N_sob = 256;     % Sobol base samples (reduced from 4,096 for tractable runtime; documented deviation)
par.unc.seed  = 20260627;

% Uncertain-input table: name, distribution, P5, P95 (Section 10.1)
par.unc.vars = struct( ...
    'p_elec',    struct('dist','lognormal','P5',0.04,'P95',0.18), ...
    'p_Li2CO3',  struct('dist','lognormal','P5',8,   'P95',35), ...
    'p_NaCl',    struct('dist','triangular','P5',0.03,'P95',0.10), ...
    'c_mem',     struct('dist','triangular','P5',8,  'P95',30), ...
    'eta_Li',    struct('dist','uniform','P5',0.5,  'P95',0.85), ...
    'EF_grid',   struct('dist','uniform','P5',0.05,  'P95',0.85), ...
    'i',         struct('dist','triangular','P5',0.05,'P95',0.12) );

%% ---------------------------------------------------------- 9.x Composite viability-index weights
par.viability.w_fin    = 0.25;   % financial (NPV/CAPEX) component weight
par.viability.w_thermo = 0.25;   % thermodynamic (EROI) component weight
par.viability.w_cost   = 0.25;   % net-LCOW-vs-baseline component weight
par.viability.w_env    = 0.25;   % environmental (net GWP) component weight

%% ---------------------------------------------------------- Plant inventory generation controls
par.inv.N_plant = 300;   % synthetic representative global inventory size (documented substitute for DesalData)
par.inv.seed    = 7;

end
