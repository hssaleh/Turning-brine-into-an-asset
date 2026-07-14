function exergy = exergy_balance(stream, energy, mineral, par)
%EXERGY_BALANCE Second-law (exergy) accounting for the brine-valorization
%plant: physical and chemical stream exergy, minimum separation work,
%subsystem exergy balances, second-law efficiencies and the exergy-
%destruction map.
%
% Purpose: Implements Eqs. (5.1)-(5.13) of the brine-valorization model.
% Inputs : stream  - struct from M00_streams.m
%          energy  - struct from energy_dispatch.m
%          mineral - struct from mineral_dispatch.m
%          par     - parameter struct
% Outputs: exergy (struct):
%             ex_ph_b, ex_ch_b [J/kg]; Exdot_b [W] (Eq. 5.5)
%             Wdot_min_sep [W] (Eq. 5.4, chemical-exergy-flux proxy)
%             eta_II_RED/PRO/ERD/th/min [-] (Eqs. 5.9-5.13)
%             Exdot_dest [N x 5 W] destruction by subsystem
%               [RED, PRO, ERD, Thermal, Minerals]
% Assumptions: dead state = receiving seawater at (T0, S_f). Physical
%          exergy uses the technology brine temperature T_b from the
%          thermal module; chemical exergy uses a binary water-salt
%          mole-fraction reduction (Eq. 5.3) with NaCl-equivalent molar
%          mass. Mixing exergy available to RED and PRO is treated as
%          the shared Gibbs-mixing rate (Eq. 3.11); reagent exergy
%          (Eq. 5.13 denominator) is neglected (documented omission -
%          reagent chemical exergy requires a full formation-enthalpy
%          database not included in this model). Loss terms in the
%          subsystem exergy balance (Eq. 5.6) are taken as zero so that
%          Exdot_dest = Exdot_in - Exdot_out (non-negative by clipping).
% Reference: model Sec. 5, Eqs. (5.1)-(5.13).

R = par.const.R; T0 = stream.T0; Mw = par.const.M_w; M_NaCl = par.const.M_NaCl;

%% Eq. 5.2 physical exergy (brine at technology temperature T_b vs dead state T0)
T_b = energy.thermal.T_b;
h_b = h_sw(stream.S_b, T_b, T0);
s_b = s_sw(stream.S_b, T_b, T0);
ex_ph_b = h_b - T0.*s_b;     % J/kg, = 0 when T_b = T0 (SWRO)

%% Eq. 5.3-5.4 chemical exergy of brine relative to the seawater dead state.
% The work-producing process for a stream MORE concentrated than the
% environment is its reversible mixing (dilution) with the receiving sea;
% the recoverable work equals the (positive) Gibbs free energy of mixing of
% the brine with seawater. This is the salinity-gradient exergy and is the
% physically correct, sign-consistent chemical-exergy resource (the earlier
% mole-fraction unmixing form returned a small negative value that clipped
% to zero).
mdot_b = stream.Q_b .* stream.rho_b;
Vw_molar = Mw/par.const.rho_w;
c_b = stream.c_NaCl_b; c_0 = stream.c_NaCl_f;          % seawater dead state
x_b = min(max(c_b.*Vw_molar,eps),0.5); x_0 = min(max(c_0.*Vw_molar,eps),0.5);
Q0 = stream.Q_b; cM = (stream.Q_b.*c_b + Q0.*c_0)./(stream.Q_b+Q0);
xM = min(max(cM.*Vw_molar,eps),0.5);
Exdot_chem_b = max(R.*T0.*( stream.Q_b.*c_b.*log(x_b./xM) + Q0.*c_0.*log(x_0./xM) ), 0);  % W
ex_ch_b = Exdot_chem_b ./ max(mdot_b, eps);            % J/kg (positive)

%% Eq. 5.5 brine physical-exergy rate (pressure + thermal components)
Exdot_b = stream.Q_b.*stream.P_b + mdot_b.*stream.cp_b.*( (T_b-T0) - T0.*log(max(T_b,T0+1e-6)./T0) );

%% Eq. 5.4 minimum separation work rate = chemical (salinity-gradient) exergy rate
Wdot_min_sep = Exdot_chem_b;

%% Eq. 5.9-5.13 pathway second-law efficiencies (each vs its own exergy fuel)
dG_RED = energy.red.dG_mix_rate; dG_PRO = energy.pro.dG_mix_rate;
eta_II_RED = energy.red.eta_II;   % already capped at par.red.eta_II_max
eta_II_PRO = energy.pro.eta_II;
eta_II_ERD = min(energy.erd.P_net ./ max(stream.Q_b.*stream.P_b, eps), 1);
eta_II_th  = min(energy.thermal.P_ORC_net ./ max(energy.thermal.W_Carnot, eps), 1);

E_min_elec_W = mineral.E_min_yr .* 1000 ./ (8760.*par.stream.f_av);  % kWh/yr -> W
eta_II_min = min(Wdot_min_sep ./ max(E_min_elec_W, eps), 5);   % can exceed 1 if E_min_elec underestimates true exergy input (documented)

%% Exergy-destruction map (Eq. 5.6, 5.8; zero-loss reduction)
dest_RED = max(dG_RED - energy.red.P_net, 0);
dest_PRO = max(dG_PRO - energy.pro.P_net, 0);
dest_ERD = max(stream.Q_b.*stream.P_b - energy.erd.P_net, 0);
dest_TH  = max(energy.thermal.W_Carnot - energy.thermal.P_ORC_net, 0);
dest_MIN = max(E_min_elec_W - Wdot_min_sep, 0);

exergy.ex_ph_b = ex_ph_b; exergy.ex_ch_b = ex_ch_b; exergy.Exdot_b = Exdot_b;
exergy.Wdot_min_sep = Wdot_min_sep;
exergy.eta_II_RED = eta_II_RED; exergy.eta_II_PRO = eta_II_PRO;
exergy.eta_II_ERD = eta_II_ERD; exergy.eta_II_th = eta_II_th; exergy.eta_II_min = eta_II_min;
exergy.Exdot_dest = [dest_RED, dest_PRO, dest_ERD, dest_TH, dest_MIN];
exergy.dest_labels = {'RED','PRO','ERD','Thermal','Minerals'};
end
