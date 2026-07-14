function out = precip_sequence(stream, par)
%PRECIP_SEQUENCE (M1) Selective sequential precipitation of Ca, Mg, Na, K.
%
% Purpose: Implements Eqs. (4.4)-(4.7): saturation-index-driven recovery
%          of gypsum/calcite, brucite, halite and sylvite from brine, in
%          ascending-solubility order, with reagent demand and stage
%          specific energy.
% Inputs : stream - struct from M00_streams.m (m_i_b, Q_b, rho_b, S_b,
%             gamma_b)
%          par    - parameter struct (par.m1.*, par.ions.*)
% Outputs: out (struct):
%             mdot_solid [N x 4 kg/s] for [CaSO4.2H2O, Mg(OH)2, NaCl, KCl]
%             mdot_reagent [N x 2 kg/s] for [Na2CO3, NaOH]
%             P_M1 [N x 1 kW] stage process power (mixing+separation+drying)
%             SI [N x 4 -] saturation indices
% Assumptions: each stage's driving ion content is evaluated on the
%          as-received brine composition (not sequentially depleted
%          across stages); the activity coefficient of every ion pair is
%          approximated by the lumped NaCl-equivalent gamma_b (Sec. 4
%          Pitzer-lite reduction); Mg(OH)2 precipitation assumes alkali
%          dosing holds pH ~ 10.5 (m_OH ~ 3.16e-4 mol/kg).
% Reference: model Sec. 4.2, Eqs. (4.4)-(4.7).

names = par.ions.name;
idxCa = find(strcmpi(names,'Ca')); idxMg = find(strcmpi(names,'Mg'));
idxNa = find(strcmpi(names,'Na')); idxK  = find(strcmpi(names,'K'));
idxSO4= find(strcmpi(names,'SO4'));idxCl = find(strcmpi(names,'Cl'));

m_Ca = stream.m_i_b(:,idxCa); m_Mg = stream.m_i_b(:,idxMg);
m_Na = stream.m_i_b(:,idxNa); m_K  = stream.m_i_b(:,idxK);
m_SO4= stream.m_i_b(:,idxSO4); m_Cl = stream.m_i_b(:,idxCl);
g = stream.gamma_b;
conv = stream.rho_b ./ (1+stream.S_b/1000);   % mol/kg-water -> mol/m^3 conversion factor

%% Eq. 4.4 saturation indices and equilibrium key-ion concentrations
% Stage 1: CaSO4.2H2O (gypsum)
IAP1 = (g.*m_Ca).*(g.*m_SO4);
SI1  = log10(IAP1./par.m1.Ksp_CaSO4);
m_Ca_eq1 = par.m1.Ksp_CaSO4 ./ max(g.^2.*m_SO4, eps);

% Stage 2: Mg(OH)2 (brucite), alkali dosing to pH ~ 10.5
m_OH = 3.16e-4;
m_Mg_eq2 = par.m1.Ksp_MgOH2 ./ max(g.^3.*m_OH.^2, eps);

% Stage 3: NaCl (halite)
m_Na_eq3 = par.m1.Ksp_NaCl ./ max(g.^2.*m_Cl, eps);
SI3 = log10((g.*m_Na).*(g.*m_Cl) ./ par.m1.Ksp_NaCl);

% Stage 4: KCl (sylvite)
m_K_eq4 = par.m1.Ksp_KCl ./ max(g.^2.*m_Cl, eps);
SI4 = log10((g.*m_K).*(g.*m_Cl) ./ par.m1.Ksp_KCl);

%% Eq. 4.5 recovered mass rates (kg/s), floored at zero (precipitation only if SI>0)
c_Ca = m_Ca.*conv; c_Ca_eq1 = m_Ca_eq1.*conv;
c_Mg = m_Mg.*conv; c_Mg_eq2 = m_Mg_eq2.*conv;
c_Na = m_Na.*conv; c_Na_eq3 = m_Na_eq3.*conv;
c_K  = m_K.*conv;  c_K_eq4  = m_K_eq4.*conv;

eta = par.m1.eta_rec; Qb_rho = stream.Q_b; % Q_b*rho_b folded via molar->mass below

mdot_CaSO4 = max(0, stream.Q_b.*max(c_Ca-c_Ca_eq1,0).*par.m1.M_CaSO4.*eta);
mdot_MgOH2 = max(0, stream.Q_b.*max(c_Mg-c_Mg_eq2,0).*par.m1.M_MgOH2.*eta);
mdot_NaCl  = max(0, stream.Q_b.*max(c_Na-c_Na_eq3,0).*par.m1.M_NaCl.*eta);
mdot_KCl   = max(0, stream.Q_b.*max(c_K -c_K_eq4 ,0).*par.m1.M_KCl.*eta);

mdot_solid = [mdot_CaSO4, mdot_MgOH2, mdot_NaCl, mdot_KCl];

%% Eq. 4.6 reagent demand.
% Gypsum (CaSO4.2H2O) precipitates by concentration/SI alone -> NO reagent.
% Mg(OH)2 is precipitated with lime Ca(OH)2 (industrial practice, 1 mol
% Ca(OH)2 per mol Mg(OH)2): Ca(OH)2 + Mg2+ -> Mg(OH)2 + Ca2+.
mdot_lime = par.m1.f_ex .* (1) .* (par.m1.M_CaOH2/par.m1.M_MgOH2) .* mdot_MgOH2;  % kg/s
mdot_reagent = [mdot_lime, zeros(size(mdot_lime))];   % col1 = lime (Mg stage), col2 = none

%% Eq. 4.7 stage specific energy -> process power (kW)
SEC_M1 = par.m1.e_mix + par.m1.e_sep + par.m1.e_dry;   % kWh/kg, identical across stages (assumption)
P_M1 = sum(mdot_solid,2) .* SEC_M1 .* 3600;             % kW

out.mdot_solid = mdot_solid; out.mdot_reagent = mdot_reagent; out.P_M1 = P_M1;
out.SI = [SI1, nan(size(SI1)), SI3, SI4];
out.labels = {'CaSO4_2H2O','Mg_OH_2','NaCl','KCl'};
end
