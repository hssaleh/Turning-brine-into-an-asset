function iv = validate_inputs(par)
%VALIDATE_INPUTS Audit, validation and verification of all model INPUT data.
%
% Purpose: Independently checks every parameter in the configuration struct
%          against physical admissibility bounds and published literature
%          ranges BEFORE the model runs, so that no result is produced from
%          an out-of-range or internally-inconsistent input. Complements
%          validate_results.m (which checks the OUTPUTS) and
%          literature_benchmarks.m (which checks computed quantities).
% Inputs : par - struct from config_parameters().
% Outputs: iv (struct): .checks (name, status PASS/WARN/FAIL, value, range,
%          detail), .n_pass, .n_warn, .n_fail, .all_ok (no FAIL).
% References: Sharqawy et al. (2010); RED/PRO/ZLD/DLE techno-economic
%          literature; standard seawater composition (Millero 2008).

C = {};   % {name, status, valuestr, rangestr}
function add(nm, ok, warnonly, val, rng)
    if ok, st='PASS'; elseif warnonly, st='WARN'; else, st='FAIL'; end
    C{end+1} = {nm, st, val, rng}; %#ok<AGROW>
end
inb  = @(x,lo,hi) all(x>=lo & x<=hi);
s = @(x) num2str(x,'%.4g');

%% --- fractions / efficiencies must lie in (0,1] ---
effs = { 'red.eta_pump',par.red.eta_pump; 'red.eta_DC',par.red.eta_DC; ...
    'red.epsilon_sp',par.red.epsilon_sp; 'red.alpha',par.red.alpha; ...
    'red.eta_II_max',par.red.eta_II_max; 'pro.eta_turb',par.pro.eta_turb; ...
    'pro.eta_ERD',par.pro.eta_ERD; 'pro.eta_II_max',par.pro.eta_II_max; ...
    'erd.eta_ERD',par.erd.eta_ERD; 'ec.eta_F',par.ec.eta_F; ...
    'thermal.eta_II_ORC',par.thermal.eta_II_ORC; 'm1.eta_rec',par.m1.eta_rec; ...
    'm2.eta_NaCl',par.m2.eta_NaCl; 'm2.eta_KCl',par.m2.eta_KCl; 'm2.f_heat',par.m2.f_heat; ...
    'm3.eta_Li',par.m3.eta_Li; 'm4.eta_Br',par.m4.eta_Br; 'm4.eta_Rb',par.m4.eta_Rb; ...
    'm5.r_ZLD',par.m5.r_ZLD; 'stream.f_av',par.stream.f_av; 'stream.r',par.stream.r; ...
    'stream.R_salt',par.stream.R_salt };
for k=1:size(effs,1)
    v=effs{k,2}; add(['eff: ' effs{k,1}], inb(v,0,1) && v>0, false, s(v), '(0, 1]');
end

%% --- prices strictly positive ---
prices = { 'p_CaSO4',par.m1.p_CaSO4; 'p_MgOH2',par.m1.p_MgOH2; 'p_NaCl',par.m1.p_NaCl; ...
    'p_KCl',par.m1.p_KCl; 'p_lime',par.m1.p_lime; 'p_Li2CO3',par.m3.p_Li2CO3; ...
    'p_Br',par.m4.p_Br; 'p_Rb',par.m4.p_Rb; 'p_water',par.m5.p_water; ...
    'p_elec',par.econ.p_elec; 'c_mem',par.econ.c_mem };
for k=1:size(prices,1)
    v=prices{k,2}; add(['price: ' prices{k,1}], v>0, false, s(v), '> 0');
end

%% --- literature-range checks (WARN if outside, not fatal) ---
add('RED stack cost (USD/m2)', inb(par.econ.cap.RED.C_ref/par.econ.cap.RED.Size_ref*par.econ.CEPCI_ratio*par.econ.f_install,80,400), true, ...
    s(par.econ.cap.RED.C_ref/par.econ.cap.RED.Size_ref*par.econ.CEPCI_ratio*par.econ.f_install), '80-400');
add('PRO stack cost (USD/m2)', inb(par.econ.cap.PRO.C_ref/par.econ.cap.PRO.Size_ref*par.econ.CEPCI_ratio*par.econ.f_install,80,500), true, ...
    s(par.econ.cap.PRO.C_ref/par.econ.cap.PRO.Size_ref*par.econ.CEPCI_ratio*par.econ.f_install), '80-500');
add('SWRO SEC_des (kWh/m3)', inb(par.tech.SEC_des.SWRO,2.5,4.5), true, s(par.tech.SEC_des.SWRO), '2.5-4.5');
add('MED SEC_des (kWh/m3 th-eq)', inb(par.tech.SEC_des.MED,6,16), true, s(par.tech.SEC_des.MED), '6-16');
add('MSF SEC_des (kWh/m3 th-eq)', inb(par.tech.SEC_des.MSF,10,25), true, s(par.tech.SEC_des.MSF), '10-25');
add('ZLD SEC (kWh/m3)', inb(par.m5.SEC_ZLD,7,25), true, s(par.m5.SEC_ZLD), '7-25');
add('Li SEC (kWh/kg Li2CO3)', inb(par.m3.SEC_Li,10,80), true, s(par.m3.SEC_Li), '10-80');
add('Li q_max (mg/g)', inb(par.m3.q_max,15,45), true, s(par.m3.q_max), '15-45');
add('discount rate i', inb(par.econ.i,0.03,0.15), true, s(par.econ.i), '0.03-0.15');
add('project life N (yr)', inb(par.econ.N,15,40), true, s(par.econ.N), '15-40');
add('six-tenths exponent n_u', inb(par.econ.n_u,0.5,0.8), true, s(par.econ.n_u), '0.5-0.8');
add('install factor f_install', inb(par.econ.f_install,1.1,2.0), true, s(par.econ.f_install), '1.1-2.0');
add('grid EF (kgCO2/kWh)', inb(par.lca.EF_grid,0.0,1.0), true, s(par.lca.EF_grid), '0-1');
add('feed salinity S_f (g/kg)', inb(par.stream.S_f,30,45), true, s(par.stream.S_f), '30-45');
add('SGE dilute salinity (g/kg)', inb(par.sge.S_dilute,0.5,15), true, s(par.sge.S_dilute), '0.5-15');
add('RED membrane permselectivity', inb(par.red.alpha,0.8,1.0), true, s(par.red.alpha), '0.8-1.0');
add('halite sat. molality (mol/kg)', inb(par.m2.m_sat_NaCl,5.5,6.5), true, s(par.m2.m_sat_NaCl), '5.5-6.5');

%% --- internal consistency / cross-checks ---
sumw = sum(par.ions.w35);
add('ion table closure (sum w35 ~ 35)', abs(sumw-35)<1.0, false, s(sumw), '35 +/- 1 g/kg');
add('ion arrays same length', numel(par.ions.name)==numel(par.ions.w35) && numel(par.ions.w35)==numel(par.ions.M) && numel(par.ions.M)==numel(par.ions.z), false, ...
    s(numel(par.ions.name)), 'equal');
add('electroneutrality (|charge imbalance|)', abs(sum((par.ions.w35./1000)./par.ions.M.*par.ions.z))<0.05, true, ...
    s(sum((par.ions.w35./1000)./par.ions.M.*par.ions.z)), '~0 eq/kg');
add('T_f = T0 (dead-state consistency)', abs(par.stream.T_f-par.stream.T0)<1e-9, true, s(par.stream.T_f-par.stream.T0), '0 K');
add('brine molality < halite saturation', (par.stream.S_f/(1-par.stream.r))/(1000)/par.const.M_NaCl < par.m2.m_sat_NaCl, true, ...
    s((par.stream.S_f/(1-par.stream.r))/(1000)/par.const.M_NaCl), ['< ' s(par.m2.m_sat_NaCl)]);
add('MC draws >= 100', par.unc.N_MC>=100, true, s(par.unc.N_MC), '>= 100');
add('uncertain P5 < P95 for all vars', all(structfun(@(z) z.P5<z.P95, par.unc.vars)), false, 'ok', 'P5<P95');
add('market shares in (0,1]', all(struct2array_local(par.market.addressable_share)>0 & struct2array_local(par.market.addressable_share)<=1), false, 'ok', '(0,1]');
add('viability weights sum to 1', abs(par.viability.w_fin+par.viability.w_thermo+par.viability.w_cost+par.viability.w_env-1)<1e-9, false, ...
    s(par.viability.w_fin+par.viability.w_thermo+par.viability.w_cost+par.viability.w_env), '= 1');

%% --- assemble & report ---
names=cellfun(@(x)x{1},C,'uni',0)'; stat=cellfun(@(x)x{2},C,'uni',0)';
valc=cellfun(@(x)x{3},C,'uni',0)'; rngc=cellfun(@(x)x{4},C,'uni',0)';
iv.checks = table(names,stat,valc,rngc,'VariableNames',{'parameter','status','value','range'});
iv.n_pass=sum(strcmp(stat,'PASS')); iv.n_warn=sum(strcmp(stat,'WARN')); iv.n_fail=sum(strcmp(stat,'FAIL'));
iv.all_ok = iv.n_fail==0;

fprintf('\n================== INPUT DATA VALIDATION ==================\n');
for k=1:numel(names)
    if ~strcmp(stat{k},'PASS')
        fprintf('  [%s] %-36s value=%-10s range=%s\n', stat{k}, names{k}, valc{k}, rngc{k});
    end
end
fprintf('  (%d PASS not listed)\n', iv.n_pass);
fprintf('----------------------------------------------------------\n');
if iv.all_ok, verdict='ADMISSIBLE'; else, verdict='HAVE FAILURES'; end
fprintf('  %d checks: %d PASS, %d WARN, %d FAIL.  Inputs %s\n', numel(names), iv.n_pass, iv.n_warn, iv.n_fail, verdict);
fprintf('==========================================================\n\n');
end

% ----------------------------------------------------------------------
function a = struct2array_local(s)
f = fieldnames(s); a = zeros(1,numel(f));
for k=1:numel(f), a(k)=s.(f{k}); end
end
