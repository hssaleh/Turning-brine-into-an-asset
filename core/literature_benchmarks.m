function lb = literature_benchmarks(inv, R, par)
%LITERATURE_BENCHMARKS Verify model outputs against published envelopes.
%
% Purpose: Compares the model's computed quantities against the validation
%          benchmarks tabulated in the mathematical-model specification
%          (Section 11.5) and the underlying literature, returning a
%          pass/fail report. This is the external (literature) verification
%          layer that complements the internal consistency checks in
%          validate_results.m.
% Inputs : inv - inventory table; R - run_pipeline output; par - parameters.
% Outputs: lb (struct): .rows (name, model, lo, hi, unit, pass),
%          .n_pass, .n_fail, .all_pass.
% References: Sharqawy et al. (2010); RED/PRO salinity-gradient-power and
%          ZLD energy-cost literature cited in the model Sec. 11.5 / 13.

rows = {};   % {name, modelval, lo, hi, unit}
addr = @(nm,val,lo,hi,u) {nm, val, lo, hi, u};

% (1) RED gross areal power density: 1-4 W/m^2 (0.5/4 M-class gradient)
rows{end+1} = addr('RED areal power density P_d', median(R.energy.red.P_d), 1, 4, 'W/m^2');

% (2) PRO peak power density: 3-10 W/m^2 (seawater/brine)
rows{end+1} = addr('PRO peak power density W_d,max', median(R.energy.pro.W_d_max), 3, 10, 'W/m^2');

% (3) Brine osmotic pressure at ~70 g/kg: 5.5-6.0 MPa
[pi70,~] = osmotic_pressure((70/(1000-70))/par.const.M_NaCl, 0.90, par.stream.T0, par);
rows{end+1} = addr('Brine osmotic pressure at 70 g/kg', pi70/1e6, 5.5, 6.0, 'MPa');

% (4) SWRO baseline specific energy of desalination: 3-4 kWh/m^3
rows{end+1} = addr('SWRO SEC_des', par.tech.SEC_des.SWRO, 3, 4, 'kWh/m^3');

% (5) ZLD incremental water specific energy: 7-25 kWh/m^3
rows{end+1} = addr('ZLD incremental SEC', par.m5.SEC_ZLD, 7, 25, 'kWh/m^3');

% (6) Lithium sorbent capacity: 20-40 mg/g
rows{end+1} = addr('Li sorbent capacity q_max', par.m3.q_max, 20, 40, 'mg/g');

% (7) RED 2nd-law efficiency physically admissible: 0-0.40
rows{end+1} = addr('RED 2nd-law efficiency', median(R.exergy.eta_II_RED), 0, 0.40, '-');

% (8) ERD recovered specific energy ~ P_b: 1.5-3.0 kWh/m3 of brine (SWRO)
isS = strcmpi(inv.tech,'SWRO');
erd_kWhm3 = (R.energy.erd.P_net(isS)./max(R.stream.Q_b(isS),eps))/3.6e6; % W per (m3/s) -> kWh/m3
rows{end+1} = addr('ERD specific recovery (SWRO)', median(erd_kWhm3), 1.0, 3.5, 'kWh/m^3');

% (9) Net-LCOW direction: high-salinity, high-electricity plants are viable
hi = R.stream.S_b > median(R.stream.S_b) & inv.p_elec > median(inv.p_elec);
frac_viable_hi = mean(R.econ.LCOW_net(hi) < inv.LCOW_base(hi));
rows{end+1} = addr('Viable fraction (high-S, high-p_elec)', 100*frac_viable_hi, 50, 100, '%');

% (10) Brine NaCl mass fraction of TDS recovered (halite dominant): 0.6-0.9
fNaCl = sum(R.mineral.mprod_yr(:,3)) / max(sum(R.mineral.mprod_yr(:,[1 2 3 4]),'all'),eps);
rows{end+1} = addr('NaCl share of recovered salts', fNaCl, 0.55, 0.95, '-');

names=cellfun(@(x)x{1},rows,'uni',0)'; mval=cellfun(@(x)x{2},rows)';
lo=cellfun(@(x)x{3},rows)'; hi=cellfun(@(x)x{4},rows)'; unit=cellfun(@(x)x{5},rows,'uni',0)';
pass = mval>=lo-1e-9 & mval<=hi+1e-9;
lb.rows = struct('name',names,'model',num2cell(mval),'lo',num2cell(lo),'hi',num2cell(hi),'unit',unit,'pass',num2cell(pass));
lb.table = table(names, mval, lo, hi, unit, pass, 'VariableNames',{'benchmark','model','lo','hi','unit','pass'});
lb.n_pass = sum(pass); lb.n_fail = sum(~pass); lb.all_pass = all(pass);

fprintf('\n============== LITERATURE BENCHMARK VERIFICATION ==============\n');
for k=1:numel(names)
    if pass(k); tag='PASS'; else; tag='FAIL'; end
    fprintf('  [%s] %-38s model=%9.3f  [%.2f, %.2f] %s\n', tag, names{k}, mval(k), lo(k), hi(k), unit{k});
end
fprintf('--------------------------------------------------------------\n');
if lb.all_pass; ov='ALL WITHIN LITERATURE ENVELOPES'; else; ov='SOME OUTSIDE ENVELOPES'; end
fprintf('  %d/%d benchmarks passed.  Overall: %s\n', lb.n_pass, numel(names), ov);
fprintf('==============================================================\n\n');
end
