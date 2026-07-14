function vr = validate_results(inv, R, mc, par)
%VALIDATE_RESULTS Automated verification and validation gate for the
%brine-valorization model outputs.
%
% Purpose: Runs a battery of physical-consistency, conservation, second-
%          law and numerical-sanity checks across every result struct and
%          returns a pass/fail report. The user requirement is that all
%          produced analysis data be validated and verified; this routine
%          is the single auditable gate that enforces it.
% Inputs : inv - inventory table; R - run_pipeline output; mc - monte_carlo
%          output (may be []); par - parameter struct.
% Outputs: vr (struct): .checks (table-like struct array name/pass/detail),
%          .n_pass, .n_fail, .all_pass (logical).
% Reference: model Sec. 2.5 / 6 closure requirements; 2nd-law admissibility.

c = {};   % each entry: {name, pass(logical), detail(string)}

add = @(nm,ok,dt) {nm, logical(ok), string(dt)};

% --- mass balance closure (Eq. 2.1-2.3) ---
mc_err = max(R.massE.mass_closure_err);
c{end+1} = add('Mass balance closure Q_f=Q_p+Q_b', mc_err < 1e-9, sprintf('max rel err = %.2e', mc_err));

% --- salinity / density / osmotic ordering (Eq. 2.6-2.13) ---
c{end+1} = add('Brine salinity > feed salinity', all(R.stream.S_b > R.stream.S_f), 'all plants');
c{end+1} = add('Brine density > feed density',  all(R.stream.rho_b > R.stream.rho_f), 'all plants');
c{end+1} = add('Osmotic pressure positive & ordered', all(R.stream.pi_b > R.stream.pi_f & R.stream.pi_f>0), 'all plants');

% --- no NaN/Inf in headline outputs ---
hd = [R.econ.LCOW_net; R.econ.NPV; R.econ.dLCOW; R.agg.VI; R.lca.CO2_mitigation_t_yr];
c{end+1} = add('Headline outputs finite (no NaN/Inf)', all(isfinite(hd)), sprintf('%d values', numel(hd)));

% --- second-law admissibility: eta_II <= 1 for all pathways ---
e = R.exergy;
okII = all(e.eta_II_RED<=1+1e-9) && all(e.eta_II_PRO<=1+1e-9) && ...
       all(e.eta_II_ERD<=1+1e-9) && all(e.eta_II_th<=1+1e-9);
c{end+1} = add('Second-law efficiencies <= 1', okII, ...
    sprintf('max etaII RED/PRO/ERD/th = %.2f/%.2f/%.2f/%.2f', ...
    max(e.eta_II_RED),max(e.eta_II_PRO),max(e.eta_II_ERD),max(e.eta_II_th)));

% --- RED/PRO net power <= eta_II_max * Gibbs ceiling (Eq. 3.11) ---
okRED = all(R.energy.red.P_net <= par.red.eta_II_max.*R.energy.red.dG_mix_rate + 1e-6);
c{end+1} = add('RED net power within thermodynamic ceiling', okRED, 'P_net <= 0.4*dG_mix');

% --- exergy destruction non-negative (Gouy-Stodola) ---
c{end+1} = add('Exergy destruction non-negative', all(e.Exdot_dest(:) >= -1e-6), ...
    sprintf('min = %.2e W', min(e.Exdot_dest(:))));

% --- energy balance sign: recovered energy non-negative ---
c{end+1} = add('Recovered energy non-negative', all(R.energy.E_rec_yr >= 0), 'all plants');

% --- portfolio selection consistency: selected => marginal NPV>0 ---
selOK = all(all( (R.econ.sel==0) | (R.econ.NPV_k > -1e-6) ));
c{end+1} = add('Portfolio selection consistent (sel => NPV_k>0)', selOK, 'all pathways');

% --- market-absorption fractions within (0,1] ---
c{end+1} = add('Market-absorption fractions in (0,1]', all(R.agg.f_absorb>0 & R.agg.f_absorb<=1+1e-9), ...
    sprintf('min f_abs = %.3f', min(R.agg.f_absorb)));

% --- economic identity: dLCOW = (ACC - CF)/Q_p ---
lhs = R.econ.dLCOW; rhs = (R.econ.ACC_total - R.econ.CF_net)./R.econ.Q_p_m3yr;
c{end+1} = add('Net-LCOW identity holds', max(abs(lhs-rhs)) < 1e-6, ...
    sprintf('max dev = %.2e', max(abs(lhs-rhs))));

% --- viability index bounded [0,1] ---
c{end+1} = add('Viability index in [0,1]', all(R.agg.VI>=0 & R.agg.VI<=1), ...
    sprintf('range [%.2f, %.2f]', min(R.agg.VI), max(R.agg.VI)));

% --- Monte Carlo (if supplied): outputs finite, P5<=P50<=P95 ---
if ~isempty(mc)
    okMC = true; dt = 'P5<=P50<=P95 for all KPIs';
    for m = 1:numel(mc.kpi_names)
        s = mc.stats.(mc.kpi_names{m});
        if ~(s.P5<=s.P50+1e-9 && s.P50<=s.P95+1e-9 && isfinite(s.mean)); okMC=false; dt=mc.kpi_names{m}; break; end
    end
    c{end+1} = add('Monte-Carlo percentile ordering', okMC, dt);
end

% --- assemble ---
names = cellfun(@(x)x{1}, c, 'uni', 0)';
pass  = cellfun(@(x)x{2}, c)';
detail= cellfun(@(x)x{3}, c, 'uni', 0)';
vr.checks = struct('name',names,'pass',num2cell(pass),'detail',detail);
vr.n_pass = sum(pass); vr.n_fail = sum(~pass); vr.all_pass = all(pass);

fprintf('\n================ VALIDATION REPORT ================\n');
for ii = 1:numel(names)
    if pass(ii); tag = 'PASS'; else; tag = 'FAIL'; end
    fprintf('  [%s] %-46s | %s\n', tag, names{ii}, detail{ii});
end
fprintf('--------------------------------------------------\n');
if vr.all_pass; overall = 'ALL PASS'; else; overall = 'FAILURES PRESENT'; end
fprintf('  %d/%d checks passed.  Overall: %s\n', vr.n_pass, numel(names), overall);
fprintf('==================================================\n\n');
end
