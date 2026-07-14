function export_results_mat(inv, R, mc, vr, lb, iv, par, matpath)
%EXPORT_RESULTS_MAT Persist the full reproducible result set to a MAT file
%and a compact CSV of headline per-plant results.
%
% Purpose: Saves every input, intermediate and output struct so that all
%          figures, tables and statistics can be regenerated without
%          re-running the model, satisfying the reproducibility
%          requirement.
% Inputs : inv,R,mc,vr,par - model artefacts; matpath - .mat output path.
% Outputs: none (writes <matpath> and a sibling _plants.csv).
% Reference: project result-export requirements (CLAUDE.md).

results.inventory = inv;
results.pipeline  = R;
results.montecarlo = mc;
results.validation = vr;
results.lit_benchmarks = lb;
results.input_validation = iv;
results.parameters = par;
results.timestamp  = datetime('now');
save(matpath, 'results', '-v7.3');

% compact per-plant CSV of headline results
csv = fullfile(fileparts(matpath), 'BrineValorization_plants.csv');
T = table(inv.id, inv.region, inv.tech, inv.S_f, R.stream.S_b, inv.r, ...
    R.econ.CAPEX_total, R.econ.OPEX_total, R.econ.Rev_total, R.econ.NPV, ...
    R.econ.IRR, R.econ.payback, R.econ.BCR, R.econ.LCOW_net, R.econ.dLCOW, ...
    R.lca.CO2_mitigation_t_yr, R.lca.EROI, R.agg.VI, ...
    'VariableNames', {'id','region','tech','S_f','S_b','r','CAPEX','OPEX','Revenue', ...
    'NPV','IRR','payback_yr','BCR','LCOW_net','dLCOW','CO2_mit_tyr','EROI','VI'});
writetable(T, csv);
fprintf('MAT export complete: %s\nCSV export complete: %s\n', matpath, csv);
end
