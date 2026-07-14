%% MAIN_DRIVER  Turning Brine into an Asset - end-to-end model orchestration
%
% Project   : Global Techno-Economic, Exergetic and Environmental Assessment
%             of Energy and Mineral Recovery from Desalination Brine
% Purpose   : Single entry point that builds the parameter set and plant
%             inventory, runs the full deterministic pipeline, propagates
%             uncertainty (Latin-hypercube Monte Carlo + regression-based
%             global sensitivity), validates every result, generates the
%             complete publication-quality figure suite and exports all
%             numeric data to Excel (one sheet per figure), MAT and CSV.
% Author    : <Hossam AbdelMeguid>   Date: <2026>   Version: 2.0
% Inputs    : none (all constants come from config_parameters.m).
% Outputs   : results/BrineValorization_Figures.xlsx, results/*.mat/.csv,
%             results/figures/*.fig, console validation + headline report.
% Reference : Brine_Valorization_Mathematical_Model.md, Sections 1-10;
%             Brine_Valorization_Paper_Design.docx; Nature concepts doc.

clear; clc; close all;
root = fileparts(mfilename('fullpath'));
addpath(genpath(root));
resdir = fullfile(root,'results'); if ~exist(resdir,'dir'), mkdir(resdir); end
figdir = fullfile(resdir,'figures'); if ~exist(figdir,'dir'), mkdir(figdir); end

fprintf('==== Turning Brine into an Asset : model run ====\n');

%% 1. Parameters and inventory (with input-data validation/verification)
par = config_parameters();
iv  = validate_inputs(par);
if ~iv.all_ok
    warning('Input-data validation reported failures - inspect the report above before running.');
end
inv = generate_inventory(par);
fprintf('Inventory built: %d plants across %d regions.\n', height(inv), numel(unique(inv.region)));

%% 2. Deterministic pipeline (Sections 2-9)
R = run_pipeline(inv, par);
fprintf('Deterministic pipeline complete.\n');

%% 3. Uncertainty propagation and global sensitivity (Section 10)
mc = monte_carlo(inv, par);
fprintf('Monte-Carlo ensemble complete (N=%d draws).\n', mc.N);

%% 4. Validation and verification gate (internal consistency + literature)
vr = validate_results(inv, R, mc, par);
lb = literature_benchmarks(inv, R, par);
if ~vr.all_pass
    warning('Internal validation reported failures - inspect the report above.');
end
if ~lb.all_pass
    warning('Some outputs fall outside literature envelopes - inspect the benchmark report.');
end

%% 5. Headline results report
k = fleet_kpis(inv, R);
fprintf('\n---------------- HEADLINE RESULTS ----------------\n');
fprintf('  Global recovered energy ........ %7.1f TWh/yr\n', k.E_rec_TWh);
fprintf('  Global CO2 mitigation .......... %7.1f Mt/yr\n', k.CO2_mit_Mt);
fprintf('  Global merchant NPV ............ %7.0f  B USD\n', k.NPV_global_B);
fprintf('  Integration premium ............ %7.0f  B USD\n', k.integ_prem_B);
fprintf('  Fleet-mean net-LCOW ............ %7.2f  USD/m3\n', k.LCOW_net_mean);
fprintf('  Plants with NPV>0 .............. %7.0f %%\n', 100*k.frac_NPV_pos);
fprintf('  Plants with net-LCOW<baseline .. %7.0f %%\n', 100*k.frac_LCOW_down);
fprintf('  Median IRR / payback ........... %5.0f %% / %.1f yr\n', 100*median(R.econ.IRR,'omitnan'), median(R.econ.payback(isfinite(R.econ.payback))));
fprintf('  Mg(OH)2 absorbable ............. %7.2f Mt/yr (of %.1f technical)\n', R.agg.global.mineral_Mt_yr_constr(2), R.agg.global.mineral_Mt_yr_tech(2));
fprintf('--------------------------------------------------\n\n');

%% 6. Figures (grouped tiledlayout, all visible, saved as .fig)
FIGS = make_figures(inv, R, mc, lb, par, figdir);
fprintf('Generated %d figure groups (saved to %s).\n', numel(FIGS), figdir);

%% 7. Exports
xlsx = fullfile(resdir,'BrineValorization_Figures.xlsx');
export_excel(FIGS, xlsx);
% Append the input-validation and literature-benchmark audit tables
writetable(iv.checks, xlsx, 'Sheet', 'Input_Validation');
writetable(lb.table,  xlsx, 'Sheet', 'Literature_Benchmarks');
export_results_mat(inv, R, mc, vr, lb, iv, par, fullfile(resdir,'BrineValorization_Results.mat'));

fprintf('\n==== Run complete. All outputs in %s ====\n', resdir);
