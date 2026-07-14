function R = run_pipeline(plant, par)
%RUN_PIPELINE Single deterministic evaluation of the full brine-
%valorization model for a given plant inventory and parameter set.
%
% Purpose: Chains every physics/economics/environment module in
%          dependency order and returns all per-plant result structs in
%          a single container. This is the deterministic core re-used by
%          both main_driver.m (once) and monte_carlo.m (N_MC times under
%          perturbed parameters), guaranteeing identical model logic in
%          the point estimate and the uncertainty propagation.
% Inputs : plant - inventory table from generate_inventory.m
%          par   - parameter struct (possibly Monte-Carlo-perturbed)
% Outputs: R (struct) with fields stream, energy, mineral, exergy,
%          massE, econ, lca, agg (each the output of its module).
% Reference: model Sections 2-9 (full integration).

R.stream  = M00_streams(plant, par);
R.energy  = energy_dispatch(R.stream, par);
R.mineral = mineral_dispatch(R.stream, par);
R.exergy  = exergy_balance(R.stream, R.energy, R.mineral, par);
R.massE   = mass_energy_integration(plant, R.stream, R.energy, R.mineral, par);
R.econ    = techno_economics(plant, R.stream, R.energy, R.mineral, par);
R.lca     = lca_environmental(plant, R.stream, R.energy, R.mineral, R.massE, R.econ, par);
R.agg     = global_aggregate(plant, R.energy, R.mineral, R.econ, R.lca, par);
end
