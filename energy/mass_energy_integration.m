function massE = mass_energy_integration(plant, stream, energy, mineral, par)
%MASS_ENERGY_INTEGRATION Plant-level mass and energy balance closure.
%
% Purpose: Implements Eqs. (6.1)-(6.5): primary desalination energy
%          demand, recovered/valorization energy roll-up, net energy
%          impact, energy-offset fraction and the mass-balance closure
%          check (Q_f = Q_p + Q_b) used as a model-validation gate.
% Inputs : plant   - inventory table (C_des, tech)
%          stream  - struct from M00_streams.m
%          energy  - struct from energy_dispatch.m
%          mineral - struct from mineral_dispatch.m
%          par     - parameter struct
% Outputs: massE (struct):
%             E_dem_yr [kWh/yr] primary desalination demand (Eq. 6.3 LHS)
%             E_val_yr [kWh/yr] total valorization energy use/recovery
%               (= E_rec_yr - E_min_yr - E_ZLD_yr - E_EC_yr, signed: +
%               harvested, - consumed)
%             E_net_yr [kWh/yr] net plant energy impact (Eq. 6.4)
%             f_offset [-] fraction of primary demand offset by net
%               recovered energy (Eq. 6.5), clipped to [-1, 1] for
%               reporting sanity
%             mass_closure_err [-] relative mass-balance closure error
%               |Q_f-(Q_p+Q_b)|/Q_f, must be ~0 (machine precision)
% Assumptions: thermal-equivalent SEC for MED/MSF (par.tech.SEC_des)
%          is used uniformly as the primary energy metric (no separate
%          electrical/thermal exergy weighting at this integration
%          stage; that distinction is carried instead in the exergy
%          module's second-law efficiencies).
% Reference: model Sec. 6, Eqs. (6.1)-(6.5).

f_av = par.stream.f_av;
N = height(plant);

%% Eq. 6.3 primary desalination energy demand
SEC = zeros(N,1);
techs = fieldnames(par.tech.SEC_des);
for k = 1:numel(techs)
    mask = strcmpi(plant.tech, techs{k});
    SEC(mask) = par.tech.SEC_des.(techs{k});
end
Q_p_m3yr = stream.Q_p .* 86400 .* 365 .* f_av;     % m3/yr
E_dem_yr = SEC .* Q_p_m3yr;                          % kWh/yr

%% Eq. 6.4 net plant energy impact
E_val_yr = energy.E_rec_yr - mineral.E_min_yr - mineral.E_ZLD_yr - energy.E_EC_yr;
E_net_yr = E_val_yr - 0;   % no further terms beyond valorization vs. recovery at this stage

%% Eq. 6.5 energy-offset fraction relative to primary desalination demand
f_offset = E_net_yr ./ max(E_dem_yr, eps);
f_offset = max(min(f_offset, 1), -1);

%% Mass-balance closure check (validation gate)
Q_check = stream.Q_p + stream.Q_b;
mass_closure_err = abs(stream.Q_f - Q_check) ./ max(stream.Q_f, eps);
assert(all(mass_closure_err < 1e-9), 'mass_energy_integration: mass-balance closure violated');

massE.E_dem_yr = E_dem_yr; massE.E_val_yr = E_val_yr; massE.E_net_yr = E_net_yr;
massE.f_offset = f_offset; massE.mass_closure_err = mass_closure_err;
end
