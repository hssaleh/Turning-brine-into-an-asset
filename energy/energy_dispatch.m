function energy = energy_dispatch(stream, par)
%ENERGY_DISPATCH Dispatcher selecting admissible energy-recovery
%pathways by plant technology and aggregating their outputs.
%
% Purpose: Implements the Sec. 3 dispatch logic and Eqs. (6.1)-(6.2):
%          RED and PRO require a low-salinity counter-stream (assumed
%          available, fed by the plant's own intake seawater); ERD
%          requires residual brine pressure (SWRO only); the thermal
%          route requires a temperature lift (MED/MSF only); the
%          electrochemical route is always scored as a conversion
%          option. Net recovered energy excludes the EC pathway, per
%          Sec. 3.5 ("not net primary energy harvested from the
%          gradient").
% Inputs : stream - struct from M00_streams.m
%          par    - parameter struct
% Outputs: energy (struct) with per-pathway sub-structs (red, pro, erd,
%          thermal, ec) and aggregate fields:
%             Edot_rec [W], E_rec_yr [kWh/yr], E_EC_yr [kWh/yr]
% Reference: model Sec. 3, Eqs. (6.1)-(6.2).

N = numel(stream.Q_b);
isSWRO = strcmpi(stream.tech, 'SWRO');
isThermalTech = ~isSWRO;

red = red_power(stream, par);
pro = pro_power(stream, par);
erd = erd_power(stream, par);
thermal = thermal_power(stream, par);
ec = electrochem_value(stream, par);

% RED/PRO admissible for all plants (intake seawater as dilute stream)
% ERD admissible only where residual pressure exists (SWRO)
erd.P_net(~isSWRO) = 0; erd.P_gross(~isSWRO) = 0;
% Thermal admissible only where a temperature lift exists (MED/MSF handled inside thermal_power via T_b<=T0)

Edot_rec = red.P_net + pro.P_net + erd.P_net + thermal.P_ORC_net;   % Eq. 6.1
Edot_rec = max(Edot_rec, 0);

f_av = par.stream.f_av;
E_rec_yr = Edot_rec .* 8760 .* f_av / 1000;        % Eq. 6.2, W->kW, kWh/yr
E_EC_yr  = ec.P_demand .* 8760 .* f_av / 1000;     % electrochemical electrical demand, kWh/yr

energy.red = red; energy.pro = pro; energy.erd = erd; energy.thermal = thermal; energy.ec = ec;
energy.Edot_rec = Edot_rec; energy.E_rec_yr = E_rec_yr; energy.E_EC_yr = E_EC_yr;
end
