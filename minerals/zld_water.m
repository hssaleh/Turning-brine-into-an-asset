function out = zld_water(stream, par)
%ZLD_WATER (M5) Near-zero-liquid-discharge incremental water recovery.
%
% Purpose: Implements Eqs. (4.15)-(4.17): added high-pressure RO/MD/
%          MED-MVC stage lifting overall recovery, incremental
%          freshwater, its specific energy and the residual salinity.
% Inputs : stream - struct from M00_streams.m (Q_b, S_b)
%          par    - parameter struct (par.m5.*)
% Outputs: out.Q_water_add [m3/s], out.P_M5 [kW], out.S_residual [g/kg],
%          out.Q_res [m3/s] residual brine after incremental recovery
% Reference: model Sec. 4.5, Eqs. (4.15)-(4.17).

Q_water_add = stream.Q_b .* par.m5.r_ZLD;
P_M5 = Q_water_add .* par.m5.SEC_ZLD .* 3600;   % kW
S_residual = stream.S_b ./ (1-par.m5.r_ZLD);
Q_res = stream.Q_b - Q_water_add;

out.Q_water_add = Q_water_add; out.P_M5 = P_M5;
out.S_residual = S_residual; out.Q_res = Q_res;
end
