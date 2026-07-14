function out = trace_recovery(stream, par)
%TRACE_RECOVERY (M4) Trace-element recovery (bromine, rubidium).
%
% Purpose: Implements Eq. (4.14): conservative concentration accounting
%          with element-specific recovery efficiency.
% Inputs : stream - struct from M00_streams.m (Q_b, rho_b, S_b)
%          par    - parameter struct (par.m4.*, par.ions.*)
% Outputs: out.mdot_Br, out.mdot_Rb [kg/s]
% Reference: model Sec. 4.5, Eq. (4.14).

names = par.ions.name;
w35_Br = par.ions.w35(strcmpi(names,'Br'));
w35_Rb = par.ions.w35(strcmpi(names,'Rb'));

w_Br_b = w35_Br .* (stream.S_b/35) / 1000;   % kg/kg brine
w_Rb_b = w35_Rb .* (stream.S_b/35) / 1000;   % kg/kg brine

mdot_Br = par.m4.eta_Br .* stream.Q_b .* stream.rho_b .* w_Br_b;
mdot_Rb = par.m4.eta_Rb .* stream.Q_b .* stream.rho_b .* w_Rb_b;

out.mdot_Br = mdot_Br; out.mdot_Rb = mdot_Rb;
end
