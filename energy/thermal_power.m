function out = thermal_power(stream, par)
%THERMAL_POWER Low-grade-heat recovery pathway (ORC and MD alternative).
%
% Purpose: Implements Eqs. (3.21)-(3.24): sensible thermal availability,
%          Carnot-limited and ORC net work, membrane-distillation
%          distillate yield (reported as an alternative, mutually-
%          exclusive use of the same heat).
% Inputs : stream - struct from M00_streams.m (Q_b, rho_b, cp_b, T0, tech)
%          par    - parameter struct (par.thermal.*)
% Outputs: out.Qdot_th, out.W_Carnot, out.P_ORC_net [W]; out.mdot_dist [kg/s]
% Assumptions: brine temperature T_b set by technology (SWRO 25C => no
%          lift => zero thermal output; MED/MSF 40C); h_fg evaluated at
%          40C for all thermal plants (engineering simplification).
% Reference: model Sec. 3.4, Eqs. (3.21)-(3.24).

N = numel(stream.Q_b);
T_b = zeros(N,1);
techs = fieldnames(par.thermal.Tb_by_tech);
for k = 1:numel(techs)
    mask = strcmpi(stream.tech, techs{k});
    T_b(mask) = par.thermal.Tb_by_tech.(techs{k});
end

mdot_b = stream.Q_b .* stream.rho_b;
Qdot_th = mdot_b .* stream.cp_b .* max(T_b - stream.T0, 0);

W_Carnot = Qdot_th .* (1 - stream.T0./max(T_b,stream.T0+1e-6));
eta_ORC  = par.thermal.eta_II_ORC .* (1 - stream.T0./max(T_b,stream.T0+1e-6));
P_ORC_net = eta_ORC .* Qdot_th;
P_ORC_net(T_b <= stream.T0) = 0;

mdot_dist = par.thermal.GOR_md .* Qdot_th ./ par.thermal.h_fg;
mdot_dist(T_b <= stream.T0) = 0;

out.T_b = T_b; out.Qdot_th = Qdot_th; out.W_Carnot = W_Carnot;
out.P_ORC_net = P_ORC_net; out.mdot_dist = mdot_dist;
end
