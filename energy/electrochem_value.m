function out = electrochem_value(stream, par)
%ELECTROCHEM_VALUE Brine (chlor-alkali) electrolysis conversion option.
%
% Purpose: Implements Eqs. (3.25)-(3.30): reversible cell potential,
%          operating voltage, hydrogen/chlorine/caustic production
%          rates, specific energy consumption and net electrochemical
%          value rate. Scored as a conversion/credit option (not net
%          primary energy harvested from the salinity gradient).
% Inputs : stream - struct from M00_streams.m (Q_b)
%          par    - parameter struct (par.ec.*)
% Outputs: out.mdot_H2, out.mdot_Cl2, out.mdot_NaOH [kg/s]; out.SEC_H2
%          [kWh/kg]; out.P_demand [W]; out.Vdot_EC [USD/s]
% Assumptions: only a fraction (par.ec.frac_brine_to_EC) of brine flow is
%          routed to electrolysis (pilot-scale option, not full-stream);
%          E_rev approximated by the standard potential E0 (chloride
%          activity near unity in concentrated brine; activation/ohmic
%          overpotentials carry the composition dependence explicitly).
% Reference: model Sec. 3.5, Eqs. (3.25)-(3.30).

F = par.const.F;
Q_EC = par.ec.frac_brine_to_EC .* stream.Q_b;
A_cell = par.ec.A_cell_per_Qb .* Q_EC;
I = par.ec.i_cd .* A_cell;

eta_ohm = par.ec.i_cd .* par.ec.R_area;
V_cell = par.ec.E0 + par.ec.eta_act + eta_ohm;   % Eq. 3.25-3.26 (E_rev ~ E0)

mdot_H2   = (I .* par.ec.eta_F .* par.ec.M_H2)   ./ (par.ec.z_H2 .* F);
SEC_H2    = (par.ec.z_H2 .* F .* V_cell) ./ (par.ec.eta_F .* par.ec.M_H2) / 3.6e6;
mdot_Cl2  = (I .* par.ec.eta_F .* par.ec.M_Cl2)  ./ (2.*F);
mdot_NaOH = (I .* par.ec.eta_F .* par.ec.M_NaOH) ./ F;

P_demand = I .* V_cell;   % W, electrical demand of the EC cell
Vdot_EC = par.ec.p_H2.*mdot_H2 + par.ec.p_Cl2.*mdot_Cl2 + par.ec.p_NaOH.*mdot_NaOH ...
          - par.econ.p_elec.*(P_demand)/3.6e6;   % USD/s (p_elec is USD/kWh; P_demand W -> kWh/s = P/3.6e6)

out.I = I; out.V_cell = V_cell; out.mdot_H2 = mdot_H2; out.SEC_H2 = SEC_H2;
out.mdot_Cl2 = mdot_Cl2; out.mdot_NaOH = mdot_NaOH; out.P_demand = P_demand;
out.Vdot_EC = Vdot_EC;
end
