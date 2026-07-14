function h = h_sw(S, T, T0)
%H_SW Specific physical enthalpy of seawater/brine relative to the dead
%state, by closed-form integration of the Sharqawy cp(S,T) polynomial.
%
% Purpose: Provides h(S,T) for the exergy module (Eq. 5.2), with the
%          reference enthalpy set to zero at the dead state (S,T0), per
%          Sec. 5.2 of the model ("baseline reference enthalpy at dead
%          state set to zero").
% Inputs : S [g/kg], T [K], T0 [K] dead-state temperature
% Outputs: h [J/kg]
% Assumptions: physical (thermal) enthalpy only; chemical part is
%          accounted separately via ex_ch (Eq. 5.3). cp is evaluated at
%          fixed S (no salinity-mixing enthalpy in this physical term).
% Reference: Sharqawy, Lienhard & Zubair (2010); model Sec. 5.

A_cp = 5328 - 9.76e1.*S + 4.04e-1.*S.^2;
B_cp = -6.913 + 7.351e-1.*S - 3.15e-3.*S.^2;
C_cp = 9.6e-3 - 1.927e-3.*S + 8.23e-6.*S.^2;
D_cp = 2.5e-6 + 1.666e-6.*S - 7.125e-9.*S.^2;

h = A_cp.*(T-T0) + (B_cp/2).*(T.^2-T0.^2) + (C_cp/3).*(T.^3-T0.^3) + (D_cp/4).*(T.^4-T0.^4);
end
