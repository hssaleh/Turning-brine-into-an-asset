function s = s_sw(S, T, T0)
%S_SW Specific physical entropy of seawater/brine relative to the dead
%state, by closed-form integration of cp(S,T)/T.
%
% Purpose: Provides s(S,T) for the exergy module (Eq. 5.2), consistent
%          with h_sw.m (same cp polynomial, same dead-state reference).
% Inputs : S [g/kg], T [K], T0 [K]
% Outputs: s [J kg^-1 K^-1]
% Assumptions: see h_sw.m.
% Reference: Sharqawy, Lienhard & Zubair (2010); model Sec. 5.

A_cp = 5328 - 9.76e1.*S + 4.04e-1.*S.^2;
B_cp = -6.913 + 7.351e-1.*S - 3.15e-3.*S.^2;
C_cp = 9.6e-3 - 1.927e-3.*S + 8.23e-6.*S.^2;
D_cp = 2.5e-6 + 1.666e-6.*S - 7.125e-9.*S.^2;

s = A_cp.*log(T./T0) + B_cp.*(T-T0) + (C_cp/2).*(T.^2-T0.^2) + (D_cp/3).*(T.^3-T0.^3);
end
