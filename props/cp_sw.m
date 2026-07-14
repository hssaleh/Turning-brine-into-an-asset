function cp = cp_sw(S, T)
%CP_SW Seawater/brine specific heat, Sharqawy et al. (2010) correlation.
%
% Purpose: Eq. (2.8b) of the brine-valorization model.
% Inputs : S - salinity [g/kg], T - temperature [K]
% Outputs: cp - specific heat [J kg^-1 K^-1]
% Assumptions: valid 0-180 g/kg, 0-180 degC.
% Reference: Sharqawy, Lienhard & Zubair (2010) Desalin. Water Treat. 16.

A_cp = 5328 - 9.76e1.*S + 4.04e-1.*S.^2;
B_cp = -6.913 + 7.351e-1.*S - 3.15e-3.*S.^2;
C_cp = 9.6e-3 - 1.927e-3.*S + 8.23e-6.*S.^2;
D_cp = 2.5e-6 + 1.666e-6.*S - 7.125e-9.*S.^2;

cp = A_cp + B_cp.*T + C_cp.*T.^2 + D_cp.*T.^3;
end
