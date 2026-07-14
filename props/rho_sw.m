function rho = rho_sw(S, T)
%RHO_SW Seawater/brine density, Sharqawy et al. (2010) correlation.
%
% Purpose: Eq. (2.7) of the brine-valorization model. Valid 0-180 g/kg,
%          0-180 degC.
% Inputs : S - salinity [g/kg], T - temperature [K]. May be arrays
%          (elementwise, broadcasting-compatible sizes).
% Outputs: rho - density [kg/m^3]
% Assumptions: Sharqawy correlation coefficients are exact as published;
%          no extrapolation warning issued outside the valid range.
% Reference: Sharqawy, Lienhard & Zubair (2010) Desalin. Water Treat. 16.

t  = T - 273.15;      % degC
Sp = S/1000;          % kg/kg

a1 = 9.9992293e2;  a2 = 2.0341179e-2; a3 = -6.1624591e-3;
a4 = 2.2614460e-5; a5 = -4.6570659e-8;
b1 = 8.020e2;  b2 = -2.001; b3 = 1.677e-2; b4 = -3.060e-5; b5 = -1.613e-5;

rho_w0 = a1 + a2.*t + a3.*t.^2 + a4.*t.^3 + a5.*t.^4;
rho = rho_w0 + Sp.*(b1 + b2.*t + b3.*t.^2 + b4.*t.^3 + b5.*Sp);
end
