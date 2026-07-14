function mu = mu_sw(S, T)
%MU_SW Seawater/brine dynamic viscosity, Sharqawy et al. (2010) correlation.
%
% Purpose: Eq. (2.8a) of the brine-valorization model.
% Inputs : S - salinity [g/kg], T - temperature [K]
% Outputs: mu - dynamic viscosity [Pa s]
% Assumptions: valid 0-150 g/kg, 0-180 degC (Sharqawy 2010).
% Reference: Sharqawy, Lienhard & Zubair (2010) Desalin. Water Treat. 16.

t  = T - 273.15;
Sp = S/1000;

mu_w = 4.2844e-5 + 1./(0.157.*(t+64.993).^2 - 91.296);
A_mu = 1.541 + 1.998e-2.*t - 9.52e-5.*t.^2;
B_mu = 7.974 - 7.561e-2.*t + 4.724e-4.*t.^2;
mu   = mu_w .* (1 + A_mu.*Sp + B_mu.*Sp.^2);
end
