function [pi_p, a_w] = osmotic_pressure(m, phi, T, par)
%OSMOTIC_PRESSURE Rigorous osmotic pressure and water activity of brine.
%
% Purpose: Eq. (2.13) of the model: pi = phi*nu*m_s*R*T*rho_w, with
%          a_w = exp(-phi*nu*m_s*M_w).
% Inputs : m   - NaCl-equivalent salt molality [mol/kg water]
%          phi - osmotic coefficient (-) from activity_coeffs.m
%          T   - temperature [K]
%          par - parameter struct
% Outputs: pi_p - osmotic pressure [Pa]
%          a_w  - water activity (-)
% Assumptions: nu = 2 (NaCl-equivalent ions per formula unit).
% Reference: model Sec. 2.5, Eq. (2.13).

R    = par.const.R;
rho_w = par.const.rho_w;
M_w   = par.const.M_w;
nu    = par.stream.nu_NaCl;

a_w  = exp(-phi.*nu.*m.*M_w);
pi_p = phi.*nu.*m.*R.*T.*rho_w;
end
