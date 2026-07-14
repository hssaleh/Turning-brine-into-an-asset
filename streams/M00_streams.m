function stream = M00_streams(plant, par)
%M00_STREAMS Stream characterization and thermophysical property module.
%
% Project   : Turning Brine into an Asset
% Purpose   : Converts each plant record (capacity, technology, feed
%             salinity, recovery) into the complete set of stream states
%             required by the recovery models, per Sec. 2 (Eqs. 2.1-2.14)
%             of the mathematical model.
% Author    : <AUTHOR> | Date: <DATE> | Version: 1.0
%
% Inputs  : plant - table with columns C_des [m3/day], S_f [g/kg],
%             r [-], tech {string}; one row per plant, vectorised over
%             rows (N plants).
%           par   - parameter struct from config_parameters().
% Outputs : stream (struct), every field is an N x 1 vector:
%             Q_p, Q_f, Q_b [m3/s]; S_b [g/kg]; rho_f, rho_b [kg/m3];
%             mu_b [Pa s]; cp_b [J/kg/K]; I_s_b, I_s_f [mol/kg];
%             m_s_b, m_s_f [mol/kg] (NaCl-equivalent molality);
%             c_NaCl_b, c_NaCl_f [mol/m3]; pi_b, pi_f [Pa]; a_w_b, a_w_f;
%             gamma_b, gamma_f; phi_b, phi_f; T_f, T0, P_b [K, K, Pa];
%             m_i_b (N x 12 ion molality matrix).
% Assumptions: S_p (permeate salinity) and R_salt as in par.stream;
%             T_f = T0 = par.stream defaults unless plant supplies its
%             own temperature; brine gauge pressure P_b = 5.5 MPa for
%             SWRO, 0 for thermal technologies (MED/MSF), per Sec. 3.3.
% Reference : model Sec. 2, Eqs. (2.1)-(2.14).

N = height(plant);
f_av   = par.stream.f_av;
r      = plant.r(:);
S_f    = plant.S_f(:);
C_des  = plant.C_des(:);
S_p    = par.stream.S_p;
R_salt = par.stream.R_salt;
T_f    = par.stream.T_f * ones(N,1);
T0     = par.stream.T0 * ones(N,1);

%% Eq. 2.1-2.3 mass/volume balance
Q_p = C_des .* f_av / 86400;          % m^3/s
Q_f = Q_p ./ r;                        % m^3/s
Q_b = Q_p .* (1-r) ./ r;               % m^3/s

%% Eq. 2.6 brine salinity (salt-rejection refinement)
S_b = S_f .* (1 - r.*(1-R_salt)) ./ (1-r);   % g/kg

%% Physical-bounds checks (Sec. 2.5 closure requirements)
assert(all(r>0 & r<1), 'M00_streams: recovery ratio out of (0,1) bounds');
assert(all(S_b > S_f), 'M00_streams: brine salinity not greater than feed');

%% Eq. 2.7 density (feed and brine)
rho_f = rho_sw(S_f, T_f);
rho_b = rho_sw(S_b, T_f);
assert(all(rho_b > rho_f), 'M00_streams: brine density not greater than feed');

%% Eq. 2.8 viscosity / specific heat (brine, used by downstream modules)
mu_b = mu_sw(S_b, T_f);
cp_b = cp_sw(S_b, T_f);

%% Brine exit pressure by technology (Sec. 3.3)
P_b = zeros(N,1);
if any(strcmpi(plant.tech,'SWRO'))
    P_b(strcmpi(plant.tech,'SWRO')) = par.stream.P_b;
end
% MED/MSF: negligible residual hydraulic pressure
P_b(~strcmpi(plant.tech,'SWRO')) = 0;

%% Eq. 2.9 ionic speciation (all 12 ions), brine and feed
nIon = numel(par.ions.name);
m_i_b = zeros(N, nIon);
m_i_f = zeros(N, nIon);
for k = 1:nIon
    w35 = par.ions.w35(k); Mi = par.ions.M(k);
    m_i_b(:,k) = (w35/1000).*(S_b/35) ./ (Mi.*(1000-S_b)/1000);
    m_i_f(:,k) = (w35/1000).*(S_f/35) ./ (Mi.*(1000-S_f)/1000);
end

%% Eq. 2.11 ionic strength
z2 = par.ions.z.^2;
I_s_b = 0.5 * (m_i_b * z2(:));
I_s_f = 0.5 * (m_i_f * z2(:));

%% NaCl-equivalent lumped salt molality and molar concentration (Eq. 2.9-2.10 reduction)
M_NaCl = par.const.M_NaCl;
m_s_b = (S_b./(1000-S_b)) / M_NaCl;     % mol/kg water
m_s_f = (S_f./(1000-S_f)) / M_NaCl;
c_NaCl_b = m_s_b .* rho_b ./ (1+S_b/1000);  % mol/m^3
c_NaCl_f = m_s_f .* rho_f ./ (1+S_f/1000);  % mol/m^3

%% Eq. 2.12 / 4.1-4.3 activity & osmotic coefficients (Davies-Pitzer blend)
[gamma_b, phi_b] = activity_coeffs(m_s_b, I_s_b, par);
[gamma_f, phi_f] = activity_coeffs(m_s_f, I_s_f, par);

%% Eq. 2.13 osmotic pressure and water activity
[pi_b, a_w_b] = osmotic_pressure(m_s_b, phi_b, T_f, par);
[pi_f, a_w_f] = osmotic_pressure(m_s_f, phi_f, T_f, par);
assert(all(pi_b > pi_f & pi_f > 0), 'M00_streams: osmotic pressure ordering/positivity violated');

%% Pack outputs
stream.Q_p = Q_p;  stream.Q_f = Q_f;  stream.Q_b = Q_b;
stream.S_f = S_f;  stream.S_b = S_b;
stream.rho_f = rho_f; stream.rho_b = rho_b;
stream.mu_b = mu_b;   stream.cp_b = cp_b;
stream.I_s_b = I_s_b; stream.I_s_f = I_s_f;
stream.m_s_b = m_s_b; stream.m_s_f = m_s_f;
stream.c_NaCl_b = c_NaCl_b; stream.c_NaCl_f = c_NaCl_f;
stream.pi_b = pi_b;   stream.pi_f = pi_f;
stream.a_w_b = a_w_b; stream.a_w_f = a_w_f;
stream.gamma_b = gamma_b; stream.gamma_f = gamma_f;
stream.phi_b = phi_b; stream.phi_f = phi_f;
stream.T_f = T_f;     stream.T0 = T0;
stream.P_b = P_b;
stream.m_i_b = m_i_b; stream.m_i_f = m_i_f;
stream.tech = plant.tech;
stream.region = plant.region;
end
