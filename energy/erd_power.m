function out = erd_power(stream, par)
%ERD_POWER Hydraulic energy recovery (pressure exchanger) pathway.
%
% Purpose: Implements Eqs. (3.19)-(3.20): recovery of the RO brine's
%          residual hydraulic energy via an isobaric pressure exchanger.
% Inputs : stream - struct from M00_streams.m (Q_b, P_b)
%          par    - parameter struct (par.erd.eta_ERD, par.red.eta_DC)
% Outputs: out.P_gross, out.P_net [W]
% Assumptions: P_b = 0 for MED/MSF brine (no residual RO pressure), set
%          in M00_streams.m; the isobaric-exchanger option (eta_ERD) is
%          used rather than the Pelton-turbine alternative.
% Reference: model Sec. 3.3, Eqs. (3.19)-(3.20).

P_gross = stream.Q_b .* stream.P_b .* par.erd.eta_ERD;
P_net   = par.red.eta_DC .* P_gross;

out.P_gross = P_gross; out.P_net = P_net;
end
