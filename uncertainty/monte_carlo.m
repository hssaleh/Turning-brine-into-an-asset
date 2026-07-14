function mc = monte_carlo(plant, par)
%MONTE_CARLO Latin-hypercube uncertainty propagation and regression-based
%global sensitivity analysis for the brine-valorization model.
%
% Purpose: Implements Eqs. (10.1)-(10.3):
%   (1) draws N_MC Latin-hypercube samples of the uncertain inputs and
%       maps them through their specified marginal distributions
%       (lognormal / triangular / uniform) by inverse-CDF transform;
%   (2) re-runs the full deterministic pipeline (run_pipeline) under each
%       draw and records the scalar fleet KPIs (fleet_kpis);
%   (3) summarises each output as P5/P50/P95 and coefficient of variation;
%   (4) computes standardized regression coefficients (SRC) as a global
%       sensitivity measure (SRC^2 approximates the first-order Sobol
%       index for near-linear/monotonic response), with the linear-model
%       R^2 reported so the user can judge when SRC under-explains a
%       nonlinear output.
% Inputs : plant - inventory table (a fixed realisation; epistemic
%             parameter uncertainty only is propagated, not inventory
%             sampling error)
%          par   - baseline parameter struct
% Outputs: mc (struct):
%             X [N_MC x K] sampled inputs, var_names {1xK}
%             Y [N_MC x M] sampled outputs, kpi_names {1xM}
%             stats (P5,P50,P95,mean,cov per KPI)
%             SRC [K x M], R2_reg [1 x M]
%             tornado (struct, for a focus KPI): names, src sorted
% Assumptions: inventory is held fixed (the design separates parametric
%          from structural/inventory uncertainty); EF_grid uncertainty is
%          applied as a fleet-wide multiplicative scaling of the per-plant
%          grid factors about their baseline mean.
% Reference: model Sec. 10, Eqs. (10.1)-(10.3).

rng(par.unc.seed, 'twister');
N = par.unc.N_MC;
vnames = fieldnames(par.unc.vars);
K = numel(vnames);

%% Eq. 10.1 Latin-hypercube stratified-uniform design, then inverse-CDF map
U = lhs_uniform(N, K);
X = zeros(N, K);
for j = 1:K
    spec = par.unc.vars.(vnames{j});
    X(:,j) = icdf_var(U(:,j), spec);
end

%% Eq. 10.2 propagate each draw through the deterministic pipeline
EFbase = par.lca.EF_grid;             % baseline fleet grid factor
EFmeanInv = 1./mean(plant.EF_grid);   % to rescale per-plant factors about a new mean

probe = fleet_kpis(plant, run_pipeline(plant, par));
kpi_names = fieldnames(probe);
M = numel(kpi_names);
Y = zeros(N, M);

for n = 1:N
    p = par;
    draw = X(n,:);
    % Map sampled inputs onto parameter struct
    for j = 1:K
        switch vnames{j}
            case 'p_elec'
                p.econ.p_elec = draw(j); p.econ.p_elec_credit = draw(j);
            case 'p_Li2CO3'
                p.m3.p_Li2CO3 = draw(j);
            case 'p_NaCl'
                p.m1.p_NaCl = draw(j);
            case 'c_mem'
                p.econ.c_mem = draw(j);
            case 'eta_Li'
                p.m3.eta_Li = draw(j);
            case 'EF_grid'
                p.lca.EF_grid = draw(j);
            case 'i'
                p.econ.i = draw(j);
        end
    end
    % Per-plant grid factor rescaled to the drawn fleet-mean EF_grid
    pl = plant;
    pl.EF_grid = plant.EF_grid .* (p.lca.EF_grid .* EFmeanInv);
    p.lca.EF_grid = EFbase;   % restore scalar default (per-plant carries the draw)

    R = run_pipeline(pl, p);
    kk = fleet_kpis(pl, R);
    for m = 1:M
        Y(n,m) = kk.(kpi_names{m});
    end
end

%% Eq. 10.3 output statistics
stats = struct();
for m = 1:M
    col = Y(:,m);
    stats.(kpi_names{m}) = struct( ...
        'P5',  prctile(col,5), 'P50', prctile(col,50), 'P95', prctile(col,95), ...
        'mean',mean(col), 'std', std(col), ...
        'cov', std(col)./max(abs(mean(col)),eps) );
end

%% Standardized regression coefficients (global sensitivity)
Xc = (X - mean(X))./std(X);
SRC = zeros(K, M); R2_reg = zeros(1, M);
for m = 1:M
    yc = Y(:,m) - mean(Y(:,m)); sy = std(Y(:,m));
    if sy < eps
        SRC(:,m) = 0; R2_reg(m) = 0; continue;
    end
    yc = yc./sy;
    b = Xc \ yc;                       % standardized OLS slopes
    SRC(:,m) = b;
    yhat = Xc*b;
    R2_reg(m) = 1 - sum((yc-yhat).^2)/sum(yc.^2);
end

%% Tornado data for a focus KPI (global NPV)
focus = 'NPV_global_B';
fm = find(strcmp(kpi_names, focus), 1);
[~, ord] = sort(abs(SRC(:,fm)), 'descend');
tornado.kpi = focus;
tornado.var_names = vnames(ord);
tornado.SRC = SRC(ord, fm);

mc.X = X; mc.var_names = vnames; mc.Y = Y; mc.kpi_names = kpi_names;
mc.stats = stats; mc.SRC = SRC; mc.R2_reg = R2_reg; mc.tornado = tornado;
mc.N = N; mc.K = K; mc.M = M;
end

% ----------------------------------------------------------------------
function U = lhs_uniform(N, K)
%LHS_UNIFORM Latin-hypercube stratified uniform sample on (0,1)^K.
U = zeros(N, K);
for j = 1:K
    edges = ((1:N)' - rand(N,1)) / N;   % one stratified point per row band
    U(:,j) = edges(randperm(N));         % independent random permutation per column
end
end

% ----------------------------------------------------------------------
function x = icdf_var(u, spec)
%ICDF_VAR Inverse-CDF map of a uniform u in (0,1) for the supported
%marginal distributions, parameterised by their P5/P95 quantiles.
P5 = spec.P5; P95 = spec.P95;
switch lower(spec.dist)
    case 'lognormal'
        z5 = -1.6448536; z95 = 1.6448536;
        mu = (log(P95)+log(P5))/2;
        sig = (log(P95)-log(P5))/(z95 - z5);
        x = exp(mu + sig.*norminv_local(u));
    case 'triangular'
        W = (P95 - P5)/0.68377;            % symmetric-triangular width from 5-95% span
        a = P5 - 0.158114.*W; b = a + W; c = (a+b)/2;
        x = tri_icdf(u, a, c, b);
    case 'uniform'
        span = (P95 - P5)/0.9;
        a = P5 - 0.05.*span; b = a + span;
        x = a + (b-a).*u;
    otherwise
        error('icdf_var:dist','Unsupported distribution %s', spec.dist);
end
end

% ----------------------------------------------------------------------
function x = tri_icdf(u, a, c, b)
%TRI_ICDF Inverse CDF of a triangular(a,c,b) distribution.
Fc = (c-a)./(b-a);
x = zeros(size(u));
lo = u < Fc;
x(lo)  = a + sqrt(u(lo).*(b-a).*(c-a));
x(~lo) = b - sqrt((1-u(~lo)).*(b-a).*(b-c));
end

% ----------------------------------------------------------------------
function z = norminv_local(u)
%NORMINV_LOCAL Standard-normal inverse CDF (Acklam rational approximation),
%avoiding a Statistics-Toolbox dependency for the lognormal map.
a = [-3.969683028665376e+01  2.209460984245205e+02 -2.759285104469687e+02 ...
      1.383577518672690e+02 -3.066479806614716e+01  2.506628277459239e+00];
b = [-5.447609879822406e+01  1.615858368580409e+02 -1.556989798598866e+02 ...
      6.680131188771972e+01 -1.328068155288572e+01];
c = [-7.784894002430293e-03 -3.223964580411365e-01 -2.400758277161838e+00 ...
     -2.549732539343734e+00  4.374664141464968e+00  2.938163982698783e+00];
d = [ 7.784695709041462e-03  3.224671290700398e-01  2.445134137142996e+00 ...
      3.754408661907416e+00];
plow = 0.02425; phigh = 1 - plow;
z = zeros(size(u));
% lower tail
lo = u < plow;
q = sqrt(-2*log(u(lo)));
z(lo) = (((((c(1)*q+c(2)).*q+c(3)).*q+c(4)).*q+c(5)).*q+c(6)) ./ ...
        ((((d(1)*q+d(2)).*q+d(3)).*q+d(4)).*q+1);
% upper tail
hi = u > phigh;
q = sqrt(-2*log(1-u(hi)));
z(hi) = -(((((c(1)*q+c(2)).*q+c(3)).*q+c(4)).*q+c(5)).*q+c(6)) ./ ...
         ((((d(1)*q+d(2)).*q+d(3)).*q+d(4)).*q+1);
% central region
mid = ~lo & ~hi;
q = u(mid) - 0.5; rr = q.*q;
z(mid) = (((((a(1)*rr+a(2)).*rr+a(3)).*rr+a(4)).*rr+a(5)).*rr+a(6)).*q ./ ...
         (((((b(1)*rr+b(2)).*rr+b(3)).*rr+b(4)).*rr+b(5)).*rr+1);
end
