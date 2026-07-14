function inv = generate_inventory(par)
%GENERATE_INVENTORY Synthetic, literature-calibrated global desalination
%plant inventory used as a documented substitute for a licensed DesalData
%extract.
%
% Project   : Turning Brine into an Asset
% Purpose   : Build a georeferenced table of N_plant desalination plants
%             (technology, capacity, feed salinity, recovery ratio,
%             region) with regional economic/environmental context
%             vectors (electricity price, grid carbon factor, baseline
%             LCOW), reproducing the heterogeneity (Gulf high-salinity
%             SWRO/MED clusters, Mediterranean/Pacific SWRO, etc.) that a
%             real global fleet exhibits, so that every downstream
%             aggregation, viability-atlas and uncertainty figure is
%             scientifically meaningful even without a proprietary
%             inventory.
% Author    : <AUTHOR>
% Date      : <DATE>
% Version   : 1.0
%
% Inputs  : par (struct) - from config_parameters(), uses par.inv.*
% Outputs : inv (table) - one row per plant with fields:
%             id, region, tech, C_des [m3/day], S_f [g/kg], r [-],
%             lat, lon, p_elec [USD/kWh], EF_grid [kgCO2e/kWh],
%             LCOW_base [USD/m3]
% Assumptions:
%   - Regional statistics (salinity range, technology mix, price/carbon
%     factors) are representative literature-centred ranges, not a
%     measured census; this is stated explicitly wherever the inventory
%     is used downstream (Brine_Valorization_Mathematical_Model.md Sec. 9
%     calls for a "georeferenced inventory ... DesalData-derived"; no
%     such license is available in this workspace).
% References: Jones et al. (2019) Sci. Total Environ.; model Sec. 9.

rng(par.inv.seed, 'twister');
N = par.inv.N_plant;

% Region definitions: name, lon/lat box, tech-mix weights [SWRO MED MSF],
% salinity range [g/kg], electricity price range [USD/kWh],
% grid carbon factor range [kgCO2e/kWh], baseline LCOW range [USD/m3],
% representative share of global capacity (sampling weight).
regions = {
    'ArabianGulf',   [46 56], [22 30], [0.45 0.30 0.25], [42 47], [0.06 0.10], [0.55 0.75], [0.9 1.3], 0.30
    'RedSea',        [34 43], [13 25], [0.55 0.30 0.15], [38 42], [0.05 0.09], [0.50 0.70], [0.85 1.2], 0.08
    'Mediterranean', [0  35], [30 42], [0.85 0.10 0.05], [36 39], [0.07 0.16], [0.25 0.55], [0.7 1.0], 0.14
    'NorthAfrica',   [-10 12], [22 36], [0.75 0.15 0.10], [37 40], [0.04 0.09], [0.55 0.80], [0.65 0.95], 0.07
    'Caribbean',     [-85 -60],[10 22], [0.80 0.10 0.10], [35 37], [0.10 0.22], [0.45 0.70], [0.9 1.4], 0.04
    'USWestCoast',   [-124 -114],[28 49],[0.95 0.03 0.02], [33 35], [0.06 0.14], [0.15 0.35], [0.8 1.2], 0.06
    'Australia',     [113 154],[-38 -12],[0.95 0.03 0.02], [35 36], [0.10 0.20], [0.40 0.70], [0.9 1.3], 0.05
    'EastAsia',      [100 145],[5 40],  [0.80 0.12 0.08], [31 35], [0.06 0.13], [0.45 0.75], [0.7 1.1], 0.10
    'SouthAsia',     [60 90],  [5 30],  [0.85 0.08 0.07], [34 36], [0.05 0.10], [0.55 0.85], [0.6 0.9], 0.09
    'SouthernAfrica',[10 35],  [-35 -10],[0.90 0.06 0.04],[34 36], [0.07 0.13], [0.50 0.80], [0.7 1.0], 0.03
    'SouthAmerica',  [-80 -60],[-40 5], [0.92 0.04 0.04], [34 36], [0.06 0.13], [0.20 0.50], [0.7 1.0], 0.04
};
nReg = size(regions,1);
w = cellfun(@(x) x, regions(:,9));
w = w / sum(w);

regID  = randsample(nReg, N, true, w);
techNames = {'SWRO','MED','MSF'};

id        = (1:N)';
region    = cell(N,1);
tech      = cell(N,1);
S_f       = zeros(N,1);
lon       = zeros(N,1);
lat       = zeros(N,1);
p_elec    = zeros(N,1);
EF_grid   = zeros(N,1);
LCOW_base = zeros(N,1);
C_des     = zeros(N,1);
r         = zeros(N,1);

for k = 1:N
    rg = regions(regID(k), :);
    region{k} = rg{1};
    lonRange = rg{2}; latRange = rg{3};
    techW = rg{4}; Srange = rg{5};
    pelecR = rg{6}; efgridR = rg{7}; lcowR = rg{8};

    lon(k) = lonRange(1) + diff(lonRange)*rand();
    lat(k) = latRange(1) + diff(latRange)*rand();

    tech{k} = techNames{randsample(3,1,true,techW)};
    S_f(k)  = Srange(1) + diff(Srange)*rand();
    p_elec(k)    = pelecR(1) + diff(pelecR)*rand();
    EF_grid(k)   = efgridR(1) + diff(efgridR)*rand();
    LCOW_base(k) = lcowR(1) + diff(lcowR)*rand();

    % Plant capacity: log-uniform between 2,000 and 600,000 m3/day,
    % spanning small island SWRO units to mega-plants (e.g. Ras Al-Khair).
    C_des(k) = 10^(log10(2000) + (log10(600000)-log10(2000))*rand());

    switch tech{k}
        case 'SWRO'
            r(k) = 0.35 + 0.15*rand();   % 0.35-0.50
        otherwise % MED, MSF (thermal)
            r(k) = 0.25 + 0.10*rand();   % 0.25-0.35
    end
end

inv = table(id, string(region), string(tech), C_des, S_f, r, lat, lon, ...
            p_elec, EF_grid, LCOW_base, ...
            'VariableNames', {'id','region','tech','C_des','S_f','r','lat','lon', ...
                               'p_elec','EF_grid','LCOW_base'});
end
