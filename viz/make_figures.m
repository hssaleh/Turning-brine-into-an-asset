function FIGS = make_figures(inv, R, mc, lb, par, figdir)
%MAKE_FIGURES Publication-quality figure suite for the brine-valorization
%study, grouped with tiledlayout and exported as .fig.
%
% Purpose: Generates the full set of scientifically-meaningful, high-
%          impact figures spanning stream/mass, energy, exergy/entropy,
%          first- and second-law efficiency, dimensionless analysis,
%          techno-economics, environment/LCA, uncertainty and the global
%          viability atlas. Every figure carries a publication caption, a
%          detailed explanation and a results interpretation (string
%          variables) and an attached numeric-data block set for one-
%          sheet-per-figure Excel export.
% Inputs : inv - inventory table; R - run_pipeline output; mc - monte_carlo
%          output; lb - literature_benchmarks output; par - parameter
%          struct; figdir - folder for .fig files.
% Outputs: FIGS - struct array with fields name,title,caption,explanation,
%          interpretation,blocks(cell of struct('label',str,'T',table)),h.
% 

if nargin<6 || isempty(figdir), figdir = fullfile(pwd,'results','figures'); end
if ~exist(figdir,'dir'), mkdir(figdir); end
set(groot,'defaultAxesFontName','Arial','defaultAxesFontSize',11);
FIGS = struct('name',{},'title',{},'caption',{},'explanation',{},'interpretation',{},'blocks',{},'h',{});

N   = height(inv);
S_b = R.stream.S_b; S_f = R.stream.S_f; r = inv.r;
regions = unique(inv.region,'stable');
col = lines(7); ok = @(ax) set(ax,'Box','on','LineWidth',0.9,'XGrid','on','YGrid','on');

% nested helper to register a figure
    function add_fig(h, name, ttl, cap, expl, interp, blocks)
        sgtitle(h, ttl, 'FontWeight','bold','FontSize',13);
        try, savefig(h, fullfile(figdir,[name '.fig'])); catch, end
        FIGS(end+1) = struct('name',name,'title',ttl,'caption',cap, ...
            'explanation',expl,'interpretation',interp,'blocks',{blocks},'h',h);
    end
    function b = blk(label, T), b = struct('label',label,'T',T); end


%% create article figures (reorganized)
%% ===================== FIG 01 — Stream characterization =================
h = figure('Name','F02_Manuscript','Color','w','Position',[60 60 1100 760],'Visible','on');
tiledlayout(h,2,2,'Padding','compact','TileSpacing','compact');

ax=nexttile; 
sz = 5 + 45*rescale(r);
scatter(ax,S_f,S_b,sz,r,'filled'); ok(ax); colormap(ax,parula); cb=colorbar(ax);
cb.Label.String='Recovery ratio r (-)'; xlabel(ax,'Feed salinity S_f (g kg^{-1})');
ylabel(ax,'Brine salinity S_b (g kg^{-1})'); title(ax,'(a) Brine Concentration');

ax=nexttile; 
sz = 5 + 45*rescale(R.stream.pi_f/1e5);
scatter(ax,S_b,R.stream.pi_b/1e5,sz,R.stream.pi_f/1e5,'filled'); ok(ax);
cb=colorbar(ax); cb.Label.String='\pi_f (bar)'; xlabel(ax,'Brine salinity S_b (g kg^{-1})');
ylabel(ax,'Brine osmotic pressure \pi_b (bar)'); title(ax,'(b) Osmotic pressure');

Pmat = [R.energy.red.P_net R.energy.pro.P_net R.energy.erd.P_net R.energy.thermal.P_ORC_net]/1e3; % kW
Pmat(Pmat<=0)=NaN;
ax=nexttile; boxchart(ax,Pmat); 

% ylim([1e-0 1e5])


 box(ax,'on');
 set(ax,'YScale','log')
grid on
xticklabels(ax,{'RED','PRO','ERD','ORC'}); ylabel(ax,'Net power (kW)');
title(ax,'(c) Energy-recovery power'); set(ax,'YScale','linear');
ax.YScale = 'log';

% Technically available = pathway delivers a positive product/output where
% admissible; economically selected = included in the optimal portfolio.
np = numel(R.econ.pathway_names);
avail = zeros(1,np);
avail(1)=mean(R.energy.red.P_net>0);          avail(2)=mean(R.energy.pro.P_net>0);
avail(3)=mean(R.energy.erd.P_net>0);          avail(4)=mean(R.energy.thermal.P_ORC_net>0);
avail(5)=mean(R.energy.ec.mdot_H2>0);         avail(6)=mean(sum(R.mineral.m1.mdot_solid(:,1:2),2)>0);
avail(7)=mean(R.mineral.m2.mdot_NaCl>0);      avail(8)=mean(R.mineral.m3.mdot_Li2CO3>0);
avail(9)=mean(R.mineral.m4.mdot_Br>0);        avail(10)=mean(R.mineral.m5.Q_water_add>0);
sf = mean(R.econ.sel,1);
ax=nexttile; bA=bar(ax,100*[avail(:) sf(:)],'grouped'); ok(ax);
bA(1).FaceColor=[0.7 0.78 0.9]; bA(2).FaceColor=col(4,:);
set(ax,'XTick',1:np,'XTickLabel',R.econ.pathway_names,'XTickLabelRotation',35);
ylabel(ax,'Fraction of fleet (%)'); legend(ax,{'Technically available','Economically selected'},'Location','north','Orientation','horizontal');
title(ax,'(d) Availability vs selection'); ylim(ax,[0 115]);


%%
%% ===================== FIG 03 — Minerals & market absorption ============
mp_tech = R.mineral.mprod_yr;          % technical (chemically recoverable) mass
mp_sel  = R.econ.mprod_selected;       % economically-selected mass
pn = R.mineral.product_names;
techG = R.agg.global.mineral_Mt_yr_tech; consG = R.agg.global.mineral_Mt_yr_constr;
h = figure('Name','F03_Manuscript','Color','w','Position',[80 80 1400 400],'Visible','on');
tiledlayout(h,1,3,'Padding','compact','TileSpacing','compact');

ax=nexttile; Yt=[sum(mp_tech,1); sum(mp_sel,1)]'/1e3; b=bar(ax,Yt,'grouped'); ok(ax);
b(1).FaceColor=col(1,:); b(2).FaceColor=col(2,:);
set(ax,'XTick',1:7,'XTickLabel',pn,'XTickLabelRotation',30); 
set(ax,'YScale','log'); ylim([1e-1 1e6])
ylabel(ax,'Fleet production (kt yr^{-1})'); 
legend(ax,{'Technical','Selected'},'Location','northeast');
title(ax,'(a) Recovered product mass');

% ax=nexttile; idxNa=3; scatter(ax,R.stream.m_i_b(:,strcmpi(par.ions.name,'Na')).*R.stream.Q_b,mp_tech(:,idxNa)/1e3,22,col(3,:),'filled'); ok(ax);
% xlabel(ax,'Brine Na^{+} molar flow (mol s^{-1})'); ylabel(ax,'NaCl production (kt yr^{-1})'); title(ax,'(b) Halite (NaCl) yield vs Na^{+}');

ax=nexttile; Y=[techG(:) consG(:)]; b=bar(ax,Y,'grouped'); ok(ax); set(ax,'YScale','log');
set(ax,'XTick',1:7,'XTickLabel',pn,'XTickLabelRotation',30); ylabel(ax,'Global potential (Mt yr^{-1})');
legend(ax,{'Technical','Market-absorbable'},'Location','northeast'); title(ax,'(b) Technical vs market-absorbable');
b(1).FaceColor=col(1,:); b(2).FaceColor=col(2,:);

ax=nexttile; fa=100*R.agg.f_absorb; bar(ax,fa,'FaceColor',col(5,:)); ok(ax);
set(ax,'XTick',1:7,'XTickLabel',pn,'XTickLabelRotation',30); ylabel(ax,'Absorbable fraction (%)');
title(ax,'(c) Market-absorption ceiling'); 
box(ax,'on');
 set(ax,'YScale','log')
 ylim([1e-1 200])

for kk=1:7, text(ax,kk,fa(kk)*1.3,sprintf('%.1f%%',fa(kk)),'HorizontalAlignment','center','FontSize',8); end




%% ===================== FIG 04 — Techno-economics =======================
% Standalone fleet CAPEX of each pathway (cost if deployed at every plant),
% so all pathways - including the technically-capable but unselected RED/PRO -
% appear; the selected-portfolio CAPEX is overlaid for comparison.
capStand = sum(R.econ.CAPEX_k, 1) / 1e9;               % standalone, B USD
capSel   = sum(R.econ.sel .* R.econ.CAPEX_k, 1) / 1e9; % selected, B USD
capLbl   = R.econ.pathway_names;
h = figure('Name','F04 - Manuscript','Color','w','Position',[70 60  1100 760],'Visible','on');
tiledlayout(h,2,2,'Padding','compact','TileSpacing','compact');

ax=nexttile; histogram(ax,R.econ.LCOW_net,30,'FaceColor',col(1,:)); ok(ax); xline(ax,0,'r--','LineWidth',1.4);
text(ax,-3,25,{'mean = -0.44';'61% < 0'})
xline(0,'--r')
xline(-0.440,'-b')

ylim([0 30]);xlim([-3.5 2])
xlabel(ax,'Net-LCOW (USD m^{-3})'); ylabel(ax,'Plant count'); title(ax,'(a) Net-LCOW distribution');


ax=nexttile; 

sz = 5 + 45*rescale(R.agg.VI);

scatter(ax,R.econ.CAPEX_total/1e6,R.econ.NPV/1e6,sz,R.agg.VI,'filled'); ok(ax);
colormap(ax,turbo); cb=colorbar(ax); cb.Label.String='Viability index (-)';
xlim([0.5 1e3*2])
xlabel(ax,'CAPEX (10^{6} USD)'); ylabel(ax,'NPV (10^{6} USD)'); title(ax,'(b) NPV vs CAPEX'); set(ax,'XScale','log','YScale','log');

ax=nexttile; %yyaxis(ax,'left'); 
histogram(ax,100*R.econ.IRR,30); ylabel(ax,'Plant count');
% yyaxis(ax,'right'); pb=R.econ.payback; pb(~isfinite(pb))=NaN; [f,x]=ecdf(pb); stairs(ax,x,f,'LineWidth',1.6); ylabel(ax,'Payback CDF (-)');
ok(ax); xlabel(ax,'Internal rate of return (%)  '); title(ax,'(c) Profitability');
xline(31,'-b')

ax=nexttile; histogram(ax,R.econ.BCR,30,'FaceColor',col(3,:)); ok(ax); xline(ax,1,'r--','LineWidth',1.4);
xlabel(ax,'Benefit-cost ratio (-)'); ylabel(ax,'Plant count'); title(ax,'(d) Benefit-cost ratio distribution');
xline(1,'--r')
xline(1.5,'-b')


%% ===================== FIG 11 — Global resource atlas ==================
h = figure('Name','F05_Manuscript','Color','w','Position',[50 50 1100 760],'Visible','on');
tiledlayout(h,2,2,'Padding','compact','TileSpacing','compact');

ax=nexttile([1 2]); hold(ax,'on');
draw_coastlines(ax);                              % world coastlines backdrop
sz = 5 + 45*rescale(R.agg.VI);
scatter(ax,inv.lon,inv.lat,sz,R.agg.VI,'filled','MarkerEdgeColor',[.25 .25 .25]);
hold(ax,'off'); ok(ax); colormap(ax,turbo); cb=colorbar(ax); cb.Label.String='Viability index (-)';
xlabel(ax,'Longitude (\circ)'); ylabel(ax,'Latitude (\circ)'); title(ax,'(a) Global brine-valorization viability atlas');
xlim(ax,[-160 160]); ylim(ax,[-50 60]);


ax=nexttile; 
sz = 5 + 45*rescale(inv.p_elec);
scatter(ax, S_b, R.agg.VI, sz, inv.p_elec, 'filled');
% scatter(ax,S_b,R.agg.VI,24,inv.p_elec,'filled'); 
ok(ax); colormap(ax,'parula');
cb=colorbar(ax); cb.Label.String='p_{elec} (USD kWh^{-1})'; 
xlabel(ax,'Brine salinity S_b (g kg^{-1})');
ylabel(ax,'Viability index (-)'); title(ax,'(b) Salinity driver'); ylim([0.28 0.92])


regNPV=zeros(numel(regions),1); regE=zeros(numel(regions),1);
for g=1:numel(regions), m=inv.region==regions(g); regNPV(g)=R.agg.GRF*sum(R.econ.NPV(m))/1e9; regE(g)=R.agg.GRF*sum(R.econ.E_rec_sel_yr(m))/1e9; end
ax=nexttile; [rs,oi]=sort(regNPV,'ascend'); barh(ax,rs,'FaceColor',col(1,:)); ok(ax);
set(ax,'YTick',1:numel(regions),'YTickLabel',regions(oi),'YTickLabelRotation',35);
xlabel(ax,'Scaled NPV (10^{9} USD)'); title(ax,'(c) Regional NPV potential');

% ax=nexttile; bar(ax,regE(oi),'FaceColor',col(2,:)); ok(ax);
% set(ax,'XTick',1:numel(regions),'XTickLabel',regions(oi),'XTickLabelRotation',35); ylabel(ax,'Recovered energy (TWh yr^{-1})'); title(ax,'(c) Regional energy potential');


%% ===================== FIG 10 — Uncertainty (Monte Carlo) ==============
h = figure('Name','F06_Manuscript','Color','w','Position',[60 60 1100 760],'Visible','on');
tiledlayout(h,2,2,'Padding','compact','TileSpacing','compact');
gi=@(nm) find(strcmp(mc.kpi_names,nm),1);

ax=nexttile; y=mc.Y(:,gi('NPV_global_B')); histogram(ax,y,30,'FaceColor',col(1,:)); ok(ax);
s=mc.stats.NPV_global_B; xline(ax,s.P5,'k--'); xline(ax,s.P50,'r-','LineWidth',1.4); xline(ax,s.P95,'k--');
xlabel(ax,'Global NPV (10^{9} USD)'); ylabel(ax,'Draw count'); title(ax,'(a) Global NPV distribution');
xlim([0 2000])


ax=nexttile; y=mc.Y(:,gi('LCOW_net_mean')); [f,x]=ecdf(y); 
stairs(ax,x,f,'LineWidth',1.8,'Color',col(2,:)); ok(ax);
xlabel(ax,'Fleet mean net-LCOW (USD m^{-3})'); 
ylabel(ax,'Cumulative probability'); xline(ax,0,'r--'); title(ax,'(b) Net-LCOW CDF');
xlim([-4.75 1.25])
ylim([0 1.1])

ax=nexttile; 
sz = 5 + 45*rescale(mc.Y(:,gi('VI_mean')));
% scatter(ax,mc.Y(:,gi('')),mc.Y(:,gi('CO2_mit_Mt')),26,mc.Y(:,gi('VI_mean')),'filled');
scatter(ax,mc.Y(:,gi('NPV_global_B')),mc.Y(:,gi('CO2_mit_Mt')),sz,mc.Y(:,gi('VI_mean')),'filled');
ok(ax); colormap(ax,turbo); cb=colorbar(ax); cb.Label.String='VI (-)';
xlabel(ax,'Global NPV (10^{9} USD)'); ylabel(ax,'CO_2 mitigation (Mt yr^{-1})'); title(ax,'(c) Joint NPV-CO_2 ensemble');


ax=nexttile; bar(ax,regE(oi),'FaceColor',col(2,:)); ok(ax);
set(ax,'XTick',1:numel(regions),'XTickLabel',regions(oi),'XTickLabelRotation',35);
ylabel(ax,'Recovered energy (TWh yr^{-1})'); title(ax,'(d) Regional energy potential');


%%
%%
%%
%% ===================== FIG 01 — Stream characterization =================
h = figure('Name','F01_Streams','Color','w','Position',[60 60 1100 400],'Visible','on');
tiledlayout(h,1,2,'Padding','compact','TileSpacing','compact');


ax=nexttile; yyaxis(ax,'left'); scatter(ax,S_b,R.stream.rho_b,18,'filled'); ylabel(ax,'Brine density \rho_b (kg m^{-3})');
yyaxis(ax,'right'); scatter(ax,S_b,R.stream.mu_b*1e3,18,'filled'); ylabel(ax,'Brine viscosity \mu_b (mPa s)');
ok(ax); xlabel(ax,'Brine salinity S_b (g kg^{-1})'); title(ax,'(a) Density & viscosity');


ax=nexttile; histogram(ax,R.stream.I_s_b,25,'FaceColor',col(1,:)); ok(ax);
xlabel(ax,'Ionic strength I_s (mol kg^{-1})'); ylabel(ax,'Plant count'); title(ax,'(b) Ionic-strength distribution');


%% ===================== FIG 02 — Energy recovery pathways ================
Pmat = [R.energy.red.P_net R.energy.pro.P_net R.energy.erd.P_net R.energy.thermal.P_ORC_net]/1e3; % kW
h = figure('Name','F02_Energy','Color','w','Position',[70 70 1100 400],'Visible','on');
tiledlayout(h,1,2,'Padding','compact','TileSpacing','compact');



ax=nexttile; cr = R.stream.c_NaCl_b./max(R.stream.c_NaCl_f,eps);
scatter(ax,cr,R.exergy.eta_II_RED,22,col(2,:),'filled'); ok(ax);
xlabel(ax,'Concentration ratio c_b/c_f (-)'); ylabel(ax,'\eta_{II,RED} (-)'); title(ax,'(a) RED 2nd-law efficiency');


ax=nexttile; hyd = R.stream.Q_b.*R.stream.P_b/1e3;
scatter(ax,hyd,R.energy.erd.P_net/1e3,22,col(3,:),'filled'); ok(ax);
xlabel(ax,'Hydraulic power Q_b P_b (kW)'); ylabel(ax,'ERD net power (kW)'); title(ax,'(b) Hydraulic energy recovery');


% Technically available = pathway delivers a positive product/output where
% admissible; economically selected = included in the optimal portfolio.
np = numel(R.econ.pathway_names);
avail = zeros(1,np);
avail(1)=mean(R.energy.red.P_net>0);          avail(2)=mean(R.energy.pro.P_net>0);
avail(3)=mean(R.energy.erd.P_net>0);          avail(4)=mean(R.energy.thermal.P_ORC_net>0);
avail(5)=mean(R.energy.ec.mdot_H2>0);         avail(6)=mean(sum(R.mineral.m1.mdot_solid(:,1:2),2)>0);
avail(7)=mean(R.mineral.m2.mdot_NaCl>0);      avail(8)=mean(R.mineral.m3.mdot_Li2CO3>0);
avail(9)=mean(R.mineral.m4.mdot_Br>0);        avail(10)=mean(R.mineral.m5.Q_water_add>0);
sf = mean(R.econ.sel,1);


%% ===================== FIG 03 — Minerals & market absorption ============
mp_tech = R.mineral.mprod_yr;          % technical (chemically recoverable) mass
mp_sel  = R.econ.mprod_selected;       % economically-selected mass
pn = R.mineral.product_names;
techG = R.agg.global.mineral_Mt_yr_tech; consG = R.agg.global.mineral_Mt_yr_constr;
h = figure('Name','F03_Minerals','Color','w','Position',[80 80 600 400],'Visible','on');
tiledlayout(h,1,1,'Padding','compact','TileSpacing','compact');


ax=nexttile; idxNa=3; scatter(ax,R.stream.m_i_b(:,strcmpi(par.ions.name,'Na')).*R.stream.Q_b,mp_tech(:,idxNa)/1e3,22,col(3,:),'filled'); ok(ax);
xlabel(ax,'Brine Na^{+} molar flow (mol s^{-1})'); ylabel(ax,'NaCl production (kt yr^{-1})'); title(ax,' Halite (NaCl) yield vs Na^{+}');


%% ===================== FIG 04 — Exergy & second law ====================
dest = R.exergy.Exdot_dest/1e3;  % kW
destByReg = zeros(numel(regions),size(dest,2));
for g=1:numel(regions), m=inv.region==regions(g); destByReg(g,:)=sum(dest(m,:),1); end
h = figure('Name','F04_Exergy','Color','w','Position',[90 60 1150 780],'Visible','on');
tiledlayout(h,2,2,'Padding','compact','TileSpacing','compact');

ax=nexttile; bar(ax,destByReg,'stacked'); ok(ax); set(ax,'XTick',1:numel(regions),'XTickLabel',regions,'XTickLabelRotation',35);
ylabel(ax,'Exergy destruction (kW)'); legend(ax,R.exergy.dest_labels,'Location','south','Orientation','horizontal'); title(ax,'(a) Destruction map by region');
set(ax,'YScale','log')
ax.YScale = 'log';

ax=nexttile; etaII=[R.exergy.eta_II_RED R.exergy.eta_II_PRO R.exergy.eta_II_ERD R.exergy.eta_II_th];
% etaII(etaII<=0) = 1e-6;
boxchart(ax,etaII); 

% ax.YScale = 'log';
% 
% box(ax,'on');

ok(ax); 
xticklabels(ax,{'RED','PRO','ERD','ORC'});
ylabel(ax,'\eta_{II} (-)'); title(ax,'(b) Second-law efficiency by pathway');
% ylim([1e-1 2e0])

ax=nexttile; yyaxis(ax,'left'); scatter(ax,S_b,R.exergy.ex_ch_b/1e3,18,'filled'); ylabel(ax,'Chemical exergy (kJ kg^{-1})');
yyaxis(ax,'right'); scatter(ax,S_b,R.exergy.ex_ph_b/1e3,18,'filled'); ylabel(ax,'Physical exergy (kJ kg^{-1})');
ok(ax); xlabel(ax,'Brine salinity S_b (g kg^{-1})'); title(ax,'(c) Brine specific exergy');

ax=nexttile; yyaxis(ax,'left'); scatter(ax,R.stream.Q_b,R.exergy.Wdot_min_sep/1e6,20,'filled');
ylabel(ax,'Chemical (SGE) exergy rate (MW)');
yyaxis(ax,'right'); scatter(ax,R.stream.Q_b,R.exergy.Exdot_b/1e6,20,'filled');
ylabel(ax,'Physical exergy rate (MW)'); ok(ax);
xlabel(ax,'Brine volumetric flow Q_b (m^3 s^{-1})'); title(ax,'(d) Brine exergy resource vs throughput');


%% ===================== FIG 05 — First/second-law efficiency ============
h = figure('Name','F05_Efficiency','Color','w','Position',[60 90 1100 760],'Visible','on');
tiledlayout(h,2,2,'Padding','compact','TileSpacing','compact');

isSWRO = strcmpi(inv.tech,'SWRO'); isTh = R.energy.thermal.T_b > R.stream.T0 + 0.1;
f_off_sel = 100*R.econ.E_net_sel_yr ./ max(R.massE.E_dem_yr,eps);   % selected-portfolio offset
f_off_disp = min(max(f_off_sel,-200),100);                          % clip rare outliers for display
ax=nexttile; histogram(ax,f_off_disp,25,'FaceColor',col(1,:)); 
ax.YScale = 'log';
box(ax,'on');
ylim([1e-2 300])
ok(ax); xline(ax,0,'r--','LineWidth',1.2);
xlabel(ax,'Selected energy-offset fraction (%, clipped \pm)'); ylabel(ax,'Plant count'); title(ax,'(a) First-law energy offset');

% Hydraulic-exergy recovery efficiency (SWRO): ERD net power / brine pressure exergy
P_rec_W = R.econ.E_rec_sel_yr.*1000./max(8760*par.stream.f_av,eps);     % W (selected recovery)
Ex_brine_W = R.exergy.Wdot_min_sep + R.exergy.Exdot_b;                  % chemical + physical
eta_ex = P_rec_W ./ max(Ex_brine_W,eps);
ax=nexttile; scatter(ax,S_b(isSWRO),100*eta_ex(isSWRO),20,inv.p_elec(isSWRO),'filled'); ok(ax); colormap(ax,parula);
cb=colorbar(ax); cb.Label.String='p_{elec} (USD/kWh)';
xlabel(ax,'Brine salinity S_b (g kg^{-1})'); ylabel(ax,'Exergetic recovery efficiency (%)'); title(ax,'(b) Exergetic recovery (SWRO plants)');

ax=nexttile; carnot = 1 - R.stream.T0./max(R.energy.thermal.T_b,R.stream.T0+1e-6);
if any(isTh)
    scatter(ax,R.stream.Q_b(isTh),R.energy.thermal.P_ORC_net(isTh)/1e3,26,S_b(isTh),'filled');
    cb=colorbar(ax); cb.Label.String='S_b (g kg^{-1})';
end
ax.YScale = 'log';
box(ax,'on');
ok(ax); xlabel(ax,'Brine flow Q_b (m^3 s^{-1})'); ylabel(ax,'ORC net power (kW)'); title(ax,'(c) Thermal recovery (MED/MSF plants)');

ax=nexttile; yyaxis(ax,'left'); histogram(ax,R.lca.EROI,25);
ax.YScale = 'log';
box(ax,'on'); 
ylabel(ax,'Plant count');ylim([1e-2 300])
yyaxis(ax,'right'); [f,x]=ecdf(R.lca.EROI); stairs(ax,x,f,'LineWidth',1.6); ylabel(ax,'CDF (-)');


ok(ax); xlabel(ax,'EROI (-)'); title(ax,'(d) Energy return on investment');


%% ===================== FIG 06 — Dimensionless analysis =================
rho=R.stream.rho_b; mu=R.stream.mu_b;
Re  = rho.*par.red.v.*par.red.d_h./mu;                  % Reynolds (RED/PRO channel)
Sc  = mu./(rho.*par.pro.D_s);                           % Schmidt
Pe  = Re.*Sc;                                           % Peclet (mass)
CR  = R.stream.c_NaCl_b./max(R.stream.c_NaCl_f,eps);    % concentration ratio (varies with r)
Pi_ratio = R.stream.pi_b./max(R.stream.pi_f,eps);       % osmotic-pressure ratio brine/feed
Ex_no = R.exergy.ex_ch_b./max(R.exergy.ex_ph_b+R.exergy.ex_ch_b,eps); % chemical-exergy fraction
SGE_no = R.energy.red.P_RED_gross./max(R.energy.red.dG_mix_rate,eps);  % RED extraction vs Gibbs ceiling
dLCOWstar = R.econ.dLCOW./max(inv.LCOW_base,eps);       % normalized net-LCOW shift
h = figure('Name','F06_Dimensionless','Color','w','Position',[60 60 1400 800],'Visible','on');
tiledlayout(h,2,3,'Padding','compact','TileSpacing','compact');

ax=nexttile; histogram(ax,CR,25,'FaceColor',col(1,:)); ok(ax);
xlabel(ax,'Concentration ratio c_b/c_f (-)'); ylabel(ax,'count'); title(ax,'(a) Concentration ratio');

ax=nexttile; scatter(ax,S_b,Pi_ratio,18,inv.r,'filled'); ok(ax); cb=colorbar(ax); cb.Label.String='r (-)';
xlabel(ax,'S_b (g kg^{-1})'); ylabel(ax,'\pi_b/\pi_f (-)'); title(ax,'(b) Osmotic-pressure ratio');

ax=nexttile; reg=[median(Re) median(Sc) median(Pe)]; bar(ax,reg,'FaceColor',col(3,:)); ok(ax);
set(ax,'YScale','log','XTick',1:3,'XTickLabel',{'Re','Sc','Pe'}); ylabel(ax,'value (-)');
title(ax,'(c) Transport regime (laminar, mass-limited)');

ax=nexttile; scatter(ax,CR,Ex_no,18,col(4,:),'filled'); ok(ax);
xlabel(ax,'Concentration ratio (-)'); ylabel(ax,'Chemical-exergy fraction (-)'); title(ax,'(d) Exergy number vs CR');

ax=nexttile; scatter(ax,CR,SGE_no,18,col(5,:),'filled'); ok(ax); set(ax,'YScale','log');
yline(ax,par.red.eta_II_max,'r--','LineWidth',1.2,'Label','2nd-law cap (0.40)');
xlabel(ax,'Concentration ratio (-)'); ylabel(ax,'P_{RED,gross}/\DeltaG_{mix} (-)'); title(ax,'(e) Kinetic vs thermodynamic limit');

ax=nexttile; histogram(ax,dLCOWstar,25,'FaceColor',col(6,:)); ok(ax); xline(ax,0,'r--','LineWidth',1.2);
xlabel(ax,'\DeltaLCOW^{*}=\DeltaLCOW/LCOW_{base} (-)'); ylabel(ax,'count'); title(ax,'(f) Dimensionless net-LCOW shift');


%% ===================== FIG 07 — Techno-economics =======================
% Standalone fleet CAPEX of each pathway (cost if deployed at every plant),
% so all pathways - including the technically-capable but unselected RED/PRO -
% appear; the selected-portfolio CAPEX is overlaid for comparison.
capStand = sum(R.econ.CAPEX_k, 1) / 1e9;               % standalone, B USD
capSel   = sum(R.econ.sel .* R.econ.CAPEX_k, 1) / 1e9; % selected, B USD
capLbl   = R.econ.pathway_names;
h = figure('Name','F07_Economics','Color','w','Position',[70 60 600 400],'Visible','on');
tiledlayout(h,1,1,'Padding','compact','TileSpacing','compact');


ax=nexttile; bC=bar(ax,[capStand(:) capSel(:)],'grouped'); ok(ax); set(ax,'YScale','log');
bC(1).FaceColor=[0.75 0.8 0.9]; bC(2).FaceColor=col(2,:);
set(ax,'XTick',1:numel(capLbl),'XTickLabel',capLbl,'XTickLabelRotation',35);
ylabel(ax,'Fleet CAPEX (10^{9} USD)'); legend(ax,{'Standalone (if deployed)','Selected portfolio'},'Location','northwest');
title(ax,'CAPEX by subsystem'); ylim(ax,[1e-3 max(capStand)*2]);




%% ===================== FIG 08 — Economic sensitivity ===================
h = figure('Name','F08_EconSensitivity','Color','w','Position',[80 70 1400 400],'Visible','on');
tiledlayout(h,1,3,'Padding','compact','TileSpacing','compact');

ax=nexttile; histogram(ax,R.econ.integration_premium/1e6,30,'FaceColor',col(4,:)); ok(ax);
xlabel(ax,'Integration premium (10^{6} USD)'); ylabel(ax,'Plant count'); title(ax,'(a) Integration premium');

ax=nexttile; scatter(ax,inv.p_elec,R.agg.p_elec_crossover,24,R.agg.VI,'filled'); ok(ax); colormap(ax,turbo);
hold(ax,'on'); pl=xlim(ax); plot(ax,pl,pl,'k--'); hold(ax,'off'); cb=colorbar(ax); cb.Label.String='VI (-)';
xlabel(ax,'Actual electricity price (USD kWh^{-1})'); ylabel(ax,'Break-even price (USD kWh^{-1})'); title(ax,'(b) Cost-crossover');

ax=nexttile; [~,ord]=sort(abs(mc.tornado.SRC),'ascend'); barh(ax,mc.tornado.SRC(ord),'FaceColor',col(5,:)); ok(ax);
set(ax,'YTick',1:numel(ord),'YTickLabel',mc.tornado.var_names(ord)); xlabel(ax,'Std. regression coeff. on NPV (-)');
title(ax,'(c) Global sensitivity (tornado)');


%% ===================== FIG 09 — Environmental / LCA ====================
h = figure('Name','F09_Environment','Color','w','Position',[70 60 1150 780],'Visible','on');
tiledlayout(h,2,2,'Padding','compact','TileSpacing','compact');

ax=nexttile; histogram(ax,R.lca.CO2_mitigation_t_yr/1e3,30,'FaceColor',col(1,:)); ok(ax);
xlabel(ax,'CO_2 mitigation (kt yr^{-1})'); ylabel(ax,'Plant count'); title(ax,'(a) CO_2 mitigation potential');
ax.YScale = 'log';
box(ax,'on'); ylim([1e-1 300])

ax=nexttile; scatter(ax,S_b,R.lca.CI_water,22,inv.EF_grid,'filled'); ok(ax); colormap(ax,parula);
cb=colorbar(ax); cb.Label.String='Grid EF (kgCO_2 kWh^{-1})'; xlabel(ax,'Brine salinity S_b (g kg^{-1})');
ylabel(ax,'Water carbon intensity (kgCO_2 m^{-3})'); yline(ax,0,'r--'); title(ax,'(b) Water carbon intensity');

ax=nexttile; scatter(ax,S_b,R.lca.EROI,22,col(3,:),'filled'); ok(ax);
xlabel(ax,'Brine salinity S_b (g kg^{-1})'); ylabel(ax,'EROI (-)'); title(ax,'(c) EROI vs salinity');

ax = nexttile;
yyaxis(ax,'left');scatter(ax,R.lca.TDS_discharge_yr/1e3,R.lca.salt_diverted_t_yr/1e3,18,'filled')
ylabel(ax,'Salt diverted (kt yr^{-1})');yyaxis(ax,'right')
[counts,edges] = histcounts(R.lca.TDS_discharge_yr/1e3,25);centers = edges(1:end-1)+diff(edges)/2;
bar(ax,centers,counts,1,'FaceAlpha',0.4,'EdgeColor','none')
ylabel(ax,'Plant count');xlabel(ax,'Brine TDS discharge (kt yr^{-1})')


%% ===================== FIG 10 — Uncertainty (Monte Carlo) ==============
h = figure('Name','F10_Uncertainty','Color','w','Position',[60 60 600 400],'Visible','on');
tiledlayout(h,1,1,'Padding','compact','TileSpacing','compact');
gi=@(nm) find(strcmp(mc.kpi_names,nm),1);


ax=nexttile; [~,ord]=sort(abs(mc.tornado.SRC),'ascend'); barh(ax,mc.tornado.SRC(ord),'FaceColor',col(3,:)); ok(ax);
set(ax,'YTick',1:numel(ord),'YTickLabel',mc.tornado.var_names(ord)); xlabel(ax,'SRC on NPV (-)'); title(ax,'(c) Sensitivity ranking');




%% ===================== FIG 12 — Viability mapping ======================
h = figure('Name','F12_Viability','Color','w','Position',[60 60 1400 400],'Visible','on');
tiledlayout(h,1,3,'Padding','compact','TileSpacing','compact');


ax=nexttile; scatter(ax,inv.p_elec,R.agg.VI,24,S_b,'filled'); ok(ax); colormap(ax,parula);
cb=colorbar(ax); cb.Label.String='S_b (g kg^{-1})'; xlabel(ax,'Electricity price (USD kWh^{-1})');
ylabel(ax,'Viability index (-)'); title(ax,'(a) Viability vs electricity price');

ax=nexttile; sub=[R.agg.VI_fin R.agg.VI_thermo R.agg.VI_cost R.agg.VI_env];
msub=zeros(numel(regions),4); for g=1:numel(regions), m=inv.region==regions(g); msub(g,:)=mean(sub(m,:),1); end
bar(ax,msub,'stacked'); ok(ax); set(ax,'XTick',1:numel(regions),'XTickLabel',regions,'XTickLabelRotation',35);
ylabel(ax,'Weighted sub-scores (-)'); legend(ax,{'Financial','Thermodynamic','Cost','Environmental'},'Location','south','NumColumns',2);
title(ax,'(b) Viability composition by region');

ax=nexttile; scatter(ax,S_b,R.econ.LCOW_net,24,R.agg.VI,'filled'); ok(ax); colormap(ax,turbo);
cb=colorbar(ax); cb.Label.String='VI (-)'; yline(ax,0,'r--'); xlabel(ax,'Brine salinity S_b (g kg^{-1})');
ylabel(ax,'Net-LCOW (USD m^{-3})'); title(ax,'(c) Net-LCOW vs salinity');


%% ===================== FIG 13 — Literature benchmark verification =======
h = figure('Name','F13_Benchmarks','Color','w','Position',[60 60 800 400],'Visible','on');
tiledlayout(h,1,1,'Padding','compact','TileSpacing','compact');
mval=[lb.rows.model]'; lo=[lb.rows.lo]'; hi=[lb.rows.hi]'; pass=logical([lb.rows.pass]');
bnames={lb.rows.name}; nb=numel(bnames);
% normalise each benchmark to its [lo,hi] band: 0=lo, 1=hi; model shown as marker
norm_m = (mval-lo)./max(hi-lo,eps);
ax=nexttile; hold(ax,'on');
for kk=1:nb
    patch(ax,[0 1 1 0],[kk-0.3 kk-0.3 kk+0.3 kk+0.3],[0.85 0.92 0.85],'EdgeColor',[.7 .7 .7]);
end
cols = repmat([0.2 0.6 0.2],nb,1); cols(~pass,:)=repmat([0.85 0.2 0.2],sum(~pass),1);
scatter(ax,min(max(norm_m,-0.15),1.15),1:nb,80,cols,'filled','MarkerEdgeColor','k');
xline(ax,0,'k--'); xline(ax,1,'k--'); hold(ax,'off'); ok(ax);
set(ax,'YTick',1:nb,'YTickLabel',bnames,'YDir','reverse'); ylim(ax,[0.4 nb+0.6]); xlim(ax,[-0.25 1.25]);
xlabel(ax,'Position within literature band  (0 = lower bound, 1 = upper bound)');
title(ax,sprintf('Literature benchmark verification: %d/%d within envelope',lb.n_pass,nb));


%% ---- nested coastline helper -----------------------------------------
    function draw_coastlines(ax)
    %DRAW_COASTLINES Plot world coastlines as a light-grey backdrop. Uses the
    %Mapping Toolbox 'coastlines' data if available, else a bundled fallback.
        clat = []; clon = [];
        try
            S = load('coastlines','coastlat','coastlon'); clat = S.coastlat; clon = S.coastlon;
        catch
            try, c = load('coast.mat'); clat = c.lat; clon = c.long; catch, end
        end
        if ~isempty(clat)
            plot(ax, clon, clat, '-', 'Color',[0.6 0.6 0.6], 'LineWidth',0.5, 'HandleVisibility','off');
        end
    end

end
