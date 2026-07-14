# Methodology — *Turning Brine into an Asset*
### Global Techno‑Economic, Exergetic and Environmental Assessment of Energy and Mineral Recovery from Desalination Brine

**Project:** Turning Brine into an Asset
**Author placeholder:** `<AUTHOR>`  ·  **Date placeholder:** `<DATE>`  ·  **Version:** 2.1
**Implementation:** MATLAB R2025b (see [README.md](README.md))
**Model specification:** `docs/Brine_Valorization_Mathematical_Model.md` (Sections 1–13)

---

## 1. Purpose and Scope

This document describes, in full detail, the methodology implemented in the MATLAB code base that solves and simulates the brine‑valorization mathematical model. The model couples **five energy‑recovery pathways** and **five resource‑recovery pathways** to a georeferenced inventory of desalination plants, and closes the analysis with second‑law (exergy) accounting, a net levelized‑cost‑of‑water (net‑LCOW) economic layer, a life‑cycle energy/carbon and avoided‑discharge layer, a global spatial aggregation, and a Monte‑Carlo uncertainty / global‑sensitivity layer.

| Group | Pathways |
|---|---|
| **Energy recovery** | (E1) Reverse electrodialysis **RED** · (E2) Pressure‑retarded osmosis **PRO** · (E3) Hydraulic energy recovery / pressure exchanger **ERD** · (E4) Thermal / low‑grade‑heat **ORC + MD** · (E5) Electrochemical **chlor‑alkali** |
| **Resource recovery** | (M1) Selective sequential precipitation (gypsum, brucite) · (M2) Evaporative crystallization (halite, sylvite) · (M3) Lithium recovery · (M4) Trace elements (Br, Rb) · (M5) Near‑ZLD incremental water |
| **Analysis layers** | (A1) Exergy / second law · (A2) Techno‑economics (net‑LCOW) · (A3) Life‑cycle assessment · (A4) Global resource aggregation · (A5) Uncertainty & sensitivity |

---

## 2. Software Architecture

The model is fully modular; every function receives an explicit input struct and returns an explicit output struct, and **every numeric constant lives in a single `config_parameters.m`** so the model is auditable and reproducible. There are no global variables. Plant‑level computations are **vectorized over the 300‑plant inventory**, so the Monte‑Carlo layer re‑evaluates the whole fleet per draw efficiently.

```
BrineValorization/
  main_driver.m                 % orchestration: inputs → pipeline → MC → V&V → figures → export
  config_parameters.m           % single source of truth for every constant / range
  data/     generate_inventory.m
  props/    rho_sw, mu_sw, cp_sw, h_sw, s_sw, activity_coeffs, osmotic_pressure
  streams/  M00_streams.m
  energy/   red_power, pro_power, erd_power, thermal_power, electrochem_value,
            energy_dispatch, mass_energy_integration
  minerals/ precip_sequence, crystallizer, lithium_recovery, trace_recovery,
            zld_water, mineral_dispatch
  exergy/   exergy_balance.m
  economics/techno_economics.m
  environment/ lca_environmental.m
  aggregation/ global_aggregate.m
  uncertainty/ monte_carlo.m
  core/     run_pipeline, fleet_kpis, validate_inputs, validate_results, literature_benchmarks
  viz/      make_figures, export_excel, export_results_mat
  results/  (generated outputs)
  docs/     model spec + formulas Word document
```

### Execution pipeline (`main_driver.m`)
1. `config_parameters()` → build parameter struct `par`.
2. **`validate_inputs(par)`** → 59‑point input‑data audit (V&V layer 1).
3. `generate_inventory(par)` → synthetic 300‑plant global fleet.
4. `run_pipeline(inv,par)` → deterministic per‑plant solve (Sections 2–9).
5. `monte_carlo(inv,par)` → Latin‑hypercube uncertainty + regression sensitivity (Section 10).
6. **`validate_results`** (14 internal checks, V&V layer 2) and **`literature_benchmarks`** (10 checks, V&V layer 3).
7. `make_figures` → 13 publication figures (all `tiledlayout`‑grouped, `.fig` saved).
8. `export_excel`, `export_results_mat` → workbook (15 sheets), MAT, CSV.

---

## 3. Detailed Methodology

### 3.1 Plant inventory (`generate_inventory.m`)
A synthetic, literature‑calibrated **300‑plant georeferenced inventory** across 11 desalination regions is generated with a fixed random seed. Each plant carries technology (SWRO/MED/MSF), capacity (log‑uniform 2 000–600 000 m³/day), feed salinity, recovery ratio, latitude/longitude, regional electricity price, grid‑carbon factor and baseline LCOW. This is a **documented, transparent substitute** for a licensed DesalData extract (which the model spec calls for but is unavailable), preserving the regional heterogeneity needed for meaningful global aggregation.

### 3.2 Stream characterization (`M00_streams.m`, `props/*`)
Implements Eqs. 2.1–2.14. Mass/volume balances give permeate, feed and brine flows and the salt‑rejection‑refined brine salinity. Thermophysical properties use the **Sharqawy et al. (2010)** seawater correlations for density, viscosity and specific heat; enthalpy/entropy use closed‑form integrals of the specified `cp(S,T)`. Ionic speciation (12 ions) yields molality, molar concentration and ionic strength. Activity coefficients use a **Davies–Pitzer‑lite blend** (Davies at low ionic strength, a NaCl‑dominant Pitzer reduction at brine‑level `I_s ≈ 1.3–1.6 mol/kg`); osmotic pressure follows the rigorous activity relation. Physical‑bounds assertions (`0<r<1`, `S_b>S_f`, `ρ_b>ρ_f`, `π_b>π_f>0`) act as an inline closure gate.

### 3.3 Energy‑recovery subsystems (`energy/*`)
* **RED** (`red_power`): cell‑pair EMF, area resistance, matched‑load gross power, areal power density, laminar pumping parasitic, and the **Gibbs free‑energy‑of‑mixing ceiling**. The draw is the brine; the dilute counter‑stream is a **low‑salinity treated/brackish effluent** (`par.sge.S_dilute = 8 g/kg`), consistent with the model's "low‑salinity counter‑stream". Net power is **capped at `η_II,max = 0.40 × ΔĠ_mix`** for thermodynamic admissibility (a finite stack cannot exceed a realizable fraction of the mixing ceiling).
* **PRO** (`pro_power`): implicit water flux (ICP/ECP) solved by damped vectorized fixed‑point iteration; areal & peak power density; membrane **physically sized as `A_PRO = Q_dilute/J_w`**; the same 0.40 mixing‑ceiling cap.
* **ERD** (`erd_power`): isobaric pressure‑exchanger recovery of residual SWRO brine pressure (5.5 MPa).
* **Thermal** (`thermal_power`): sensible availability, Carnot‑limited and ORC net work, membrane‑distillation distillate — only where a temperature lift exists (MED/MSF).
* **Electrochemical** (`electrochem_value`): chlor‑alkali cell voltage, H₂/Cl₂/NaOH production and net priced value.
* **Dispatcher** (`energy_dispatch`): applies admissibility masks and rolls up recovered energy.

### 3.4 Resource‑recovery subsystems (`minerals/*`)
* **M1 precipitation** (`precip_sequence`): saturation‑index‑driven **gypsum (CaSO₄·2H₂O)** and **brucite (Mg(OH)₂)** from raw brine; Mg(OH)₂ is precipitated with **lime Ca(OH)₂** (industrial practice), gypsum needs no reagent.
* **M2 crystallization** (`crystallizer`): evaporative concentration to halite saturation recovers **NaCl (Na‑limited)** and, from the residual liquor, **KCl (K‑limited)**; carries the evaporation energy split between low‑grade heat and electricity. Because desal brine is ~83 % NaCl by recovered‑salt mass, halite is the single largest technical resource.
* **M3 lithium** (`lithium_recovery`): Langmuir uptake, recovered Li → battery‑grade Li₂CO₃, with a capital‑intensive DLE cost and processing OPEX.
* **M4 trace** (`trace_recovery`): conservative Br/Rb accounting with element‑specific efficiency (Rb down‑rated to a bulk‑equivalent value).
* **M5 near‑ZLD** (`zld_water`): incremental RO/MD water recovery and residual‑salinity uplift.
* **Dispatcher** (`mineral_dispatch`): assembles the 7‑product mass vector, energy demand, reagent cost and gross revenue.

### 3.5 Exergy / second‑law (`exergy_balance.m`)
Implements Eqs. 5.1–5.13 against the local‑seawater dead state. Physical exergy uses the enthalpy/entropy departures; **chemical (salinity‑gradient) exergy is the positive brine‑vs‑seawater Gibbs‑mixing exergy** (the physically correct, sign‑consistent form). Subsystem second‑law efficiencies and a **Gouy–Stodola** exergy‑destruction map (per subsystem, per region) are produced.

### 3.6 Mass–energy integration (`mass_energy_integration.m`)
Eqs. 6.1–6.5: primary desalination energy demand, recovered/consumed valorization energy, net plant energy impact and the first‑law offset fraction, with a machine‑precision mass‑balance closure gate.

### 3.7 Techno‑economics (`techno_economics.m`)
Implements Eqs. 7.1–7.18 with the **central economic innovation — marginal per‑pathway portfolio selection**: each of the ten pathways is evaluated on its own attributable CAPEX (six‑tenths scaling, CEPCI‑ and installation‑adjusted), revenue and OPEX, annualized via the capital‑recovery factor, and **included in a plant's retrofit only if its standalone marginal NPV > 0**. Reported plant net‑LCOW, NPV, IRR (vectorized Newton–Raphson), payback and BCR therefore describe the **optimal valorization portfolio**; the difference between portfolio NPV and the best single‑pathway NPV defines the **integration premium**. RED/PRO stacks are **right‑sized to the gradient‑limited power** (`A_econ = P_capped/P_d`) and costed at recalibrated literature stack prices (~250–300 USD/m²).

### 3.8 Life‑cycle assessment (`lca_environmental.m`)
Eqs. 8.1–8.13: operational GWP from the selected‑portfolio electricity balance and regional grid factor; embodied GWP from membranes and process energy; **avoided‑burden credits** for displaced virgin mineral and water production; net GWP, water carbon intensity, EROI, and the **avoided brine‑discharge salinity‑load** indicator.

### 3.9 Global aggregation (`global_aggregate.m`)
Eqs. 9.1–9.7 with a deliberate **two‑tier logic**: (i) **plant‑level viability** is assessed at **merchant prices** (a single plant is a price‑taker), producing NPV, net‑LCOW, the electricity‑price cost‑crossover, and a bounded **composite viability index** (financial + thermodynamic + cost + environmental sub‑scores); (ii) the **global resource potential** is reported both as technical potential and as **market‑absorption‑constrained** potential, where each product's global recovered mass is capped at its addressable global‑market share. The ratio exposes which products are demand‑limited (bulk NaCl, Mg(OH)₂, Br) versus chemistry/energy‑limited (Li, KCl, gypsum, Rb).

### 3.10 Uncertainty & sensitivity (`monte_carlo.m`)
Eqs. 10.1–10.3: a **Latin‑hypercube** ensemble (default N = 600) over seven uncertain inputs mapped through their marginal distributions (lognormal/triangular/uniform via inverse‑CDF, with a bundled Acklam `norminv`), re‑running the full deterministic pipeline per draw. Outputs are summarized (P5/P50/P95, CoV) and **standardized regression coefficients (SRC)** provide the global‑sensitivity tornado (with the linear‑model R² reported so the user can judge SRC adequacy; SRC² ≈ first‑order Sobol for near‑linear response).

---

## 4. Verification and Validation (three independent layers)

| Layer | Module | Checks | Result |
|---|---|---|---|
| **Input data** | `validate_inputs.m` | 59 — efficiencies ∈ (0,1], prices > 0, literature ranges, ion‑table closure (Σw ≈ 35), electroneutrality, viability‑weight normalization, distribution ordering, cross‑checks | **59/59 PASS** |
| **Internal consistency** | `validate_results.m` | 14 — mass closure, salinity/density/osmotic ordering, finiteness, second‑law admissibility (`η_II ≤ 1`, RED ≤ 0.40 ceiling), non‑negative exergy destruction, selection consistency, net‑LCOW identity, VI bounds, MC percentile ordering | **14/14 PASS** |
| **Literature benchmarks** | `literature_benchmarks.m` | 10 — RED/PRO power density, brine osmotic pressure, SWRO & ZLD specific energy, Li capacity, RED η_II, ERD specific recovery, high‑salinity viability, NaCl share | **10/10 within envelope** |

---

## 5. Key Modeling Decisions (documented deviations from the spec)

1. **Marginal portfolio selection** replaces "capitalize every pathway everywhere," which otherwise buries viable high‑value recovery under uneconomic bulk‑evaporation cost.
2. **Two‑tier market logic** — per‑plant merchant viability vs a *global* market‑absorption ceiling — separates plant bankability from world‑demand limits.
3. **Thermodynamic admissibility caps** on RED/PRO (`η_II ≤ 0.40 × ΔĠ_mix`) prevent second‑law violations.
4. **Low‑salinity effluent counter‑stream** (8 g/kg) for RED/PRO, matching literature power densities.
5. **Lime (not caustic)** for Mg(OH)₂; **evaporative‑crystallizer coupling** for NaCl/KCl.
6. **Right‑sized RED/PRO membranes** + literature stack costs (fixes an 5–9× area over‑sizing).
7. **Reduced Monte‑Carlo/Sobol sample sizes** (N = 600) vs the spec's 10 000 for tractable runtime — documented in `config_parameters.m`.
8. The plant inventory is a **synthetic, documented substitute** for licensed DesalData.

---

## 6. Generated Figures (13 grouped `tiledlayout` figures)

All figures are saved as MATLAB `.fig` in `results/figures/`; their numeric data + metadata are exported one‑sheet‑per‑figure to `results/BrineValorization_Figures.xlsx`.

| # | Name / Title | Caption | Description (panels) | Data interpretation |
|---|---|---|---|---|
| 1 | **F01_Streams** — Stream Characterization & Thermophysical Properties | Thermophysical and chemical characterization of desalination brine across the 300‑plant global inventory. | (a) Brine vs feed salinity coloured by recovery (Eq. 2.6); (b) Sharqawy density & viscosity vs brine salinity, dual axes (Eqs. 2.7–2.8); (c) brine osmotic pressure coloured by feed osmotic pressure (Eq. 2.13); (d) ionic‑strength distribution (Eq. 2.11). | Brine salinity rises steeply with recovery to 60–80 g/kg; density and viscosity increase quasi‑linearly, raising pumping parasitics; osmotic pressure spans 45–70 bar, bounding the gradient and minimum separation work; ionic strength clusters at 1.2–1.6 mol/kg, the Pitzer‑lite regime. |
| 2 | **F02_Energy** — Energy‑Recovery Pathways & Portfolio Selection | Salinity‑gradient, hydraulic and thermal energy‑recovery performance and economic‑portfolio selection. | (a) Net‑power box‑plots RED/PRO/ERD/ORC; (b) RED η_II vs concentration ratio; (c) ERD net power vs residual hydraulic power; (d) technically‑available vs economically‑selected fraction per pathway. | RED/PRO generate net power within the 0.40 second‑law ceiling but, at literature stack costs, have negative marginal NPV so are **not selected**; the light‑vs‑dark bar gap shows salinity‑gradient power is not yet competitive with cheap, proven ERD (selected at ~every SWRO plant); ORC only at the ~26 % thermal plants. |
| 3 | **F03_Minerals** — Mineral Recovery & Market‑Absorption Constraint | Mineral‑recovery potential and the market‑absorption constraint separating technical from deployable yield. | (a) Fleet technical vs selected mass per product; (b) halite yield vs brine Na⁺ flow; (c) global technical vs market‑absorbable potential (log); (d) absorbable fraction `f_absorb`. | Halite (NaCl) is by far the largest technical resource (~83 % of recovered salt), then Mg(OH)₂, gypsum, Br; yet the world market absorbs only a few percent of NaCl/Mg(OH)₂/Br while Li, KCl, gypsum, Rb are fully absorbable — **demand, not chemistry, is the binding global constraint**. |
| 4 | **F04_Exergy** — Exergy & Second‑Law Analysis | Exergy accounting and second‑law performance (Gouy–Stodola framework). | (a) Regional exergy‑destruction map by subsystem; (b) η_II by pathway; (c) chemical vs physical specific exergy; (d) chemical (SGE) & physical exergy rate vs throughput. | Destruction concentrates in mineral processing and the high‑throughput Gulf/Red‑Sea regions; ERD is most exergy‑efficient (~0.9), RED bounded at 0.40; brine chemical (salinity) exergy is the principal recoverable resource for SWRO, physical (thermal) for MED/MSF. |
| 5 | **F05_Efficiency** — First‑ & Second‑Law Efficiency Analysis | First‑ and second‑law efficiency indicators of the integrated retrofit. | (a) Selected first‑law energy‑offset distribution; (b) exergetic recovery efficiency vs salinity (SWRO); (c) ORC net power vs brine flow (MED/MSF); (d) EROI histogram + CDF. | The selected portfolio is a net energy consumer at many plants (energy‑intensive salt/Li trains) offset by ERD credit; exergetic recovery efficiency rises with salinity; thermal recovery scales with throughput; EROI is bimodal — <1 where mineral processing dominates, ≫1 where high‑value recovery leads. |
| 6 | **F06_Dimensionless** — Dimensionless Analysis | Dimensionless characterization of transport, thermodynamic and economic behaviour. | (a) Concentration ratio; (b) osmotic‑pressure ratio vs salinity; (c) Re/Sc/Pe regime (log); (d) chemical‑exergy fraction vs CR; (e) kinetic RED power / Gibbs ceiling with 0.40 cap; (f) normalized net‑LCOW shift. | Channel flow is laminar (Re≈3), very high Schmidt (~600) and mass‑transfer‑limited (Pe≈2000) → concentration polarization governs RED/PRO; the kinetic RED power exceeds the Gibbs ceiling 3–4×, so the **second law (not kinetics) limits extraction**; the dimensionless net‑LCOW shift is predominantly negative. |
| 7 | **F07_Economics** — Techno‑Economic Performance | Techno‑economic performance of the optimal valorization portfolio (merchant prices). | (a) Net‑LCOW distribution vs break‑even; (b) CAPEX by subsystem — standalone vs selected (log); (c) NPV vs CAPEX coloured by viability; (d) IRR histogram + payback CDF. | Net‑LCOW is negative for most plants (fleet‑mean ≈ −1.4 USD/m³) — water becomes a profitable by‑product; right‑sized RED/PRO stacks are an order of magnitude cheaper than the mineral trains yet still excluded (revenue too small); median IRR ≈ 31 %, ~3‑yr payback. |
| 8 | **F08_EconSensitivity** — Economic Robustness & Sensitivity | Benefit‑cost ratio, integration premium, cost‑crossover and global sensitivity. | (a) BCR distribution; (b) integration premium; (c) break‑even vs actual electricity price coloured by viability; (d) SRC tornado for global NPV. | BCR clusters above unity (~1.7); positive integration premium confirms co‑deployment synergy; almost all plants sit below the crossover line (break‑even price exceeds the price paid); the discount rate dominates NPV sensitivity, then electricity price. |
| 9 | **F09_Environment** — Environmental & Life‑Cycle Assessment | Life‑cycle environmental performance and avoided‑discharge benefit. | (a) CO₂‑mitigation distribution; (b) water carbon intensity vs salinity coloured by grid factor; (c) EROI vs salinity; (d) salt diverted vs brine TDS load + avoided‑fraction histogram. | Most plants achieve net CO₂ mitigation via avoided virgin production; many reach negative water carbon intensity (a climate credit) where grid carbon is high; several percent of the salinity load is diverted from discharge — a co‑benefit beyond carbon and revenue. |
| 10 | **F10_Uncertainty** — Uncertainty Propagation & Global Sensitivity | Monte‑Carlo uncertainty propagation and global sensitivity of headline outcomes. | LHS ensemble (N=600). (a) Global‑NPV distribution P5/P50/P95; (b) fleet‑mean net‑LCOW CDF; (c) SRC ranking; (d) joint NPV–CO₂ ensemble coloured by viability. | Global NPV is strongly positive across the entire ensemble and net‑LCOW lies almost entirely below zero — the profitability conclusion is robust to parameter uncertainty; discount rate and electricity price govern the spread; NPV and CO₂ mitigation are positively associated (economic and climate goals aligned). |
| 11 | **F11_Atlas** — Global Resource Atlas & Regional Potential | Georeferenced viability atlas and regional resource potential. | (a) Global map of plant viability with world **coastlines**; (b,c) region‑aggregated NPV and recovered‑energy potential scaled to the global fleet. | High‑viability clusters coincide with high‑salinity, high‑throughput Gulf/Red‑Sea SWRO fleets and high‑grid‑carbon regions; the Arabian Gulf dominates both NPV and recoverable‑energy potential — the priority deployment geography. |
| 12 | **F12_Viability** — Viability Mapping & Drivers | Viability mapping: drivers of plant‑level attractiveness. | (a,b) Viability index vs salinity and vs electricity price; (c) regional sub‑score decomposition (financial/thermodynamic/cost/environmental); (d) net‑LCOW vs salinity coloured by viability. | Viability rises with brine salinity and falls with electricity price; financial and cost sub‑scores dominate, the environmental sub‑score lifts high‑grid‑carbon regions; the most concentrated brines are simultaneously the most profitable to valorize. |
| 13 | **F13_Benchmarks** — Literature Benchmark Verification | Verification of model outputs against published literature envelopes (model Sec. 11.5). | Each row = a benchmark; the green band is its literature range normalised to [0,1]; the marker is the model value (green = inside, red = outside). | All 10/10 benchmarks fall within their literature envelopes — RED 3.8 W/m², PRO 5.5 W/m², brine π 5.73 MPa at 70 g/kg, NaCl‑dominant (83 %) recovered‑salt spectrum — externally validating the calibrated model. |

---

## 7. Generated Tables and Output Files

Data are exported to a single Excel workbook (one sheet per figure + two audit sheets), a MAT archive and a CSV.

| Table / sheet | Name | Caption | Description | Data interpretation |
|---|---|---|---|---|
| Excel `Figure_001…013` | Per‑figure data sheets | Numeric data behind each figure. | Each sheet leads with the figure **name, title, caption, explanation, interpretation**, followed by every plotted data block (per‑plant series, regional aggregates, global potentials, MC samples, SRC matrix, etc.). | Provides full reproducibility of every plotted point; column headers carry units. |
| Excel `Input_Validation` | Input‑data audit | 59‑parameter admissibility & literature‑range audit. | Columns: parameter, status (PASS/WARN/FAIL), value, admissible range. | Confirms every input is physically admissible and literature‑consistent before the model runs (59/59 PASS). |
| Excel `Literature_Benchmarks` | Output verification | Model outputs vs published envelopes. | Columns: benchmark, model value, lower bound, upper bound, unit, pass. | Confirms 10/10 computed quantities fall within literature ranges. |
| `BrineValorization_plants.csv` | Per‑plant headline results | Compact per‑plant result table. | id, region, tech, S_f, S_b, r, CAPEX, OPEX, revenue, NPV, IRR, payback, BCR, net‑LCOW, ΔLCOW, CO₂ mitigation, EROI, viability index. | One row per plant for external analysis / mapping / statistics. |
| `BrineValorization_Results.mat` | Full reproducible archive | Complete result structs (`-v7.3`). | `inventory`, `pipeline` (all module outputs), `montecarlo`, `validation`, `lit_benchmarks`, `input_validation`, `parameters`, `timestamp`. | Enables regeneration of every figure/table without re‑running the model. |
| `results/figures/*.fig` | 13 MATLAB figures | Editable vector figures. | One `.fig` per grouped figure (F01–F13). | Publication‑ready, re‑openable and editable in MATLAB. |

### Console V&V reports
Running `main_driver.m` prints three verification reports to the Command Window — **Input Data Validation (59)**, **Validation Report (14)** and **Literature Benchmark Verification (10)** — plus a headline‑results summary.

---

## 8. Reference Basis
Sharqawy, Lienhard & Zubair, *Desalin. Water Treat.* **16** (2010) 354–380 (seawater properties). Pitzer, *J. Phys. Chem.* **77** (1973) 268 (ion‑interaction model). Salinity‑gradient‑power, ZLD and selective‑recovery benchmarks per the model concept sources. Full equation set with native Word equations: `docs/Brine_Valorization_Formulas.docx`.
