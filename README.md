# Turning Brine into an Asset — MATLAB Model

Global techno‑economic, exergetic and environmental assessment of energy and mineral recovery from desalination brine. This repository solves and simulates the full mathematical model in `docs/Brine_Valorization_Mathematical_Model.md`, produces 13 publication‑quality figures with one‑sheet‑per‑figure data export, and runs three independent verification & validation layers.

> **Documentation:** see [`METHODOLOGY.md`](METHODOLOGY.md) for the full methodology and the complete figure/table catalogue, and `docs/Brine_Valorization_Formulas.docx` for every equation in native Word format.

---

## 1. Software Environment

| Item | Requirement |
|---|---|
| **MATLAB version** | Developed and verified on **MATLAB R2025b** (25.2). Requires **R2020b or newer** (uses `tiledlayout`, `boxchart`, `xline`/`yline`, `exportgraphics`, `turbo`). R2024a+ recommended per the project standard. |
| **Operating system** | Developed and tested on **Windows 11 Pro (64‑bit)**. Platform‑independent MATLAB code — also runs on macOS and Linux. The only OS‑specific optional step is opening the Word formulas file. |
| **Hardware** | Any modern CPU; < 1 GB RAM. Full run (deterministic + 600‑draw Monte Carlo + 13 figures + exports) completes in roughly **1–2 minutes**. |

### Required MATLAB toolboxes

| Toolbox | Required? | Used for |
|---|---|---|
| **MATLAB** (base) | **Yes** | Core computation, `tiledlayout`/`boxchart` figures, `writetable`/`writecell` Excel export, `savefig`. |
| **Statistics and Machine Learning Toolbox** | **Yes** | `prctile` (percentiles), `ecdf` (CDF plots), `randsample` (inventory sampling). |
| **Mapping Toolbox** | **Optional** | World coastlines on the Figure 11 atlas (`load coastlines`). If absent, the code degrades gracefully — the map is drawn without coastlines (a `try/catch` fallback is built in). |

The model deliberately avoids other toolbox dependencies: Latin‑hypercube sampling and the normal inverse‑CDF are implemented locally in `uncertainty/monte_carlo.m`, so no Global Optimization or additional toolboxes are needed.

> Check your installation with `ver` (or `matlab.addons.installedAddons`). The Statistics and Machine Learning Toolbox is the only hard third‑party dependency.

---

## 2. How to Run

1. **Open MATLAB** and set the current folder to the project root (the folder containing `main_driver.m`), e.g.
   ```matlab
   cd 'D:\MATLAB_codes\ClaudeCodeWS\BrineValorization'
   ```
2. **Run the driver.** It adds all subfolders to the path automatically, so no manual `addpath` is needed:
   ```matlab
   main_driver
   ```
3. **Watch the Command Window.** Three V&V reports print in sequence — **Input Data Validation (59 checks)**, **Validation Report (14 checks)** and **Literature Benchmark Verification (10 checks)** — followed by a headline‑results summary and the export confirmations.
4. **Collect the outputs** from the `results/` folder (see §4).

### Running individual components
Every module is a standalone function and can be called directly for inspection, for example:
```matlab
par = config_parameters();          % all constants and ranges
iv  = validate_inputs(par);         % 59-point input audit
inv = generate_inventory(par);      % 300-plant synthetic fleet
R   = run_pipeline(inv, par);       % full deterministic solve (all module outputs in R)
k   = fleet_kpis(inv, R);           % scalar headline KPIs
lb  = literature_benchmarks(inv, R, par);   % 10 literature checks
```

### Adjusting the run
* **Reproducibility:** random seeds are fixed in `config_parameters.m` (`par.inv.seed`, `par.unc.seed`); results are identical run‑to‑run.
* **Monte‑Carlo size:** change `par.unc.N_MC` (default 600) in `config_parameters.m` — larger = smoother distributions, longer runtime.
* **Any physical/economic assumption:** edit the single `config_parameters.m`; no constant is hard‑coded elsewhere.

---

## 3. Input Files

The model is **self‑contained** — it generates its own inputs from parameters; no external data files are required to run.

| Input | Type | Description |
|---|---|---|
| `config_parameters.m` | MATLAB function | **The master input.** Every constant, default value, cost, efficiency, price, emission factor and uncertainty range (model Sections 1–10). Single source of truth. |
| `data/generate_inventory.m` | MATLAB function | Generates the synthetic **300‑plant georeferenced inventory** (a documented, literature‑calibrated substitute for a licensed DesalData extract). Controlled by `par.inv.*`. |
| `docs/Brine_Valorization_Mathematical_Model.md` | Reference | The model specification the code implements (not read at runtime). |

To use a **real plant inventory**, replace `generate_inventory.m` with a function returning a table with the same columns (`id, region, tech, C_des, S_f, r, lat, lon, p_elec, EF_grid, LCOW_base`); the rest of the pipeline is unchanged.

---

## 4. Output Files

All outputs are written to `results/`.

| File | Description |
|---|---|
| `BrineValorization_Figures.xlsx` | **15 sheets.** `Figure_001`…`Figure_013` — one per figure, each with the figure name/title/caption/explanation/interpretation followed by all plotted numeric data blocks. Plus `Input_Validation` (59‑parameter audit) and `Literature_Benchmarks` (10 checks). |
| `BrineValorization_Results.mat` | Full reproducible archive (`-v7.3`): inventory, all pipeline structs, Monte‑Carlo ensemble, all three V&V structs, parameters, timestamp. |
| `BrineValorization_plants.csv` | Compact per‑plant headline table (NPV, IRR, payback, BCR, net‑LCOW, CO₂ mitigation, EROI, viability index, …). |
| `figures/F01…F13.fig` | 13 editable MATLAB figures (grouped `tiledlayout`, publication styling). |

Console output additionally prints the three verification reports and the headline‑results block.

---

## 5. Project Structure

```
BrineValorization/
├── main_driver.m              # run this
├── config_parameters.m        # all inputs / constants
├── README.md   METHODOLOGY.md
├── data/         # synthetic plant inventory
├── props/        # seawater properties + activity model
├── streams/      # stream characterization
├── energy/       # RED, PRO, ERD, thermal, electrochemical + dispatch
├── minerals/     # M1–M5 recovery + dispatch
├── exergy/       # second-law / exergy balance
├── economics/    # net-LCOW techno-economics + portfolio selection
├── environment/  # life-cycle assessment + avoided discharge
├── aggregation/  # global aggregation + viability + market cap
├── uncertainty/  # Latin-hypercube Monte Carlo + sensitivity
├── core/         # pipeline, KPIs, input/output/literature validation
├── viz/          # figures + Excel/MAT export
├── results/      # generated outputs (figures, xlsx, mat, csv)
└── docs/         # model spec + native-Word formulas document
```

---

## 6. Verification & Validation

Every run executes three independent V&V layers (all currently passing):

| Layer | Module | Result |
|---|---|---|
| Input‑data audit | `core/validate_inputs.m` | **59/59 PASS** |
| Internal consistency | `core/validate_results.m` | **14/14 PASS** |
| Literature benchmarks | `core/literature_benchmarks.m` | **10/10 within envelope** |

---

## 7. Notes and Assumptions
* The plant inventory is a **synthetic, documented substitute** for a licensed DesalData extract; all downstream global aggregates inherit this caveat and are stated as order‑of‑magnitude potentials.
* RED/PRO are shown to be **technically capable but not yet economically selected** from desalination brine — a deliberate, literature‑consistent result, not a modeling error.
* Cost references, recovery efficiencies and the RED/PRO stack prices are documented engineering assumptions in `config_parameters.m`; substitute site‑specific values there as needed.

---

## 8. Contact

**Author:** `<AUTHOR>`  ·  Prof. Hossam S. S. AbdelMeguid
**Affiliation:** Thermo‑fluid Engineering, Water Systems and Sustainable Energy
**Email:** `hssaleh@mans.edu.eg`
**Date:** `<DATE>`  ·  **Version:** 2.1

---

## 9. Citation

If you use this model, please cite the accompanying study:

> *Turning brine into an asset: a global techno‑economic and environmental assessment of energy and mineral recovery from desalination brine.* `<AUTHOR>`, `<YEAR>`.

*Generated with MATLAB R2025b. Model specification: `docs/Brine_Valorization_Mathematical_Model.md`.*
