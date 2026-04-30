# Surface Segregation and Precipitation-Driven Self-Healing in Fe–X Alloys

This repository contains the Python code, input datasets and derived numerical outputs used to reproduce the surface-segregation calculations and robustness analyses reported in the manuscript:

**A Thermodynamic Framework for Precipitation-Driven Self-Healing: Linking Surface Segregation and Precipitate Nucleation at Cavity Surfaces in Fe–X Alloys**

The calculations evaluate equilibrium surface segregation in Fe–Au, Fe–Cu, Fe–Mo and Fe–W alloys using a multilayer thermodynamic model coupled to Monte Carlo minimisation. The model outputs layer-resolved composition, adsorption/surface excess, surface energy and robustness metrics used to interpret surface-conditioned nucleation at cavity surfaces.

---

## Repository structure

```text
repository/
│
├── README.md
├── requirements.txt
├── LICENSE
├── .gitignore
│
├── src/
│   ├── seg_model.py
│   └── single_mc.py
│
├── scripts/
│   ├── sweep_A_run.py
│   └── run_robustness_tests.py
│
├── data/
│   ├── thermodynamic/
│   │   └── steel_database_fix.tdb
│   │
│   └── surface_energy/
│       └── surfaces.json
│
├── outputs/
│   ├── main_sweep/
│   ├── robustness_outputs/
│   └── A_sweep_outputs/
│
├── figures/
│   └── manuscript_figures/
│
└── docs/
    └── code_inventory.md
```

---

## Code components

### Main model

`src/seg_model.py`

Analytical multilayer segregation model coupled to Monte Carlo minimisation. It performs temperature and crystallographic-plane sweeps for Fe–X alloys and exports wide, long and optional Monte Carlo trace datasets.

Main functions:
- evaluate binary Fe–X thermodynamic quantities;
- compute pure-element and alloy surface-energy terms;
- minimise the multilayer surface free energy;
- export layer composition, adsorption and surface-energy outputs;
- optionally include an elastic screening contribution.

### Single-condition model

`src/single_mc.py`

Reusable single-case implementation of the same segregation model. It is designed for robustness testing where alloy, plane, temperature, seed, number of layers, interaction-parameter scaling and elastic scaling are varied independently.

Main entry point:
- `run_single_case(...)`

### Elastic-factor sweep

`scripts/sweep_A_run.py`

Driver script for sweeping the elastic scaling factor `A`.

Default behaviour:
- evaluates `A` values from `A_MIN` to `A_MAX`;
- optionally runs a no-elastic baseline;
- parallelises independent cases;
- writes results to `A_sweep_outputs/`;
- writes a manifest file describing generated outputs.

### Robustness testing

`scripts/run_robustness_tests.py`

Driver script for robustness analyses.

Implemented tests:
1. seed robustness;
2. number-of-layers robustness (`L`);
3. interaction-parameter robustness (`omega_scale`).

Default outputs:
- `robustness_results__long.csv`;
- `robustness_results__summary.csv`;
- optional representative Monte Carlo trace files.

---

## Input datasets

### CALPHAD thermodynamic database

`data/thermodynamic/steel_database_fix.tdb`

This file is the MatCalc Fe-based thermodynamic database used for CALPHAD calculations. It is required by `seg_model.py` and `single_mc.py` to evaluate thermodynamic quantities for the Fe–X binary systems.

Source:

```text
MatCalc open databases
https://www.matcalc.at/index.php/databases/open-databases
```

Before redistributing this file in a public repository, check the applicable MatCalc database licence. If redistribution is not permitted, do not include the `.tdb` file directly. Instead, instruct users to download it from MatCalc and place it at:

```text
data/thermodynamic/steel_database_fix.tdb
```

### Pure-element surface-energy dataset

`data/surface_energy/surfaces.json`

This file contains pure-element surface-energy data used to define elemental surface-energy terms in the segregation model. The data are derived from the surface-energy database reported by:

```text
Tran, R., Xu, Z., Radhakrishnan, B., Winston, D., Sun, W., Persson, K. A. & Ong, S. P.
Surface energies of elemental crystals.
Scientific Data 3, 160080 (2016).
https://doi.org/10.1038/sdata.2016.80
```

The database reports calculated surface energies for elemental crystals, including facet-resolved values, Wulff-weighted surface energies, surface anisotropy and shape factors. The complete dataset was made available as JSON files through the associated data repository.

---

## Model scope

The repository supports calculations for the following binary Fe–X systems:

| Alloy label in code | Solute | Nominal solute content |
|---|---:|---:|
| `Fe-Au_wt` | Au | 2.87 wt.% |
| `Fe-1.15Cu_wt` | Cu | 1.15 wt.% |
| `Fe-6.207Mo_wt` | Mo | 6.207 wt.% |
| `Fe-3.8W_wt` | W | 3.8 wt.% |

Default crystallographic planes:

```text
{100}, {110}, {111}
```

Default temperature range:

```text
400–700 °C, in 25 °C increments
```

---

## Main outputs

The scripts generate CSV datasets containing:

- surface energy, `gamma_Jm2`;
- adsorption/surface excess, `ads_monolayers`;
- first-two-layer adsorption, `ads_first2_monolayers`;
- layer-resolved solute concentrations;
- Monte Carlo convergence traces, when enabled;
- sensitivity metrics for `A`, `L`, seed and `omega_scale`.

Recommended output organisation:

```text
outputs/
├── main_sweep/
│   ├── results_wide.csv
│   ├── results_long.csv
│   └── mc_trace.csv
│
├── robustness_outputs/
│   ├── robustness_results__long.csv
│   ├── robustness_results__summary.csv
│   └── trace__*.csv
│
└── A_sweep_outputs/
    ├── manifest.csv
    ├── EL1__A*__wide.csv
    ├── EL1__A*__long.csv
    ├── EL1__A*__trace.csv
    ├── EL0__A*__wide.csv
    └── EL0__A*__long.csv
```

Large trace files should only be retained if they are required to reproduce convergence figures or representative robustness checks.

---

## Relationship to manuscript figures and tables

| Manuscript item | Calculation type | Primary script | Expected output |
|---|---|---|---|
| Figure 3 | Monte Carlo convergence | `src/seg_model.py` or `src/single_mc.py` | `mc_trace.csv` / `trace__*.csv` |
| Figure 4 | Sensitivity to interaction parameter | `scripts/run_robustness_tests.py` | `robustness_results__summary.csv` |
| Figure 5 | Sensitivity to elastic factor `A` | `scripts/sweep_A_run.py` | `A_sweep_outputs/manifest.csv` and associated CSV files |
| Table 1 | Elastic contribution summary | `scripts/sweep_A_run.py` | processed `A_sweep_outputs/*.csv` |
| Figure 6 | Surface energy versus temperature | `src/seg_model.py` | `results_long.csv` or `results_wide.csv` |
| Figure 7 | Adsorption versus temperature | `src/seg_model.py` | `results_long.csv` or `results_wide.csv` |

Figure-generation scripts may be added separately if the final plotting workflow is stored independently from the numerical calculation scripts.

---

## Installation

Create and activate a clean Python environment. Example using `conda`:

```bash
conda create -n fe-surface-segregation python=3.11
conda activate fe-surface-segregation
pip install -r requirements.txt
```

Minimum expected dependencies:

```text
numpy
pandas
scipy
matplotlib
tqdm
pycalphad
```

A typical `requirements.txt` is:

```text
numpy
pandas
scipy
matplotlib
tqdm
pycalphad
```

For exact reproducibility, replace the unpinned requirements with the package versions used in the final analysis environment.

---

## Path configuration

The original scripts used local absolute paths for the thermodynamic database and the surface-energy JSON file. For public release, these paths should be replaced by repository-relative paths.

Recommended implementation inside `src/seg_model.py` and `src/single_mc.py`:

```python
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

TDB_PATH = PROJECT_ROOT / "data" / "thermodynamic" / "steel_database_fix.tdb"
SURFACES_JSON_PATH = PROJECT_ROOT / "data" / "surface_energy" / "surfaces.json"
```

If the scripts are executed from a different directory layout, adjust `PROJECT_ROOT` accordingly.

---

## Running the calculations

### 1. Main segregation sweep

Run the main multilayer segregation sweep:

```bash
python src/seg_model.py
```

Expected outputs, depending on script configuration:

```text
results_wide.csv
results_long.csv
mc_trace.csv
```

For repository consistency, move or write these files into:

```text
outputs/main_sweep/
```

### 2. Robustness analyses

Run seed, layer-number and interaction-parameter robustness tests:

```bash
python scripts/run_robustness_tests.py
```

Expected outputs:

```text
outputs/robustness_outputs/robustness_results__long.csv
outputs/robustness_outputs/robustness_results__summary.csv
```

### 3. Elastic-factor sweep

Run the sweep over the elastic scaling factor `A`:

```bash
python scripts/sweep_A_run.py
```

Expected outputs:

```text
outputs/A_sweep_outputs/manifest.csv
outputs/A_sweep_outputs/EL1__A*__wide.csv
outputs/A_sweep_outputs/EL1__A*__long.csv
outputs/A_sweep_outputs/EL0__A*__wide.csv
outputs/A_sweep_outputs/EL0__A*__long.csv
```

---

## Methodological assumptions

The calculations are thermodynamic and equilibrium-based. They do not simulate time-dependent cavity filling.

Main assumptions:
- substitutional binary Fe–X alloys;
- BCC_A2 matrix phase;
- multilayer surface region with finite number of atomic layers;
- regular-solution-like surface thermodynamic formulation;
- Monte Carlo minimisation used as an optimisation procedure;
- planar surface approximation for segregation calculations;
- no explicit cavity curvature in the segregation model;
- no diffusion kinetics or nucleation-rate calculation;
- no externally applied stress;
- optional elastic term treated only as a sensitivity parameter.

The outputs should therefore be interpreted as equilibrium descriptors of surface-conditioned thermodynamic states, not as kinetic predictions.

---

## Reproducibility notes

1. Monte Carlo results depend on the random seed. Robustness against seed variation is evaluated using `run_robustness_tests.py`.
2. The number of surface layers `L` affects adsorption estimates and should be reported when comparing calculations.
3. The interaction parameter scaling `omega_scale` is a sensitivity parameter and should not be interpreted as an independently fitted thermodynamic model.
4. The elastic scaling factor `A` is used to test sensitivity to a screened elastic contribution. It is not a validated kinetic parameter.
5. CALPHAD results depend on the `.tdb` database version. The exact file used should be archived or described with sufficient version information.
6. Surface-energy values depend on the specific `surfaces.json` dataset and the selected key, especially `weighted_surface_energy`.

---

## Recommended `.gitignore`

```gitignore
# Python
__pycache__/
*.pyc
*.pyo
*.pyd
.ipynb_checkpoints/

# Environments
.venv/
venv/
env/
*.egg-info/

# OS/editor files
.DS_Store
Thumbs.db
.vscode/
.idea/

# Temporary files
*.tmp
*.bak

# Large or non-essential Monte Carlo traces
outputs/**/mc_trace*.csv
outputs/**/trace*.csv

# Keep selected representative traces by explicitly unignoring them if needed
# !outputs/robustness_outputs/trace__representative_case.csv
```

Do not ignore CSV files that are required to reproduce manuscript figures and tables.

---

## Citation

If this repository is used, cite the associated manuscript once published.

For the pure-element surface-energy dataset, cite:

```text
Tran, R., Xu, Z., Radhakrishnan, B., Winston, D., Sun, W., Persson, K. A. & Ong, S. P.
Surface energies of elemental crystals.
Scientific Data 3, 160080 (2016).
https://doi.org/10.1038/sdata.2016.80
```

For the MatCalc thermodynamic database, cite the MatCalc open database source and any specific database documentation associated with the downloaded `.tdb` file.

---

## Licence

Specify the licence for the code in `LICENSE`.

Recommended options:
- MIT Licence for code, if permissive reuse is intended;
- CC BY 4.0 for documentation and derived non-code text, if desired.

Input datasets may have separate licence terms. In particular, verify redistribution permissions for the MatCalc `.tdb` file before publishing it in a public repository.
