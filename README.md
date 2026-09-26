# FEPCert

[![CI](https://github.com/sircalch/fepcert/actions/workflows/test.yml/badge.svg)](https://github.com/sircalch/fepcert/actions)
[![PyPI version](https://img.shields.io/pypi/v/fepcert.svg?color=blue)](https://pypi.org/project/fepcert/)
[![Python versions](https://img.shields.io/pypi/pyversions/fepcert.svg)](https://pypi.org/project/fepcert/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234575.svg)](https://doi.org/10.5281/zenodo.1234575)

> **Automated Quality-Control, Phase Space Overlap, and Thermodynamic Cycle Closure Certification for Alchemical Free Energy Simulations (FEP, TI, BAR; overlap diagnostics for BAR/MBAR workflows).**

---

## Overview

**FEPCert** is an open-source scientific software package designed to systematically audit, validate, and certify the statistical convergence, phase space connectivity, and internal consistency of alchemical free energy calculations (Free Energy Perturbation, Thermodynamic Integration, Bennett Acceptance Ratio, and Multistate BAR) from simulation engines such as **GROMACS**, **AMBER**, **NAMD**, and **OpenMM**.

In computational drug discovery and biophysics, assessing whether an alchemical relative binding free energy ($\Delta\Delta G_{\text{bind}}$) calculation has converged is essential to prevent false predictions:

- 📊 **Phase Space Overlap Matrix ($\Pi_{ij}$)**:
  - Computes pairwise probability density overlap between adjacent $\lambda$-states ($\Pi_{k, k+1}$).
  - Enforces the Klimovich criterion ($\min \Pi_{k, k+1} \ge 0.03$, 3%) to guarantee unbiased BAR/MBAR reweighting.
  - Automatically flags phase space bottlenecks where intermediate $\lambda$-windows are required.
- 📐 **Thermodynamic Integration (TI) & BAR**:
  - Exact numerical quadrature of $\langle \partial H/\partial\lambda \rangle_\lambda$ with standard error propagation.
  - Self-consistent Bennett Acceptance Ratio (BAR) solution with asymptotic variance error bounds.
- ⏳ **Time-Series Convergence & Hysteresis**:
  - Compares forward cumulative $\Delta G(t)$ and reverse cumulative $\Delta G(t)$ trajectories.
  - Quantifies dissipated work $W_{\text{diss}} = \frac{1}{2}|\Delta G_{\text{fwd}} + \Delta G_{\text{rev}}|$ and flags simulations trapped in metastable conformations.
- 🔄 **Thermodynamic Cycle Closure ($\oint \Delta\Delta G \approx 0$)**:
  - Graph-based cycle detection (e.g. triangles $L_1 \rightarrow L_2 \rightarrow L_3 \rightarrow L_1$ and quadrilaterals).
  - Computes cycle closure error $\sum_{(i,j) \in C} \Delta\Delta G_{ij}$ and overall cycle $\text{RMSE} \le 0.50$ kcal/mol.
- 📑 **Publication Deliverables**:
  - Interactive self-contained `report.html` dashboard.
  - Publication vector plots (Overlap matrix heatmap $\Pi_{ij}$, $\langle \partial H/\partial\lambda \rangle$ integration curve, forward/reverse convergence lines) in SVG, PDF, PNG (300 DPI).
  - Ready-to-compile LaTeX summary tables (`.tex`).
  - Draft **Methods** text snippet and BibTeX citation (`citation.bib`).

```
   GROMACS dhdl*.xvg / Alchemical CSV Table
                      │
                      ▼
  ┌───────────────────────────────────────────────────────────┐
  │                         FEPCert                           │
  │  ├── Phase Space Overlap Matrix (\Pi_ij >= 3% Klimovich)  │
  │  ├── Thermodynamic Integration & BAR Estimator            │
  │  ├── Forward-Reverse Time Convergence & Dissipated Work   │
  │  └── Perturbation Network Cycle Closure Audit (\oint \Delta\Delta G)│
  └───────────────────────────────────────────────────────────┘
                      │
                      ▼
  ┌───────────────────────────────────────────────────────────┐
  │                   Publication Deliverables                │
  │  ├── report.html (Interactive Dashboard & Badges)         │
  │  ├── fepcert_overlap_matrix.pdf/svg/png                   │
  │  ├── fepcert_ti_gradient_curve.pdf/svg/png                │
  │  ├── fepcert_time_convergence.pdf/svg/png                 │
  │  ├── fepcert_summary_table.tex / .csv                     │
  │  ├── methods_snippet.txt (Ready for Manuscript)           │
  │  └── citation.bib (BibTeX Reference)                      │
  └───────────────────────────────────────────────────────────┘
```

---

## Installation

### From PyPI
> **Note:** PyPI release pending. Until then, install from the tagged GitHub release:

```bash
pip install "git+https://github.com/sircalch/fepcert@v1.0.0"
```

### From Source
```bash
git clone https://github.com/sircalch/fepcert.git
cd fepcert
pip install -e .[dev]
```

---

## Quickstart (CLI)

### 1. Run Benchmark Demo (Instant Alchemical Simulation & 4-Ligand Cycle Audit)
```bash
fepcert demo -o my_fep_audit/
```
Open `my_fep_audit/report.html` in any browser to inspect the interactive report!

### 2. Assess GROMACS dhdl Files
```bash
fepcert assess -d gromacs_fep_run/ --pattern "dhdl*.xvg*" -o fep_report/
# one sub-directory per window (e.g. 0000/dhdl.xvg.bz2):
fepcert assess -d gromacs_fep_run/ --pattern "dhdl.xvg*" --recursive -o fep_report/
```

Lambda values, temperature and the energy unit are read from each file header (`fep-lambda = 0.2000`, `(coul-lambda, vdw-lambda) = (...)`), and dH/dλ columns are identified by their legends. Windows are ordered by their header λ, never by file name; if λ cannot be read, FEPCert stops with an error instead of guessing. GROMACS energies (kJ/mol) are converted to the reporting unit (`--unit`, default kcal/mol). Analyse each leg (e.g. Coulomb and van der Waals) in its own directory.

### 3. Audit Thermodynamic Cycle Closure in a Chemical Network
```bash
fepcert cycle -i perturbation_network.csv -o cycle_report/
```

---

## Validation

`validation/validate_alchemtest.py` compares FEPCert with alchemlyb 2.5 on the real GROMACS benzene
decoupling data set of `alchemtest` (2 legs, 21 windows of 40 ns). The raw `.xvg.bz2` files are read directly.

| Leg | TI FEPCert (kJ/mol) | TI alchemlyb | BAR FEPCert | BAR alchemlyb |
|---|---|---|---|---|
| Coulomb | 7.705 ± 0.054 | 7.705 ± 0.054 | 7.594 ± 0.041 | 7.594 ± 0.041 |
| VDW | −7.622 ± 0.123 | −7.622 ± 0.121 | −7.565 ± 0.086 | −7.565 ± 0.086 |

The TI standard error of each window is corrected for autocorrelation (s·√(g/n)). This is why the VDW TI
error is slightly larger than alchemlyb's default.

Version 1.1.0 corrects the following defects found by this validation:

- λ was not read from real GROMACS headers and was then assigned by lexicographic file order. VDW gave
  −11.66 instead of −7.62 kJ/mol.
- Values in kJ/mol were labelled as kcal/mol without conversion.
- The trapezoid error ignored that interior windows enter two intervals.
- The TI standard error ignored autocorrelation.
- The BAR variance and the hysteresis metric were corrected earlier; see the git history.

## Python API Usage

```python
from fepcert import assess_fep_quality
from fepcert.parsers import parse_gromacs_dhdl_directory, parse_perturbation_network_csv
from fepcert.reporters import generate_fepcert_figures, generate_fepcert_manuscript_assets, generate_fepcert_html_report

# 1. Parse GROMACS directory & perturbation network
from fepcert.cli import convert_energy_unit

data = parse_gromacs_dhdl_directory("fep_simulations/", file_pattern="dhdl*.xvg*")
grads = convert_energy_unit(data["gradients_list"], data["unit"], "kcal/mol")  # GROMACS writes kJ/mol
network = parse_perturbation_network_csv("network.csv")

# 2. Assess free energy convergence
report = assess_fep_quality(
    metadata={"transformation": "Lig1 -> Lig2", "engine": "GROMACS 2024"},
    lambda_values=data["lambda_values"],
    gradients_list=grads,
    network_edges=network,
    unit="kcal/mol"
)

print(f"Overall Certification: {report.overall_status}")
print(f"Free Energy: {report.free_energy.delta_g:.2f} +/- {report.free_energy.delta_g_error:.2f} kcal/mol")
print(f"Min Adjacent Overlap: {report.overlap_result.min_adjacent_overlap*100:.1f}%")

# 3. Export all publication assets
generate_fepcert_figures(report, "output_dir/")
generate_fepcert_manuscript_assets(report, "output_dir/")
generate_fepcert_html_report(report, "output_dir/report.html")
```

---

## Citation

If you use FEPCert in your research, please cite:

```bibtex
@software{monreal2026fepcert,
  author = {Monreal-Hern{\'a}ndez, Andre},
  title = {{FEPCert: An Open-Source Toolkit for Quality-Control, Phase Space Overlap, and Thermodynamic Cycle Closure Certification of Alchemical Free Energy Simulations}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/sircalch/fepcert}
}
```

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

