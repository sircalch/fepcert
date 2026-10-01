# FEPCert

[![CI](https://github.com/sircalch/fepcert/actions/workflows/test.yml/badge.svg)](https://github.com/sircalch/fepcert/actions)
[![PyPI version](https://img.shields.io/pypi/v/fepcert.svg?color=blue)](https://pypi.org/project/fepcert/)
[![Python versions](https://img.shields.io/pypi/pyversions/fepcert.svg)](https://pypi.org/project/fepcert/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22217579.svg)](https://doi.org/10.5281/zenodo.22217579)

> **Quality checks for alchemical free-energy calculations: TI and BAR, phase-space overlap, time convergence and cycle closure.**

---

## Overview

**FEPCert** is an open-source Python package that checks alchemical free-energy calculations (TI and BAR, GROMACS dhdl.xvg or tables) and perturbation networks, and reports PASS / WARNING / FAIL / not applicable for each check.

- **TI and BAR**: trapezoidal TI with autocorrelation-corrected errors, including paths in which several lambda components change (bonded, Coulomb, van der Waals); BAR with the asymptotic variance of Shirts et al. (2003). Both match alchemlyb.
- **Phase-space overlap**: two-state MBAR overlap of every pair of neighbouring states from the forward and reverse work (read from the Delta H columns of GROMACS files), with the 0.03 threshold of Klimovich et al. (2015). With dH/dlambda only, the check is reported as not applicable.
- **Time convergence**: TI from the first and second halves of the data, and forward/backward estimates.
- **Cycle closure**: every simple cycle up to four edges of a perturbation network, its closure error and, when edge uncertainties are given, its z-score against the combined uncertainty.

> **Version 1.2.0** was validated on real GROMACS data (alchemtest, 28 analyses) against alchemlyb and pymbar, and on the 92 maps of the public FEP+ benchmark (1404 cycles). It corrects an overlap metric that did not measure overlap, cycles missed by the enumeration and thresholds that ignored the energy unit; see the CHANGELOG.

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

