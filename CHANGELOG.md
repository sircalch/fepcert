# Changelog

## 1.2.0 (2026-10-01)

Validated on real GROMACS data (alchemtest: benzene and ethanol hydration, the ligand and complex
legs of an absolute binding free energy in T4 lysozyme, a water particle; full, truncated and
sparse-window variants, 28 analyses) against alchemlyb 2.5 and pymbar 4, and on the 92
perturbation maps of the public FEP+ benchmark (1404 cycles) against networkx (`validation/`).
TI, BAR, its uncertainty, the forward/backward convergence estimates and the new overlap agree
with the references to rounding error.

### Fixed
- **Phase-space overlap was not an overlap.** It was the overlap of histograms of dH/dlambda in
  neighbouring windows, compared with the 0.03 threshold that Klimovich et al. defined for the
  MBAR overlap matrix. On the benchmark it was uncorrelated with the MBAR overlap and passed a set
  of windows whose MBAR overlap was 0.006, where TI was off by 2 kJ/mol. It is now the two-state
  MBAR overlap of each neighbouring pair, computed from the forward and reverse work (agrees with
  pymbar to 4e-12); with dH/dlambda only, the check is NOT_APPLICABLE.
- **Cycle enumeration missed cycles.** Four-node cycles were identified by their node set, so the
  distinct cycles on the same four ligands collapsed into one: 104 of 1404 cycles (18 of 92 maps)
  were missed and the largest closure error was underestimated (1.25 instead of 1.74 kcal/mol for
  TYK2). Cycles are now identified by their edge set; repeated edges are combined.
- **Thresholds ignored the unit.** Convergence, TI and BAR thresholds are in kcal/mol but were
  applied to values in kJ/mol unchanged (API, or the CLI with `--unit kJ/mol`), giving false
  FAIL/WARNING verdicts. They are now converted.
- Methods paragraph: overlap was "certified ... exceeding the threshold" and cycle closure
  "confirmed" whatever the values; it now states the results, and the version number is read from
  the package.

### Added
- **Multi-component lambda paths** (bonded, Coulomb and van der Waals changing in one set of
  windows, as in GROMACS absolute binding free energies): `calculate_ti_free_energy_multi`; the
  directory parser returns `lambda_matrix` and `gradients_by_component` instead of raising.
- **Neighbour work from GROMACS files**: the parser reads the Delta H columns and returns
  `forward_work_list` / `reverse_work_list`, so BAR and the overlap are computed from dhdl.xvg
  directly; the CLI reports BAR next to TI.
- **TI/BAR agreement check** (`check_estimator_agreement`): trapezoidal TI carries a quadrature
  error when the windows are too sparse, which neither its error bar nor the overlap shows; a
  difference from BAR on the same windows above three combined standard errors and 0.25 (1.0)
  kcal/mol gives WARNING (FAIL).
- **Cycle closure against the edge uncertainties**: `expected_sigma` and `z_score` per cycle and
  `frac_cycles_z_above_2` per network.

## 1.1.0 (2026-09-26)

### Fixed
- GROMACS parser: lambda and temperature read from the headers, dH/dlambda columns selected by
  legend, compressed files, recursive search; windows were ordered lexicographically by file name
  (VDW leg of benzene: -11.66 instead of -7.62 kJ/mol).
- GROMACS kJ/mol values converted to the reporting unit (they were labelled kcal/mol).
- TI error propagation of interior windows; per-window standard error corrected for autocorrelation.
- BAR asymptotic variance (Shirts et al. 2003; matches pymbar 4.0.3); time-convergence hysteresis
  compares disjoint halves (it was identically zero).

## 1.0.0 (2026-08-31)

Initial release.
