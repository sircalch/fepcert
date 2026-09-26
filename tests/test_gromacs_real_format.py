"""
Regression tests for real GROMACS dhdl.xvg headers, lambda ordering, units and TI error
propagation (bugs found by validation against alchemlyb on the alchemtest benzene data).
"""
import bz2
import os

import numpy as np
import pytest

from fepcert.core.ti_bar import calculate_ti_free_energy, statistical_inefficiency
from fepcert.parsers.gromacs_dhdl import parse_gromacs_dhdl, parse_gromacs_dhdl_directory
from fepcert.cli import convert_energy_unit

SINGLE = r"""# GROMACS dhdl
@    title "dH/d\xl\f{} and \xD\f{}H"
@    xaxis  label "Time (ps)"
@    yaxis  label "dH/d\xl\f{} and \xD\f{}H (kJ/mol [\xl\f{}]\S-1\N)"
@ subtitle "T = 300 (K) \xl\f{} state {state}: fep-lambda = {lam:.4f}"
@ s0 legend "dH/d\xl\f{} fep-lambda = {lam:.4f}"
@ s1 legend "\xD\f{}H \xl\f{} to 0.0000"
@ s2 legend "pV (kJ/mol)"
"""

MULTI = r"""@    yaxis  label "dH/d\xl\f{} and \xD\f{}H (kJ/mol [\xl\f{}]\S-1\N)"
@ subtitle "T = 298 (K) \xl\f{} state 7: (coul-lambda, vdw-lambda) = (1.0000, 0.3000)"
@ s0 legend "Total Energy (kJ/mol)"
@ s1 legend "dH/d\xl\f{} coul-lambda = 1.0000"
@ s2 legend "dH/d\xl\f{} vdw-lambda = 0.3000"
0.0 -1000.0 1.5 -4.0
2.0 -1001.0 1.7 -3.5
"""


def _write_window(path, state, lam, grad, n=50, compress=False):
    text = SINGLE.replace("{state}", str(state)).replace("{lam:.4f}", f"{lam:.4f}")
    rows = "".join(f"{2.0 * i:.1f} {grad:.4f} 0.1 0.77\n" for i in range(n))
    data = (text + rows).encode()
    if compress:
        with bz2.open(path + ".bz2", "wb") as fh:
            fh.write(data)
    else:
        with open(path, "wb") as fh:
            fh.write(data)


def test_real_single_component_header(tmp_path):
    p = tmp_path / "dhdl.xvg"
    _write_window(str(p), 3, 0.2, 5.0)
    d = parse_gromacs_dhdl(str(p))
    assert d["lambda"] == pytest.approx(0.2)
    assert d["temperature_k"] == pytest.approx(300.0)
    assert d["unit"] == "kJ/mol"
    assert np.allclose(d["gradients"], 5.0)   # dH/dl column, not a Delta H or pV column


def test_multi_component_header_selects_columns_by_legend(tmp_path):
    p = tmp_path / "dhdl.xvg"
    p.write_text(MULTI)
    d = parse_gromacs_dhdl(str(p))
    assert d["temperature_k"] == pytest.approx(298.0)
    assert set(d["components"]) == {"coul", "vdw"}
    assert d["components"]["vdw"]["lambda"] == pytest.approx(0.3)
    assert np.allclose(d["components"]["vdw"]["gradients"], [-4.0, -3.5])  # not the energy column
    assert np.allclose(d["components"]["coul"]["gradients"], [1.5, 1.7])


def test_windows_are_ordered_by_header_lambda_not_file_name(tmp_path):
    # Lexicographic order of these names is 0, 1000, 250, 500, 750.
    lams = {"dhdl.0.xvg": 0.0, "dhdl.250.xvg": 0.25, "dhdl.500.xvg": 0.5,
            "dhdl.750.xvg": 0.75, "dhdl.1000.xvg": 1.0}
    for i, (name, lam) in enumerate(lams.items()):
        _write_window(str(tmp_path / name), i, lam, grad=10.0 * lam)   # dH/dl = 10 lambda
    d = parse_gromacs_dhdl_directory(str(tmp_path), "dhdl*.xvg")
    assert d["lambda_values"] == pytest.approx([0.0, 0.25, 0.5, 0.75, 1.0])
    r = calculate_ti_free_energy(d["lambda_values"], d["gradients_list"], unit="kJ/mol")
    assert r.delta_g == pytest.approx(5.0)   # integral of 10 lambda from 0 to 1


def test_compressed_windows_in_subdirectories(tmp_path):
    for i, lam in enumerate([0.0, 0.5, 1.0]):
        sub = tmp_path / f"{int(lam * 1000):04d}"
        sub.mkdir()
        _write_window(str(sub / "dhdl.xvg"), i, lam, grad=2.0, compress=True)
    d = parse_gromacs_dhdl_directory(str(tmp_path), "dhdl.xvg*", recursive=True)
    assert d["n_windows"] == 3
    assert d["lambda_values"] == pytest.approx([0.0, 0.5, 1.0])


def test_missing_lambda_raises_instead_of_guessing(tmp_path):
    for name in ("dhdl.a.xvg", "dhdl.b.xvg"):
        (tmp_path / name).write_text("0.0 1.0\n2.0 1.1\n")
    with pytest.raises(ValueError, match="lambda"):
        parse_gromacs_dhdl_directory(str(tmp_path), "dhdl*.xvg")


def test_kj_to_kcal_conversion():
    out = convert_energy_unit([np.array([4.184, 8.368])], "kJ/mol", "kcal/mol")
    assert np.allclose(out[0], [1.0, 2.0])
    assert convert_energy_unit([np.array([1.0])], "kJ/mol", "kJ/mol")[0][0] == 1.0


def test_trapezoid_error_counts_interior_windows_in_both_intervals():
    # Three equally spaced windows with known standard errors: weights are 1/4, 1/2, 1/4.
    rng = np.random.default_rng(0)
    lams = [0.0, 0.5, 1.0]
    grads = [rng.normal(0.0, 1.0, 4000) for _ in lams]
    r = calculate_ti_free_energy(lams, grads, unit="kJ/mol")
    sig = np.array(r.gradient_errors)
    expected = np.sqrt((0.25 * sig[0]) ** 2 + (0.5 * sig[1]) ** 2 + (0.25 * sig[2]) ** 2)
    assert r.delta_g_error == pytest.approx(expected)


def test_ti_standard_error_is_corrected_for_autocorrelation():
    rng = np.random.default_rng(1)
    phi, n = 0.9, 20000
    x = np.empty(n)
    x[0] = rng.normal()
    e = rng.normal(0, np.sqrt(1 - phi * phi), n)
    for t in range(1, n):
        x[t] = phi * x[t - 1] + e[t]
    assert statistical_inefficiency(x) == pytest.approx(19.0, rel=0.15)
    r = calculate_ti_free_energy([0.0, 1.0], [x, x], unit="kJ/mol")
    naive = x.std(ddof=1) / np.sqrt(n)
    assert r.gradient_errors[0] == pytest.approx(naive * np.sqrt(19.0), rel=0.1)
