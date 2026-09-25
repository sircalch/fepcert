"""
Regression tests against reference implementations and known-answer cases.
"""
import numpy as np
import pytest

from fepcert.core.ti_bar import calculate_bar_free_energy
from fepcert.core.convergence import evaluate_time_convergence

KB = 0.001987204259
T = 298.15
BETA = 1.0 / (KB * T)


def test_bar_matches_pymbar_value_and_uncertainty():
    # Gaussian work, Crooks-consistent, true dF = 2.0 kT, unequal sample sizes.
    rng = np.random.default_rng(2024)
    wf = rng.normal(2.0 + 2.0, 2.0, 800)
    wr = rng.normal(-2.0 + 2.0, 2.0, 600)
    res = calculate_bar_free_energy([0.0, 1.0], [wf / BETA], [wr / BETA], temperature_k=T)
    # Reference: pymbar 4.0.3, pymbar.other_estimators.bar(wf, wr)
    assert res.delta_g * BETA == pytest.approx(1.983380106152966, abs=1e-6)
    assert res.delta_g_error * BETA == pytest.approx(0.05947453156592188, rel=1e-4)


def test_hysteresis_detects_drift_between_halves():
    # Gradients drift by +1 kcal/mol in the second half of every window -> dG shifts by ~1.
    rng = np.random.default_rng(0)
    lambdas = np.linspace(0.0, 1.0, 6)
    series = []
    for _ in lambdas:
        g = rng.normal(0.0, 0.1, 1000)
        g[500:] += 1.0
        series.append(g)
    res = evaluate_time_convergence(lambdas.tolist(), series)
    assert res.final_hysteresis == pytest.approx(1.0, abs=0.05)
    assert res.status == "FAIL"


def test_hysteresis_small_for_stationary_data():
    rng = np.random.default_rng(1)
    lambdas = np.linspace(0.0, 1.0, 6)
    series = [rng.normal(5.0, 0.5, 2000) for _ in lambdas]
    res = evaluate_time_convergence(lambdas.tolist(), series)
    assert 0.0 < res.final_hysteresis < 0.1
