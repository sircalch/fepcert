"""
Tests for phase space overlap matrix calculation.
"""

import numpy as np
import pytest
from fepcert.core.overlap import calculate_overlap_matrix, _calculate_1d_distribution_overlap


def test_distribution_overlap_identical():
    s = np.random.normal(0, 1, 500)
    ov = _calculate_1d_distribution_overlap(s, s)
    assert np.isclose(ov, 1.0, atol=0.05)


def test_distribution_overlap_disjoint():
    s1 = np.random.normal(0, 0.5, 500)
    s2 = np.random.normal(10, 0.5, 500)
    ov = _calculate_1d_distribution_overlap(s1, s2)
    assert ov < 0.01


def test_gradients_alone_are_not_applicable():
    rng = np.random.default_rng(42)
    data = [rng.normal(i * 1.0, 1.5, 300) for i in range(5)]
    res = calculate_overlap_matrix(data)
    assert res.status == "NOT_APPLICABLE"


def test_two_state_overlap_limits():
    from fepcert.core.overlap import two_state_overlap
    zeros = np.zeros(400)
    # identical states: every sample is equally likely in both, O_01 = N1/(N0+N1) = 0.5
    assert two_state_overlap(zeros, zeros) == pytest.approx(0.5, abs=1e-12)
    # well separated states (large work in both directions): overlap vanishes
    assert two_state_overlap(np.full(400, 40.0), np.full(400, 40.0)) < 1e-10


def test_overlap_from_work_detects_a_bottleneck():
    rng = np.random.default_rng(0)
    # Gaussian work with variance s^2 and mean s^2/2 (linear response): overlap falls as s grows
    fw, rw = [], []
    for s in (0.5, 0.5, 6.0, 0.5):
        fw.append(rng.normal(s * s / 2, s, 2000))
        rw.append(rng.normal(s * s / 2, s, 2000))
    kt = 0.008314462618 * 300.0  # work above is in kT; convert to kJ/mol
    fw = [f * kt for f in fw]
    rw = [r * kt for r in rw]
    res = calculate_overlap_matrix(forward_work_list=fw, reverse_work_list=rw, temperature_k=300.0, unit="kJ/mol")
    assert res.n_states == 5
    assert res.adjacent_overlaps[2] < 0.03 < res.adjacent_overlaps[0]
    assert res.status == "FAIL"
