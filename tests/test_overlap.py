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


def test_overlap_matrix_pipeline():
    rng = np.random.default_rng(42)
    # 5 lambda states with gradual shift
    data = [rng.normal(i * 1.0, 1.5, 300) for i in range(5)]
    res = calculate_overlap_matrix(data, min_pass_overlap=0.05)
    
    assert res.n_states == 5
    assert res.overlap_matrix.shape == (5, 5)
    assert len(res.adjacent_overlaps) == 4
    assert res.status == "PASS"
    assert res.min_adjacent_overlap > 0.05
