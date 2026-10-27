import numpy as np
import pytest

from mc_engine import PricingResult


def test_from_samples_hand_computed():
    # mean = 3, sample variance = 2.5, so SE = sqrt(2.5 / 5) = sqrt(0.5)
    r = PricingResult.from_samples([1, 2, 3, 4, 5], runtime_sec=0.1)
    se = np.sqrt(0.5)
    assert r.price == pytest.approx(3.0)
    assert r.std_error == pytest.approx(se)
    assert r.ci_95 == pytest.approx((3 - 1.959964 * se, 3 + 1.959964 * se))
    assert r.n_paths == 5
    assert r.runtime_sec == 0.1


def test_estimate_within_3_standard_errors_of_true_mean():
    # Same check every pricer will get: |MC - truth| < 3 SE.
    rng = np.random.default_rng(42)
    r = PricingResult.from_samples(rng.normal(loc=1.0, scale=2.0, size=100_000))
    assert abs(r.price - 1.0) < 3 * r.std_error


def test_rejects_single_sample():
    with pytest.raises(ValueError):
        PricingResult.from_samples([1.0])
