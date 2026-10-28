import numpy as np
import pytest

from mc_engine import PricingResult
from mc_engine.analytics import black_scholes as bs
from mc_engine.engine import price
from mc_engine.models.gbm import GBM
from mc_engine.payoffs.european import European


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


# ---- price(): GBM + European payoff vs Black-Scholes ------------------------

GBM_MODEL = GBM(S0=100, r=0.05, sigma=0.2)


@pytest.mark.parametrize("n_steps", [1, 12])
@pytest.mark.parametrize("kind", ["call", "put"])
@pytest.mark.parametrize("T", [0.25, 1.0, 2.0])
@pytest.mark.parametrize("K", [80, 90, 100, 110, 120])
def test_price_matches_black_scholes_within_3_se(K, T, kind, n_steps):
    # The exact GBM step has no discretisation error, so 1 step and 12 steps must both
    # reproduce the closed form; only Monte Carlo noise (the SE) separates them.
    res = price(GBM_MODEL, European(K, kind), T, n_paths=200_000, n_steps=n_steps)
    exact = bs.price(100, K, T, 0.05, 0.2, kind)
    assert abs(res.price - exact) < 3 * res.std_error


def test_price_result_fields_and_reproducibility():
    a = price(GBM_MODEL, European(100), 1.0, n_paths=1000, seed=3)
    b = price(GBM_MODEL, European(100), 1.0, n_paths=1000, seed=3)
    assert a.price == b.price
    assert a.n_paths == 1000 and a.runtime_sec > 0
    assert a.ci_95[0] < a.price < a.ci_95[1]


def test_standard_error_shrinks_like_one_over_sqrt_n():
    small = price(GBM_MODEL, European(100), 1.0, n_paths=10_000)
    large = price(GBM_MODEL, European(100), 1.0, n_paths=1_000_000)
    assert large.std_error / small.std_error == pytest.approx(0.1, rel=0.05)


def test_confidence_interval_covers_about_95_percent():
    # Over many independent runs, the 95% CI should contain the true price ~95% of the time.
    # With 200 runs, the hit count is Binomial(200, 0.95): sd = 1.5%, so 3 sd = [90.4%, 99.6%].
    exact = bs.price(100, 100, 1.0, 0.05, 0.2)
    hits = [
        r.ci_95[0] < exact < r.ci_95[1]
        for r in (price(GBM_MODEL, European(100), 1.0, n_paths=10_000, seed=s) for s in range(200))
    ]
    assert 0.904 < np.mean(hits) < 0.996
