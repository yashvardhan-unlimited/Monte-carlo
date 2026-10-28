import numpy as np
import pytest

from mc_engine.models.gbm import GBM

MODEL = GBM(S0=100, r=0.05, sigma=0.2)
N_PATHS, N_STEPS, T = 200_000, 12, 1.0


@pytest.fixture(scope="module")
def paths():
    return MODEL.simulate(N_PATHS, N_STEPS, T, np.random.default_rng(1))


def test_shape_start_and_positivity(paths):
    assert paths.shape == (N_PATHS, N_STEPS + 1)
    assert np.all(paths[:, 0] == 100)
    assert np.all(paths > 0)


def test_same_seed_same_paths():
    a = MODEL.simulate(10, 5, T, np.random.default_rng(7))
    b = MODEL.simulate(10, 5, T, np.random.default_rng(7))
    c = MODEL.simulate(10, 5, T, np.random.default_rng(8))
    np.testing.assert_array_equal(a, b)
    assert not np.array_equal(a, c)


def test_zero_volatility_grows_at_risk_free_rate():
    p = GBM(S0=100, r=0.05, sigma=0.0).simulate(3, 4, T, np.random.default_rng(0))
    t = np.linspace(0, T, 5)
    np.testing.assert_allclose(p, np.broadcast_to(100 * np.exp(0.05 * t), (3, 5)))


def test_discounted_price_is_a_martingale(paths):
    # Under Q, E[e^{-rt} S_t] = S0 at every grid time t, not just at maturity.
    t = np.linspace(0, T, N_STEPS + 1)
    disc = np.exp(-MODEL.r * t) * paths
    mean = disc.mean(axis=0)
    se = disc.std(axis=0, ddof=1) / np.sqrt(N_PATHS)
    assert np.all(np.abs(mean[1:] - 100) < 3 * se[1:])


def test_log_price_mean_and_variance(paths):
    # ln(S_t/S0) ~ N((r - sigma^2/2) t, sigma^2 t): the -sigma^2/2 is what this checks.
    t = np.linspace(0, T, N_STEPS + 1)[1:]
    x = np.log(paths[:, 1:] / 100)
    mean, var = x.mean(axis=0), x.var(axis=0, ddof=1)
    true_mean, true_var = (0.05 - 0.02) * t, 0.04 * t
    assert np.all(np.abs(mean - true_mean) < 3 * np.sqrt(true_var / N_PATHS))
    # For normal data, the sample variance has relative standard error sqrt(2/(n-1)).
    assert np.all(np.abs(var / true_var - 1) < 3 * np.sqrt(2 / (N_PATHS - 1)))


def test_increments_are_uncorrelated(paths):
    # Independent increments: neighbouring log-returns should have ~0 correlation.
    inc = np.diff(np.log(paths), axis=1)
    corr = np.corrcoef(inc[:, 0], inc[:, 1])[0, 1]
    assert abs(corr) < 3 / np.sqrt(N_PATHS)


@pytest.mark.parametrize("bad", [dict(S0=0), dict(sigma=-0.1)])
def test_rejects_bad_parameters(bad):
    with pytest.raises(ValueError):
        GBM(**(dict(S0=100, r=0.05, sigma=0.2) | bad))
