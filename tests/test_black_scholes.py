import numpy as np
import pytest

from mc_engine import PricingResult
from mc_engine.analytics import black_scholes as bs

# A random grid of inputs (fixed seed) used for the identity checks below.
_rng = np.random.default_rng(0)
N = 200
GRID = dict(
    S=_rng.uniform(50, 150, N),
    K=_rng.uniform(50, 150, N),
    T=_rng.uniform(0.05, 3, N),
    r=_rng.uniform(-0.01, 0.10, N),
    sigma=_rng.uniform(0.05, 0.8, N),
)


# ---- known textbook values ------------------------------------------------

def test_standard_example():
    # S=K=100, T=1, r=5%, sigma=20%: the most-quoted Black-Scholes example.
    assert bs.price(100, 100, 1, 0.05, 0.2, "call") == pytest.approx(10.4506, abs=1e-4)
    assert bs.price(100, 100, 1, 0.05, 0.2, "put") == pytest.approx(5.5735, abs=1e-4)


def test_hull_bsm_example():
    # Hull, Options Futures & Other Derivatives, BSM chapter: S=42, K=40, T=0.5, r=10%, sigma=20%.
    assert bs.price(42, 40, 0.5, 0.10, 0.2, "call") == pytest.approx(4.76, abs=0.005)
    assert bs.price(42, 40, 0.5, 0.10, 0.2, "put") == pytest.approx(0.81, abs=0.005)


def test_hull_greeks_example():
    # Hull's Greeks chapter: S=49, K=50, T=20 weeks, r=5%, sigma=20%, call.
    args = (49, 50, 20 / 52, 0.05, 0.2)
    assert bs.price(*args) == pytest.approx(2.40, abs=0.005)
    assert bs.delta(*args) == pytest.approx(0.522, abs=5e-4)
    assert bs.gamma(*args) == pytest.approx(0.066, abs=5e-4)
    assert bs.vega(*args) == pytest.approx(12.1, abs=0.05)
    assert bs.theta(*args) == pytest.approx(-4.31, abs=0.005)
    assert bs.rho(*args) == pytest.approx(8.91, abs=0.005)


# ---- put-call parity: C - P = S - K e^{-rT}, and its derivatives -----------

def test_put_call_parity_price():
    S, K, T, r = GRID["S"], GRID["K"], GRID["T"], GRID["r"]
    diff = bs.price(**GRID, kind="call") - bs.price(**GRID, kind="put")
    np.testing.assert_allclose(diff, S - K * np.exp(-r * T), atol=1e-10)


def test_put_call_parity_greeks():
    # Differentiate both sides of parity: each Greek difference is known exactly.
    K, T, r = GRID["K"], GRID["T"], GRID["r"]
    disc_K = K * np.exp(-r * T)
    np.testing.assert_allclose(bs.delta(**GRID) - bs.delta(**GRID, kind="put"), 1.0, atol=1e-12)
    np.testing.assert_allclose(bs.theta(**GRID) - bs.theta(**GRID, kind="put"), -r * disc_K, atol=1e-10)
    np.testing.assert_allclose(bs.rho(**GRID) - bs.rho(**GRID, kind="put"), T * disc_K, atol=1e-10)


@pytest.mark.parametrize("kind", ["call", "put"])
def test_greeks_satisfy_black_scholes_pde(kind):
    # Theta + 1/2 sigma^2 S^2 Gamma + r S Delta = r V ties four Greeks and the price together.
    S, r, sig = GRID["S"], GRID["r"], GRID["sigma"]
    lhs = (bs.theta(**GRID, kind=kind) + 0.5 * sig**2 * S**2 * bs.gamma(**GRID)
           + r * S * bs.delta(**GRID, kind=kind))
    np.testing.assert_allclose(lhs, r * bs.price(**GRID, kind=kind), atol=1e-9)


# ---- every Greek equals a finite difference of the price ------------------

@pytest.mark.parametrize("kind", ["call", "put"])
def test_greeks_match_finite_differences(kind):
    S, K, T, r, sig = GRID["S"], GRID["K"], GRID["T"], GRID["r"], GRID["sigma"]
    h = 1e-4

    def V(S=S, T=T, r=r, sig=sig):
        return bs.price(S, K, T, r, sig, kind)

    fd = {
        "delta": (V(S=S + h) - V(S=S - h)) / (2 * h),
        "gamma": (V(S=S + h) - 2 * V() + V(S=S - h)) / h**2,
        "vega": (V(sig=sig + h) - V(sig=sig - h)) / (2 * h),
        "theta": -(V(T=T + h) - V(T=T - h)) / (2 * h),  # time passing = T shrinking
        "rho": (V(r=r + h) - V(r=r - h)) / (2 * h),
    }
    # Central differences carry O(h^2) truncation error (about 4e-5 relative for vega at
    # sigma=5%), so compare relatively. A wrong formula would be off by O(1), not 1e-4.
    tol = dict(rtol=1e-4, atol=1e-6)
    np.testing.assert_allclose(bs.delta(S, K, T, r, sig, kind), fd["delta"], **tol)
    np.testing.assert_allclose(bs.gamma(S, K, T, r, sig), fd["gamma"], rtol=1e-4, atol=1e-4)
    np.testing.assert_allclose(bs.vega(S, K, T, r, sig), fd["vega"], **tol)
    np.testing.assert_allclose(bs.theta(S, K, T, r, sig, kind), fd["theta"], **tol)
    np.testing.assert_allclose(bs.rho(S, K, T, r, sig, kind), fd["rho"], **tol)


# ---- the formula really is the discounted risk-neutral expectation --------

@pytest.mark.parametrize("K", [80, 100, 120])
@pytest.mark.parametrize("kind", ["call", "put"])
def test_monte_carlo_within_3_standard_errors(K, kind):
    S0, T, r, sigma = 100, 1.0, 0.05, 0.2
    Z = np.random.default_rng(42).standard_normal(1_000_000)
    ST = S0 * np.exp((r - 0.5 * sigma**2) * T + sigma * np.sqrt(T) * Z)
    payoff = np.maximum(ST - K, 0) if kind == "call" else np.maximum(K - ST, 0)
    mc = PricingResult.from_samples(np.exp(-r * T) * payoff)
    assert abs(mc.price - bs.price(S0, K, T, r, sigma, kind)) < 3 * mc.std_error


# ---- input handling -------------------------------------------------------

@pytest.mark.parametrize("bad", [dict(S=0), dict(K=-1), dict(T=0), dict(sigma=0)])
def test_rejects_non_positive_inputs(bad):
    args = dict(S=100, K=100, T=1, r=0.05, sigma=0.2) | bad
    with pytest.raises(ValueError):
        bs.price(**args)


def test_rejects_unknown_kind():
    with pytest.raises(ValueError):
        bs.price(100, 100, 1, 0.05, 0.2, kind="straddle")
