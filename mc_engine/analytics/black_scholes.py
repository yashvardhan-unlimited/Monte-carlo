"""Black-Scholes closed-form prices and Greeks for European options (no dividends).

These are the benchmark every Monte Carlo pricer in this project is tested against.
Every function accepts scalars or NumPy arrays (broadcast together), so one call can
price a whole grid of strikes or spots at once.

Units and conventions (explained in docs/02_black_scholes.md):
    S, K   prices, same currency.
    T      time to maturity in years.
    r      continuously compounded risk-free rate, decimal (0.05 = 5%).
    sigma  volatility per sqrt(year), decimal (0.2 = 20%).
    vega, rho   change in price per 1.00 change in sigma or r. Divide by 100 for "per 1%".
    theta       change in price per year of calendar time passing. Divide by 365 for "per day".

Call and put share one formula through a sign phi = +1 (call) or -1 (put):
    price = phi * [S N(phi d1) - K e^{-rT} N(phi d2)]
Writing it once means a fix to one cannot silently miss the other.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import norm

Kind = Literal["call", "put"]
FloatOrArray = float | NDArray[np.float64]


def _phi(kind: Kind) -> int:
    """+1 for a call, -1 for a put."""
    if kind == "call":
        return 1
    if kind == "put":
        return -1
    raise ValueError(f"kind must be 'call' or 'put', got {kind!r}")


def _prepare(S: ArrayLike, K: ArrayLike, T: ArrayLike, r: ArrayLike, sigma: ArrayLike):
    """Convert inputs to float arrays, validate them, and compute d1 and d2.

    The formula divides by sigma*sqrt(T) and takes log(S/K), so zero or negative
    values would silently produce inf/NaN. Fail loudly instead.
    """
    S, K, T, r, sigma = (np.asarray(x, dtype=float) for x in (S, K, T, r, sigma))
    if np.any(S <= 0) or np.any(K <= 0) or np.any(T <= 0) or np.any(sigma <= 0):
        raise ValueError("S, K, T and sigma must all be strictly positive")
    vol_sqrt_t = sigma * np.sqrt(T)  # total volatility until maturity
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / vol_sqrt_t
    d2 = d1 - vol_sqrt_t
    return S, K, T, r, sigma, d1, d2


def price(
    S: ArrayLike, K: ArrayLike, T: ArrayLike, r: ArrayLike, sigma: ArrayLike, kind: Kind = "call"
) -> FloatOrArray:
    """Black-Scholes price of a European call or put.

    Call: S N(d1) - K e^{-rT} N(d2).   Put: K e^{-rT} N(-d2) - S N(-d1).
    """
    phi = _phi(kind)
    S, K, T, r, sigma, d1, d2 = _prepare(S, K, T, r, sigma)
    return phi * (S * norm.cdf(phi * d1) - K * np.exp(-r * T) * norm.cdf(phi * d2))


def delta(
    S: ArrayLike, K: ArrayLike, T: ArrayLike, r: ArrayLike, sigma: ArrayLike, kind: Kind = "call"
) -> FloatOrArray:
    """dV/dS: shares of stock needed to hedge one option. Call N(d1), put N(d1) - 1."""
    phi = _phi(kind)
    *_, d1, _d2 = _prepare(S, K, T, r, sigma)
    return phi * norm.cdf(phi * d1)


def gamma(S: ArrayLike, K: ArrayLike, T: ArrayLike, r: ArrayLike, sigma: ArrayLike) -> FloatOrArray:
    """d^2V/dS^2: how fast delta changes. Same for calls and puts (parity is linear in S)."""
    S, K, T, r, sigma, d1, _d2 = _prepare(S, K, T, r, sigma)
    return norm.pdf(d1) / (S * sigma * np.sqrt(T))


def vega(S: ArrayLike, K: ArrayLike, T: ArrayLike, r: ArrayLike, sigma: ArrayLike) -> FloatOrArray:
    """dV/dsigma, per 1.00 of volatility. Same for calls and puts."""
    S, K, T, r, sigma, d1, _d2 = _prepare(S, K, T, r, sigma)
    return S * norm.pdf(d1) * np.sqrt(T)


def theta(
    S: ArrayLike, K: ArrayLike, T: ArrayLike, r: ArrayLike, sigma: ArrayLike, kind: Kind = "call"
) -> FloatOrArray:
    """dV/dt per year as calendar time passes (= -dV/dT). Usually negative ("time decay")."""
    phi = _phi(kind)
    S, K, T, r, sigma, d1, d2 = _prepare(S, K, T, r, sigma)
    decay = -S * norm.pdf(d1) * sigma / (2 * np.sqrt(T))  # loss of optionality, same for both
    return decay - phi * r * K * np.exp(-r * T) * norm.cdf(phi * d2)


def rho(
    S: ArrayLike, K: ArrayLike, T: ArrayLike, r: ArrayLike, sigma: ArrayLike, kind: Kind = "call"
) -> FloatOrArray:
    """dV/dr, per 1.00 of interest rate. Call K T e^{-rT} N(d2), put -K T e^{-rT} N(-d2)."""
    phi = _phi(kind)
    S, K, T, r, sigma, _d1, d2 = _prepare(S, K, T, r, sigma)
    return phi * K * T * np.exp(-r * T) * norm.cdf(phi * d2)
