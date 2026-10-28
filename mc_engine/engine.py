"""Core engine: the result type and the generic Monte Carlo pricer.

Every Monte Carlo pricer in this project returns a ``PricingResult``. Having one
result type means tests, benchmarks and reports can treat all pricers the same way.

``price`` knows nothing about any particular model or payoff. It only relies on the
contract from docs/00_architecture.md:
    model.simulate(n_paths, n_steps, T, rng) -> paths, shape (n_paths, n_steps + 1)
    payoff(paths) -> one cash flow per path, shape (n_paths,)
"""

from __future__ import annotations

import time
from collections.abc import Callable
from dataclasses import dataclass
from typing import Protocol

import numpy as np
from numpy.typing import ArrayLike, NDArray

# z such that P(-z < Z < z) = 0.95 for a standard normal Z, i.e. Phi^{-1}(0.975).
Z_95 = 1.959963984540054


@dataclass(frozen=True)
class PricingResult:
    """Outcome of one Monte Carlo pricing run.

    Attributes:
        price: The Monte Carlo estimate (sample mean of discounted payoffs).
        std_error: Standard error of that estimate, s / sqrt(n).
        ci_95: 95% confidence interval (low, high) = price -/+ 1.96 * std_error.
        n_paths: Number of simulated paths (samples) used.
        runtime_sec: Wall-clock time the pricer took, in seconds.
    """

    price: float
    std_error: float
    ci_95: tuple[float, float]
    n_paths: int
    runtime_sec: float

    @classmethod
    def from_samples(cls, samples: ArrayLike, runtime_sec: float = 0.0) -> PricingResult:
        """Build a result from per-path discounted payoffs.

        Every pricer ends the same way (average, standard error, CI), so that
        arithmetic lives here once instead of being copy-pasted into each pricer.

        Args:
            samples: 1-D array of discounted payoffs, one per simulated path.
            runtime_sec: Time the simulation took, measured by the caller.

        Returns:
            The filled-in PricingResult.

        Raises:
            ValueError: If fewer than 2 samples are given (sample std is undefined).
        """
        x = np.asarray(samples, dtype=float).ravel()
        n = x.size
        if n < 2:
            raise ValueError(f"need at least 2 samples to estimate a standard error, got {n}")

        price = float(x.mean())
        # ddof=1 gives the unbiased sample variance (divide by n-1, not n).
        std_error = float(x.std(ddof=1) / np.sqrt(n))
        ci_95 = (price - Z_95 * std_error, price + Z_95 * std_error)
        return cls(price, std_error, ci_95, n, runtime_sec)


class Model(Protocol):
    """Anything with a risk-free rate and a path simulator counts as a model."""

    r: float

    def simulate(
        self, n_paths: int, n_steps: int, T: float, rng: np.random.Generator
    ) -> NDArray[np.float64]: ...


Payoff = Callable[[NDArray[np.float64]], NDArray[np.float64]]


def price(
    model: Model,
    payoff: Payoff,
    T: float,
    n_paths: int,
    n_steps: int = 1,
    seed: int = 42,
) -> PricingResult:
    """Monte Carlo price = discounted average payoff over simulated paths.

    Args:
        model: Path generator under the risk-neutral measure (e.g. GBM).
        payoff: Maps paths (n_paths, n_steps + 1) to cash flows at T (n_paths,).
        T: Maturity in years; payoffs are discounted with exp(-model.r * T).
        n_paths: Number of simulated paths.
        n_steps: Time steps per path. 1 is enough for payoffs that only use S_T.
        seed: Seed for the random generator, so every run is reproducible.

    Returns:
        PricingResult with price, standard error, 95% CI, n_paths and runtime.
    """
    start = time.perf_counter()
    rng = np.random.default_rng(seed)
    paths = model.simulate(n_paths, n_steps, T, rng)
    discounted = np.exp(-model.r * T) * payoff(paths)
    return PricingResult.from_samples(discounted, runtime_sec=time.perf_counter() - start)
