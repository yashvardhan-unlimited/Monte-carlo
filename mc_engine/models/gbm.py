"""Geometric Brownian motion under the risk-neutral measure.

    dS = r S dt + sigma S dW

simulated with the exact log-space step (derived in docs/03_gbm_and_first_mc.md):

    S(t+dt) = S(t) * exp((r - sigma^2/2) dt + sigma sqrt(dt) Z),   Z ~ N(0, 1).

"Exact" means the simulated prices at the grid times have exactly the right joint
distribution for any step size, so there is no discretisation error to manage.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class GBM:
    """GBM model with constant rate and volatility.

    Attributes:
        S0: Spot price today.
        r: Continuously compounded risk-free rate (also the drift under Q).
        sigma: Volatility per sqrt(year), decimal (0.2 = 20%).
    """

    S0: float
    r: float
    sigma: float

    def __post_init__(self) -> None:
        if self.S0 <= 0 or self.sigma < 0:
            raise ValueError("need S0 > 0 and sigma >= 0")

    def simulate(
        self, n_paths: int, n_steps: int, T: float, rng: np.random.Generator
    ) -> NDArray[np.float64]:
        """Simulate price paths on an equally spaced grid t_j = j * T / n_steps.

        Args:
            n_paths: Number of independent paths (rows).
            n_steps: Number of time steps; output has n_steps + 1 columns (column 0 is S0).
            T: Horizon in years.
            rng: NumPy random generator; pass a seeded one for reproducible results.

        Returns:
            Array of shape (n_paths, n_steps + 1).
        """
        if n_steps < 1 or T <= 0:
            raise ValueError("need n_steps >= 1 and T > 0")
        dt = T / n_steps
        Z = rng.standard_normal((n_paths, n_steps))
        # Work with log-prices: the log-increments are i.i.d. normal, so a cumulative sum
        # gives log(S_t / S0) at every grid time in one vectorised call. Adding logs
        # instead of multiplying prices also avoids accumulating rounding error.
        log_increments = (self.r - 0.5 * self.sigma**2) * dt + self.sigma * np.sqrt(dt) * Z
        paths = np.empty((n_paths, n_steps + 1))
        paths[:, 0] = self.S0
        paths[:, 1:] = self.S0 * np.exp(np.cumsum(log_increments, axis=1))
        return paths
