"""European call and put payoffs: depend only on the final price S_T."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class European:
    """Pays max(S_T - K, 0) (call) or max(K - S_T, 0) (put) at maturity.

    Call it on a paths array of shape (n_paths, n_steps + 1); it reads only the last
    column, so it works with paths from any model and any number of steps.
    """

    K: float
    kind: Literal["call", "put"] = "call"

    def __post_init__(self) -> None:
        if self.kind not in ("call", "put"):
            raise ValueError(f"kind must be 'call' or 'put', got {self.kind!r}")

    def __call__(self, paths: NDArray[np.float64]) -> NDArray[np.float64]:
        ST = paths[:, -1]  # last column = price at maturity on every path (not paths[-1]!)
        if self.kind == "call":
            return np.maximum(ST - self.K, 0.0)
        return np.maximum(self.K - ST, 0.0)
