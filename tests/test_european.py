import numpy as np
import pytest

from mc_engine.payoffs.european import European

# Three paths of two steps; only the last column (80, 100, 120) should matter.
PATHS = np.array([[100.0, 150.0, 80.0], [100.0, 50.0, 100.0], [100.0, 90.0, 120.0]])


def test_call_payoff():
    np.testing.assert_array_equal(European(K=100)(PATHS), [0, 0, 20])


def test_put_payoff():
    np.testing.assert_array_equal(European(K=100, kind="put")(PATHS), [20, 0, 0])


def test_rejects_unknown_kind():
    with pytest.raises(ValueError):
        European(K=100, kind="straddle")
