# Monte Carlo Pricing & Risk Engine — Table of Contents

Each chapter explains a concept, derives the math, walks through the code, shows
results, and ends with interview questions and exercises.

| #  | Chapter | What you learn | Code | Status |
|----|---------|----------------|------|--------|
| 00 | [Architecture](00_architecture.md) | How the engine is organised: models × payoffs × engine | `mc_engine/engine.py` | Done |
| 01 | Monte Carlo foundations | Law of large numbers, central limit theorem, standard error, confidence intervals, seeds | `engine.py` | Planned |
| 02 | Geometric Brownian Motion & Black-Scholes | Simulating stock paths; pricing a European call; checking against the Black-Scholes formula | `models/`, `payoffs/`, `analytics/` | Planned |
| 03 | Variance reduction | Antithetic variates, control variates, importance sampling | `variance_reduction/` | Planned |
| 04 | Quasi-Monte Carlo | Sobol sequences, Brownian bridge, why QMC converges faster | `variance_reduction/` | Planned |
| 05 | Path-dependent payoffs | Asian, barrier and lookback options; discretisation bias | `payoffs/` | Planned |
| 06 | Greeks | Bump-and-revalue, pathwise derivative, likelihood ratio method | `greeks/` | Planned |
| 07 | Multiple correlated assets | Cholesky decomposition, basket and spread options | `models/`, `payoffs/` | Planned |
| 08 | Stochastic volatility: Heston | Simulating variance that moves; Euler vs. better schemes; Fourier benchmark | `models/`, `analytics/` | Planned |
| 09 | Jumps: Merton jump-diffusion | Poisson jumps in prices; Merton series benchmark | `models/`, `analytics/` | Planned |
| 10 | American options: Longstaff-Schwartz | Early exercise via regression; lower/upper bounds | `payoffs/`, `engine.py` | Planned |
| 11 | Performance | Vectorisation vs. loops vs. Numba; benchmarking | `benchmarks/` | Planned |
| 12 | Portfolio risk: VaR & Expected Shortfall | Simulating P&L; risk measures; backtesting | `risk/` | Planned |

Figures live in [`figures/`](figures/). Generated reports go in `../reports/`.
