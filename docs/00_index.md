# Monte Carlo Pricing & Risk Engine — Table of Contents

Each chapter explains a concept, derives the math, walks through the code, shows
results, and ends with interview questions and exercises.

| #  | Chapter | What you learn | Code | Status |
|----|---------|----------------|------|--------|
| 00 | [Architecture](00_architecture.md) | How the engine is organised: models × payoffs × engine | `mc_engine/engine.py` | Done |
| 01 | [Foundations](01_foundations.md) | Stocks, options, payoffs, discounting; expectation, variance, normal distribution; LLN, CLT, standard error; Brownian motion; risk-neutral pricing | none (theory only) | Done |
| 02 | [Black-Scholes formula & Greeks](02_black_scholes.md) | Deriving the formula; meaning of each input and term; delta, gamma, vega, theta, rho; the benchmark for everything else | `analytics/black_scholes.py` | Done |
| 03 | [GBM simulation & first Monte Carlo pricer](03_gbm_and_first_mc.md) | GBM, volatility drag and Itô's lemma; exact vs Euler discretisation; the model/payoff/engine interfaces; matching Black-Scholes on a grid | `models/gbm.py`, `payoffs/european.py`, `engine.py` | Done |
| 04 | [Convergence: why error shrinks like 1/√N](04_convergence.md) | Standard error and confidence intervals; the 1/√N law measured on a log-log plot; CI coverage; bias vs noise; the cost of accuracy | `docs/figures/make_ch04.py` | Done |
| 05 | Variance reduction | Antithetic variates, control variates, importance sampling | `variance_reduction/` | Planned |
| 06 | Quasi-Monte Carlo | Sobol sequences, Brownian bridge, why QMC converges faster | `variance_reduction/` | Planned |
| 07 | Path-dependent payoffs | Asian, barrier and lookback options; discretisation bias | `payoffs/` | Planned |
| 08 | Greeks by Monte Carlo | Bump-and-revalue, pathwise derivative, likelihood ratio method | `greeks/` | Planned |
| 09 | Multiple correlated assets | Cholesky decomposition, basket and spread options | `models/`, `payoffs/` | Planned |
| 10 | Stochastic volatility: Heston | Simulating variance that moves; Euler vs. better schemes; Fourier benchmark | `models/`, `analytics/` | Planned |
| 11 | Jumps: Merton jump-diffusion | Poisson jumps in prices; Merton series benchmark | `models/`, `analytics/` | Planned |
| 12 | American options: Longstaff-Schwartz | Early exercise via regression; lower/upper bounds | `payoffs/`, `engine.py` | Planned |
| 13 | Performance | Vectorisation vs. loops vs. Numba; benchmarking | `benchmarks/` | Planned |
| 14 | Portfolio risk: VaR & Expected Shortfall | Simulating P&L; risk measures; backtesting | `risk/` | Planned |

Figures live in [`figures/`](figures/); from Chapter 02 on, `figures/make_chNN.py` regenerates each chapter's figures. Generated reports go in `../reports/`.
