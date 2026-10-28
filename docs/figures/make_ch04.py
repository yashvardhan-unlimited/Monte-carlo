"""Chapter 04 convergence study: how Monte Carlo error shrinks with the number of paths.

Prices an at-the-money European call for N = 100 ... 1,000,000 paths, 200 independent
runs per N, and produces:
    04_error_vs_n.png         log(error) vs log(N) with fitted slope (exact GBM and a biased Euler pricer)
    04_confidence_interval.png the 95% CI shrinking with N
    04_coverage.png            50 confidence intervals at N = 1,000: how many contain the truth
plus the Markdown tables quoted in docs/04_convergence.md.

Run from the repo root:  python -m docs.figures.make_ch04     (takes about 30 s)
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mc_engine.analytics import black_scholes as bs
from mc_engine.engine import price
from mc_engine.models.gbm import GBM
from mc_engine.payoffs.european import European

OUT = Path(__file__).parent
sys.stdout.reconfigure(encoding="utf-8")  # Windows consoles default to cp1252 (no √, σ)
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 10,
    "lines.linewidth": 2, "legend.frameon": False,
})

S0, K, r, sigma, T = 100.0, 100.0, 0.05, 0.2, 1.0
MODEL, CALL = GBM(S0, r, sigma), European(K)
EXACT = bs.price(S0, K, T, r, sigma)
NS = np.logspace(2, 6, 9).round().astype(int)  # 100, 316, 1000, ..., 1,000,000
RUNS = 200  # independent repetitions per N


def seed_for(k: int, run: int) -> int:
    """A different seed for every (N, run) pair, so no two runs share random numbers."""
    return 10_000 * k + run


def euler_one_step(n: int, seed: int) -> float:
    """Call price with one Euler step S_T = S0 (1 + rT + sigma sqrt(T) Z): biased on purpose."""
    Z = np.random.default_rng(seed).standard_normal(n)
    ST = S0 * (1 + r * T + sigma * np.sqrt(T) * Z)
    return float(np.exp(-r * T) * np.maximum(ST - K, 0).mean())


# ---- run the experiment ---------------------------------------------------------
est = np.empty((len(NS), RUNS))
se = np.empty((len(NS), RUNS))
euler = np.empty((len(NS), RUNS))
runtime = np.empty(len(NS))
for k, n in enumerate(NS):
    results = [price(MODEL, CALL, T, int(n), seed=seed_for(k, i)) for i in range(RUNS)]
    est[k] = [x.price for x in results]
    se[k] = [x.std_error for x in results]
    runtime[k] = np.median([x.runtime_sec for x in results])
    euler[k] = [euler_one_step(int(n), seed_for(k, i) + 5_000) for i in range(RUNS)]

err = est - EXACT
rmse = np.sqrt((err**2).mean(axis=1))
rmse_euler = np.sqrt(((euler - EXACT) ** 2).mean(axis=1))
covered = np.abs(err) < 1.96 * se
sigma_hat = (se * np.sqrt(NS)[:, None]).mean()  # std dev of one discounted payoff
euler_bias = (euler[-1] - EXACT).mean()

slope, intercept = np.polyfit(np.log10(NS), np.log10(rmse), 1)
single = np.abs(err[:, 0])  # one run per N, as you'd get without repetitions
slope_single, _ = np.polyfit(np.log10(NS), np.log10(single), 1)

# ---- figure 1: log error vs log N -------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
show = 20  # individual runs per N to draw (all 200 would be a solid smear)
ax.scatter(np.repeat(NS, show), np.abs(err[:, :show]).ravel(), s=8, color=BLUE, alpha=0.25,
           linewidths=0, label="|error| of single runs (exact GBM)")
ax.plot(NS, rmse, "o", color=BLUE, ms=7, label="RMSE over 200 runs (exact GBM)")
ax.plot(NS, 10**intercept * NS**slope, color=BLUE, lw=1.5,
        label=f"fitted line: slope = {slope:.3f}")
ax.plot(NS, sigma_hat / np.sqrt(NS), color=INK, ls="--", lw=1.2, label="theory  σ_X / √N")
ax.plot(NS, rmse_euler, "s", color=ORANGE, ms=6, label="RMSE, one-step Euler (biased)")
ax.axhline(abs(euler_bias), color=ORANGE, ls=":", lw=1.2, label=f"Euler bias {abs(euler_bias):.3f}")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set(xlabel="Number of paths N (log scale)", ylabel="Pricing error (log scale)",
       title="Monte Carlo error vs number of paths: ATM call, Black-Scholes 10.4506")
ax.legend(loc="lower left", fontsize=9)
fig.tight_layout()
fig.savefig(OUT / "04_error_vs_n.png", dpi=150)
plt.close(fig)

# ---- figure 2: the confidence interval shrinking -----------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
ax1.errorbar(NS, est[:, 0], yerr=1.96 * se[:, 0], fmt="o", color=BLUE, ms=5, capsize=4,
             elinewidth=1.5, label="estimate ± 1.96 SE (one run per N)")
ax1.axhline(EXACT, color=INK, ls="--", lw=1.2, label="Black-Scholes")
ax1.set_xscale("log")
ax1.set(xlabel="Number of paths N (log scale)", ylabel="Call price",
        title="95% confidence interval vs N")
ax1.legend(loc="upper right", fontsize=9)
half = 1.96 * se.mean(axis=1)
ax2.plot(NS, half, "o", color=BLUE, ms=6, label="measured half-width 1.96·SE")
ax2.plot(NS, 1.96 * sigma_hat / np.sqrt(NS), color=INK, ls="--", lw=1.2,
         label="1.96 σ_X / √N")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set(xlabel="Number of paths N (log scale)", ylabel="CI half-width (log scale)",
        title="Half-width falls 10× per 100× paths")
ax2.legend(loc="upper right", fontsize=9)
fig.tight_layout()
fig.savefig(OUT / "04_confidence_interval.png", dpi=150)
plt.close(fig)

# ---- figure 3: coverage at N = 1,000 -----------------------------------------------
k = int(np.where(NS == 1000)[0][0])
fig, ax = plt.subplots(figsize=(8, 6))
for i in range(50):
    lo, hi = est[k, i] - 1.96 * se[k, i], est[k, i] + 1.96 * se[k, i]
    hit = covered[k, i]
    ax.plot([lo, hi], [i, i], color=BLUE if hit else ORANGE, lw=2 if hit else 3)
    ax.plot(est[k, i], i, "o" if hit else "x", color=BLUE if hit else ORANGE, ms=4 if hit else 9, mew=2.5)
ax.axvline(EXACT, color=INK, ls="--", lw=1.2)
ax.plot([], [], color=BLUE, lw=2, label="contains the true price")
ax.plot([], [], color=ORANGE, lw=3, marker="x", ms=9, mew=2.5, label="misses the true price")
ax.set(xlabel="Call price", ylabel="Run number", yticks=[],
       title=f"50 independent 95% CIs at N = 1,000 "
             f"({covered[k, :50].sum()} of 50 contain 10.4506; {covered[k].mean():.1%} of all 200)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2)
fig.tight_layout()
fig.savefig(OUT / "04_coverage.png", dpi=150)
plt.close(fig)

# ---- tables --------------------------------------------------------------------
print(f"Black-Scholes {EXACT:.4f}; std dev of one discounted payoff sigma_X = {sigma_hat:.3f}")
print(f"fitted slope (RMSE over {RUNS} runs): {slope:.4f};  fitted slope (single runs): {slope_single:.4f}")
print(f"Euler one-step bias {euler_bias:+.4f}\n")
print("| N | mean SE | RMSE | RMSE·√N | 95% CI coverage | Euler RMSE | runtime (s) |")
print("|---:|---:|---:|---:|---:|---:|---:|")
for k, n in enumerate(NS):
    print(f"| {n:,} | {se[k].mean():.4f} | {rmse[k]:.4f} | {rmse[k] * np.sqrt(n):.2f} "
          f"| {covered[k].mean():.1%} | {rmse_euler[k]:.4f} | {runtime[k]:.4f} |")

sec_per_path = runtime[-1] / NS[-1]
print("\n| Target SE | Paths needed (σ_X/SE)² | Time at measured speed |\n|---:|---:|---:|")
for target in [0.1, 0.01, 0.001, 0.0001]:
    n_needed = (sigma_hat / target) ** 2
    secs = n_needed * sec_per_path
    when = f"{secs:.3f} s" if secs < 60 else f"{secs / 60:.1f} min" if secs < 3600 else f"{secs / 3600:.1f} h"
    print(f"| {target} | {n_needed:,.0f} | {when} |")
