"""Regenerate the Chapter 03 figures and results tables.

Run from the repo root:  python -m docs.figures.make_ch03
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import lognorm

from mc_engine.analytics import black_scholes as bs
from mc_engine.engine import price
from mc_engine.models.gbm import GBM
from mc_engine.payoffs.european import European

OUT = Path(__file__).parent
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 10,
    "lines.linewidth": 2, "legend.frameon": False,
})

S0, r, sigma, T = 100.0, 0.05, 0.2, 1.0
model = GBM(S0, r, sigma)
drift = r - 0.5 * sigma**2  # median growth rate of the log-price

# ---- 50 sample paths ----------------------------------------------------------
n_steps = 252
t = np.linspace(0, T, n_steps + 1)
paths = model.simulate(50, n_steps, T, np.random.default_rng(2024))
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(t, paths.T, color=BLUE, lw=0.7, alpha=0.45)
ax.plot([], [], color=BLUE, lw=0.7, alpha=0.6, label="50 simulated paths")
ax.plot(t, S0 * np.exp(r * t), color=INK, lw=2, label="mean  S0·exp(r t)")
ax.plot(t, S0 * np.exp(drift * t), color=INK, lw=1.5, ls="--", label="median  S0·exp((r − σ²/2) t)")
for sign in (1, -1):
    ax.plot(t, S0 * np.exp(drift * t + sign * 1.96 * sigma * np.sqrt(t)), color=ORANGE, lw=1.5,
            ls=":", label="95% band" if sign == 1 else None)
ax.set(xlabel="Time t (years)", ylabel="Price S_t",
       title="GBM paths: S0 = 100, r = 5%, σ = 20%, 252 daily steps")
ax.legend(loc="upper left")
fig.tight_layout()
fig.savefig(OUT / "03_paths.png", dpi=150)
plt.close(fig)

# ---- terminal distribution vs lognormal ------------------------------------
ST = model.simulate(100_000, 1, T, np.random.default_rng(7))[:, -1]
theory = lognorm(s=sigma * np.sqrt(T), scale=S0 * np.exp(drift * T))
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.hist(ST, bins=120, range=(40, 220), density=True, color=BLUE, edgecolor=SURF, linewidth=0.5,
        label="simulated S_T (100,000 paths)")
x = np.linspace(40, 220, 600)
ax.plot(x, theory.pdf(x), color=INK, lw=2, label="lognormal density (theory)")
mean, median = S0 * np.exp(r * T), S0 * np.exp(drift * T)
ax.axvline(median, color=ORANGE, ls="--", lw=1.5, label=f"median {median:.2f}")
ax.axvline(mean, color=ORANGE, ls="-", lw=1.5, label=f"mean {mean:.2f}")
ax.set(xlabel="Terminal price S_T", ylabel="Density",
       title="Terminal prices at T = 1 vs the theoretical lognormal")
ax.legend(loc="upper right")
fig.tight_layout()
fig.savefig(OUT / "03_terminal_hist.png", dpi=150)
plt.close(fig)

# ---- convergence of the price estimate ------------------------------------
exact = bs.price(S0, 100, T, r, sigma)
ns = np.unique(np.logspace(2, 6, 25).astype(int))
res = [price(model, European(100), T, n_paths=int(n), seed=0) for n in ns]
est = np.array([x.price for x in res])
half = 1.96 * np.array([x.std_error for x in res])
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.fill_between(ns, est - half, est + half, color=BLUE, alpha=0.2, linewidth=0, label="95% CI")
ax.plot(ns, est, color=BLUE, marker="o", ms=4, label="Monte Carlo estimate")
ax.axhline(exact, color=INK, ls="--", lw=1.5, label=f"Black-Scholes {exact:.4f}")
ax.set_xscale("log")
ax.set(xlabel="Number of paths (log scale)", ylabel="Call price",
       title="ATM call (K = 100, T = 1): estimate and 95% CI vs number of paths")
ax.legend(loc="upper right")
fig.tight_layout()
fig.savefig(OUT / "03_convergence.png", dpi=150)
plt.close(fig)

# ---- table: grid vs Black-Scholes (same settings as the tests) --------------
print("| T | K | Call BS | Call MC | z | Put BS | Put MC | z |")
print("|---:|---:|---:|---:|---:|---:|---:|---:|")
zs = []
for T_ in [0.25, 1.0, 2.0]:
    for K in [80, 90, 100, 110, 120]:
        row = f"| {T_} | {K} |"
        for kind in ["call", "put"]:
            for steps in [1, 12]:
                mc = price(model, European(K, kind), T_, n_paths=200_000, n_steps=steps)
                z = (mc.price - bs.price(S0, K, T_, r, sigma, kind)) / mc.std_error
                zs.append(z)
                if steps == 1:
                    row += f" {bs.price(S0, K, T_, r, sigma, kind):.4f} | {mc.price:.4f} | {z:+.2f} |"
        print(row)
zs = np.abs(zs)
print(f"\n{len(zs)} checks: within 1 SE {np.mean(zs < 1):.0%}, 2 SE {np.mean(zs < 2):.0%}, "
      f"3 SE {np.mean(zs < 3):.0%}, max |z| = {zs.max():.2f}")

# ---- table: exact step vs Euler step -----------------------------------------
print("\n| Steps | Exact price | z vs BS | Euler price | z vs BS |\n|---:|---:|---:|---:|---:|")
n = 1_000_000
for steps in [1, 2, 4, 12, 52, 252]:
    rng = np.random.default_rng(11)
    dt = T / steps
    s_exact = np.full(n, S0)
    s_euler = np.full(n, S0)
    for _ in range(steps):  # loop over time, vectorised over paths: memory stays O(n)
        Z = rng.standard_normal(n)  # same Z for both schemes, so only the scheme differs
        s_exact *= np.exp(drift * dt + sigma * np.sqrt(dt) * Z)
        s_euler *= 1 + r * dt + sigma * np.sqrt(dt) * Z
    row = f"| {steps} |"
    for s in (s_exact, s_euler):
        x = np.exp(-r * T) * np.maximum(s - 100, 0)
        p, se = x.mean(), x.std(ddof=1) / np.sqrt(n)
        row += f" {p:.4f} | {(p - exact) / se:+.1f} |"
    print(row)
print(f"(Black-Scholes {exact:.4f}; SE about 0.0147 for every row)")

# ---- table: runtime --------------------------------------------------------
print("\n| n_paths | n_steps | runtime (s) |\n|---:|---:|---:|")
for n_paths, steps in [(100_000, 1), (1_000_000, 1), (100_000, 52), (100_000, 252)]:
    rt = min(price(model, European(100), T, n_paths, steps).runtime_sec for _ in range(3))
    print(f"| {n_paths:,} | {steps} | {rt:.3f} |")
