"""Regenerate the Chapter 02 figures and results tables.

Run from the repo root:  python -m docs.figures.make_ch02
(or `pip install -e .` once, then `python docs/figures/make_ch02.py` works too)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mc_engine import PricingResult
from mc_engine.analytics import black_scholes as bs

OUT = Path(__file__).parent
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]  # one per maturity
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 10,
    "lines.linewidth": 2, "legend.frameon": False,
})

K, r, sigma = 100.0, 0.05, 0.2
MATURITIES = [1.0, 0.5, 0.1]
S = np.linspace(50, 150, 401)

# ---- prices vs spot ---------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
for ax, kind, payoff in [(axes[0], "call", np.maximum(S - K, 0)),
                         (axes[1], "put", np.maximum(K - S, 0))]:
    for T, c in zip(MATURITIES, COLORS):
        ax.plot(S, bs.price(S, K, T, r, sigma, kind), color=c, label=f"T = {T}")
    ax.plot(S, payoff, color=INK2, ls="--", lw=1, label="payoff at maturity")
    ax.set(title=f"Black-Scholes {kind} price (K = 100, r = 5%, σ = 20%)",
           xlabel="Spot price S", ylabel="Option value")
    ax.legend(loc="upper left" if kind == "call" else "upper right")
fig.tight_layout()
fig.savefig(OUT / "02_prices.png", dpi=150)
plt.close(fig)

# ---- Greeks vs spot (call) ----------------------------------------------------
panels = [
    ("Delta  ∂V/∂S", lambda T: bs.delta(S, K, T, r, sigma)),
    ("Gamma  ∂²V/∂S²", lambda T: bs.gamma(S, K, T, r, sigma)),
    ("Vega  ∂V/∂σ (per 1.00)", lambda T: bs.vega(S, K, T, r, sigma)),
    ("Theta  ∂V/∂t (per year)", lambda T: bs.theta(S, K, T, r, sigma)),
    ("Rho  ∂V/∂r (per 1.00)", lambda T: bs.rho(S, K, T, r, sigma)),
]
fig, axes = plt.subplots(2, 3, figsize=(11, 6.2))
for ax, (title, f) in zip(axes.flat, panels):
    for T, c in zip(MATURITIES, COLORS):
        ax.plot(S, f(T), color=c, label=f"T = {T}")
    ax.axvline(K, color=INK2, ls=":", lw=1)
    ax.set(title=title, xlabel="Spot price S")
legend_ax = axes.flat[-1]
legend_ax.axis("off")
legend_ax.legend(*axes.flat[0].get_legend_handles_labels(), loc="center", title="Maturity (years)")
fig.suptitle("Call Greeks vs spot (K = 100 dotted, r = 5%, σ = 20%)")
fig.tight_layout()
fig.savefig(OUT / "02_greeks.png", dpi=150)
plt.close(fig)

# ---- results tables (printed as Markdown) ------------------------------------
args = (100, 100, 1.0, 0.05, 0.2)
print("| | Call | Put |\n|---|---:|---:|")
print(f"| Price | {bs.price(*args):.4f} | {bs.price(*args, kind='put'):.4f} |")
for name, f in [("Delta", bs.delta), ("Theta (per year)", bs.theta), ("Rho", bs.rho)]:
    print(f"| {name} | {f(*args):.4f} | {f(*args, kind='put'):.4f} |")
for name, f in [("Gamma", bs.gamma), ("Vega", bs.vega)]:
    print(f"| {name} | {f(*args):.4f} | {f(*args):.4f} |")

print("\n| Option | K | Black-Scholes | MC (10^6 paths) | SE | (MC - BS)/SE |")
print("|---|---:|---:|---:|---:|---:|")
S0, T = 100.0, 1.0
Z = np.random.default_rng(42).standard_normal(1_000_000)
ST = S0 * np.exp((r - 0.5 * sigma**2) * T + sigma * np.sqrt(T) * Z)
for kind in ["call", "put"]:
    for strike in [80, 100, 120]:
        payoff = np.maximum(ST - strike, 0) if kind == "call" else np.maximum(strike - ST, 0)
        mc = PricingResult.from_samples(np.exp(-r * T) * payoff)
        exact = bs.price(S0, strike, T, r, sigma, kind)
        print(f"| {kind} | {strike} | {exact:.4f} | {mc.price:.4f} | {mc.std_error:.4f} "
              f"| {(mc.price - exact) / mc.std_error:+.2f} |")
