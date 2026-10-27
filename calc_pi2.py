"""
Five ways to estimate pi with Monte Carlo, from simplest to cleverest.

Every method estimates the SAME thing: the area of a quarter circle of
radius 1, which is pi/4. Multiply that area by 4 and you get pi.

The methods differ only in HOW they pick points:
  1. plain        - random darts, count hits          (error ~ 1/sqrt(N))
  2. averaging    - repeat method 1 many times         (same as more darts)
  3. conditional  - pick only x, use the circle height (less noise)
  4. stratified   - one point in each of N slices      (error ~ 1/N^1.5)
  5. sobol        - evenly spread "smart" points       (error ~ 1/N)

Run:  python pi_methods.py
Needs: numpy and scipy only for method 5  (pip install scipy)
"""
import math
import random


# ---------------------------------------------------------------------------
# 1. PLAIN MONTE CARLO
# ---------------------------------------------------------------------------
# Throw darts at the square [0,1] x [0,1]. A dart is inside the quarter
# circle when x^2 + y^2 <= 1. The quarter circle covers pi/4 of the square,
# so  (hits / darts) ~ pi/4   ->   pi ~ 4 * hits / darts.
#
# Weakness: every dart only answers yes (1) or no (0). That's very little
# information per dart, so you need lots of darts.
def pi_plain(n):
    hits = 0
    for _ in range(n):
        x = random.random()          # random number in [0, 1)
        y = random.random()
        if x * x + y * y <= 1:       # no sqrt needed: compare squared distance
            hits += 1
    return 4 * hits / n


# ---------------------------------------------------------------------------
# 2. AVERAGING REPEATED RUNS
# ---------------------------------------------------------------------------
# Run method 1 `runs` times, each with n darts, and take the average.
# This reduces error by sqrt(runs) -- but it uses runs * n darts in total,
# so it is exactly as good as ONE run with runs * n darts. Not a new trick,
# just another way to spend more samples. (Arithmetic mean is the right
# choice; geometric/harmonic means give practically the same number here.)
def pi_averaged(n, runs):
    total = 0
    for _ in range(runs):
        total += pi_plain(n)
    return total / runs


# ---------------------------------------------------------------------------
# 3. CONDITIONAL (SAMPLE ONLY x)
# ---------------------------------------------------------------------------
# For a given x, the circle reaches up to height sqrt(1 - x^2).
# In method 1, y is a coin flip: "is y below that height?"
# Instead of flipping the coin, use the exact answer -- the height itself.
#
# The average height of the curve over [0,1] equals the area = pi/4,
# so  pi ~ 4 * average of sqrt(1 - x^2).
#
# Why better: the randomness from y is gone, and we need only one random
# number per sample instead of two.
def pi_conditional(n):
    total = 0
    for _ in range(n):
        x = random.uniform(0, 1)
        total += math.sqrt(1 - x * x)    # height of circle at this x
    return 4 * total / n


# ---------------------------------------------------------------------------
# 4. STRATIFIED SAMPLING
# ---------------------------------------------------------------------------
# Random points clump together in some places and leave gaps in others.
# Fix: cut [0,1] into n equal slices and put exactly ONE random point in
# each slice. Slice i covers [i/n, (i+1)/n), so a random point inside it is
#       x = (i + random()) / n
#
# Then use the circle height, as in method 3.
#
# Why so good: inside one thin slice the curve is nearly flat, so each
# slice's estimate is almost exact, and coverage is always even.
def pi_stratified(n):
    total = 0
    for i in range(n):
        x = (i + random.uniform(0, 1)) / n    # one random point in slice i
        total += math.sqrt(1 - x * x)
    return 4 * total / n


# ---------------------------------------------------------------------------
# 5. SOBOL (LOW-DISCREPANCY SEQUENCE)
# ---------------------------------------------------------------------------
# A Sobol sequence is a special list of points designed so each new point
# lands in the biggest remaining gap. The square gets filled evenly, much
# more evenly than with random points. "scramble=True" adds a little
# randomness so different seeds give slightly different (still even) points.
#
# We use the plain hit-counting test from method 1 -- only the points change.
# Tip: use a power of two for n (1024, 16384, 131072, ...).
#
# Why useful: unlike stratifying by hand, this works in many dimensions,
# which is why it's used in finance, graphics and physics.
def pi_sobol(n, seed=0):
    from scipy.stats import qmc
    points = qmc.Sobol(d=2, scramble=True, seed=seed).random(n)  # n x 2 array
    hits = 0
    for x, y in points:
        if x * x + y * y <= 1:
            hits += 1
    return 4 * hits / n


# ---------------------------------------------------------------------------
# Compare them all
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    random.seed(42)
    n = 2 ** 17          # 131,072 samples (a power of two, for Sobol)

    results = [
        ("1. plain",          pi_plain(n)),
        ("2. averaged (x10)", pi_averaged(n // 10, 10)),   # same total samples
        ("3. conditional",    pi_conditional(n)),
        ("4. stratified",     pi_stratified(n)),
        ("5. sobol",          pi_sobol(n)),
    ]

    print(f"true pi = {math.pi:.10f}   (N = {n:,} samples each)\n")
    print(f"{'method':<20}{'estimate':>16}{'error':>14}")
    for name, est in results:
        print(f"{name:<20}{est:>16.10f}{abs(est - math.pi):>14.2e}")