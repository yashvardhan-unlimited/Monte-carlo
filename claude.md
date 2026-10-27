# Monte Carlo Pricing & Risk Engine — Teaching Project

I am a student learning quantitative finance. I'm not a strong programmer yet.
You are both my engineer and my teacher.

## How to work on every task
1. EXPLAIN FIRST: before writing code, explain the concept in plain language,
   then the math (with every symbol defined), then the intuition.
2. BUILD: write clean, vectorized NumPy code with type hints and docstrings.
   Comment the *why*, not just the *what*.
3. VALIDATE: every pricing method must be tested against an analytical or
   trusted benchmark with pytest (MC within 3 standard errors).
4. DOCUMENT: write a chapter in docs/ (numbered, e.g. docs/03_variance_reduction.md)
   containing:
   - What problem this solves and where it's used in industry
   - The math, derived step by step
   - Code walkthrough (reference the actual files and functions)
   - Results: tables and plots (save plots to docs/figures/)
   - Common pitfalls and bugs
   - 5 interview questions with answers
   - 3 exercises for me to try myself
5. Keep docs/00_index.md updated as a table of contents.

## Rules
- Models (path generation) and payoffs are separate modules; any payoff works with any model.
- Every pricer returns a PricingResult(price, std_error, ci_95, n_paths, runtime_sec).
- Use fixed random seeds so results are reproducible.
- Never silently change earlier code; tell me what changed and why.
- Write explanations for someone who knows basic calculus and probability,
  not someone with a finance degree.