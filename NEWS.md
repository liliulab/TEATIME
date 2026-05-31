# TEATIME 2.5.7

Author: Hai Chen.

## Execution modes

`TEATIME.run()` now exposes three execution modes:

1. **default** (`fast_version = FALSE`) — the exact reference pipeline (identical
   logic and step order to the original production code).
2. **fast** (`fast_version = TRUE`) — the same algorithm with C++/vectorised
   kernels. **Bit-identical** to `default` given the same seed (verified
   `max|diff| = 0` across MMRF samples and the bundled example data); ~6× faster.
3. **approximation** (`fast_version = TRUE` + `options(teatime.approx = TRUE)`) —
   a **deterministic**, faster mode. Speed comes from a *mathematical*
   substitution, not from running fewer simulations: the Monte-Carlo
   `round(rbinom(n, depth, vaf) / depth)` in `peak_test` is replaced by its exact
   analytic distribution `dbinom(0:depth, depth, vaf)` (the n → ∞ limit, noise
   free). The residual stochastic tie-break choices are pinned with a fixed
   internal seed, so the result is repeatable run-to-run with no user seed.

## Behaviour changes

- **Unseeded by default.** `seed` now defaults to `NA`, matching the original
  production runs (which set no seed). Each unseeded run may differ slightly, as
  before. Pass an explicit `seed` for reproducibility or for `fast == default`
  bit-identity checks.

## Bug fixes

- **R ≥ 4.2 inter-branch crash.** In `normalrun`, the guard `if(!is.na(fit.check$mu))`
  evaluated a length > 1 condition. Under R < 4.2 this only warned (using the
  first element); under R ≥ 4.2 it is a fatal error, which was caught and
  silently collapsed the inter branch to `NA`. Fixed to `if(!is.na(fit.check$mu[1]))`,
  preserving the pre-4.2 first-element semantics so behaviour matches the original
  code on modern R.

## Notes

- The `mu_real > 3` candidate filter is kept strict (`>`), faithful to the
  original production code.
