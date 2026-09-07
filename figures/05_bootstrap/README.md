# Iteration 05 — study-level bootstrap *(planned, not implemented)*

> **Nothing in this folder is computed yet.** This is a design note, kept here so the
> reasoning is visible before any code exists — the same way iterations 01 and 04 keep their
> failures. There is no `analysis_bootstrap.R`.

---

## Why this is the next step

Iterations 03 and 04 left a genuine impasse, and it is not one more modelling choice will
resolve:

| | autocorrelation handled | residuals normal |
| --- | :---: | :---: |
| [iteration 03](../03_natural_spline_hlme/) — `random = ~1` | ✗ ACF 0.978 | ✓ kurtosis 1.18 |
| [iteration 04](../04_natural_spline_hlme/) — `random = ~ns(…)` | mostly ✓ | ✗ kurtosis 12.55 |

Both were tested exhaustively against the obvious escapes and none worked: `cor = AR(TIME)`
is aliased with the random intercept; the spline df only slides you along the trade-off
(kurtosis 5.70 → 12.55 as ACF improves 0.938 → 0.836); trimming to the x-overlap halves the
kurtosis to 5.51 but no further.

The impasse only matters because **both models' standard errors are model-based**, and so
depend on assumptions the diagnostics reject. The way out is not a better likelihood — it is
**standard errors that do not depend on the residual distribution at all**.

There is a second, larger reason. The reviewer's first point
([`../../docs/STATISTICAL_METHODS_REVIEW.md`](../../docs/STATISTICAL_METHODS_REVIEW.md)) is
that **condition is perfectly confounded with source study** — no study measured both in vivo
and ex vivo. The effective replication is therefore the number of **studies (17)**, not
shoulders (44) and certainly not rows (3,724). Every Wald test computed so far is
pseudo-replicated. A study-level resampling scheme fixes the distributional problem and the
pseudo-replication problem with one mechanism, which is why it is worth doing before
anything else.

This is outstanding item 2 in
[`../../docs/STATISTICAL_METHODS_RESPONSE.md`](../../docs/STATISTICAL_METHODS_RESPONSE.md).

## What to estimate

**Not** the individual coefficients. The by-df Wald grid in iteration 04 showed those are
unstable: only `ns1:condin vivo` passes at every df, `condin vivo` fails at every df, and the
higher-order terms swap significance as the knots move. Bootstrapping an unstable quantity
would produce an honest interval for a number nobody should report.

Bootstrap the quantity that **has** been stable across every configuration tried:

$$
\text{ES} = \frac{1}{|\mathcal{X}|}\int_{\mathcal{X}} \bigl| \gamma_0 + \sum_k \gamma_k B_k(x) \bigr|\,dx
$$

the mean absolute in-vivo − ex-vivo difference over the overlap. It came out **9.9–10.5°**
across every df in iteration 03, **9.89–10.36°** across every df in iteration 04, and
unchanged between a random intercept and a random curve. That robustness is the finding; the
bootstrap supplies the interval it currently lacks.

Worth also carrying, since they are cheap once the loop exists: the **signed** mean
difference (the review objected to reporting only the absolute value), and the difference
curve at a few fixed elevations (say 30°, 60°, 90°, 120°) so the interval can be read
locally rather than only as one average.

## The resampling unit is the study

**Resample the 17 studies with replacement, not the 44 shoulders and never the rows.** This
is the whole point. Resampling rows would reproduce exactly the pseudo-replication the
exercise exists to remove; resampling shoulders would still treat two shoulders from one
study as independent evidence, which they are not when condition is a study-level property.

Practical consequences to handle:

- Studies contribute wildly unequal numbers of shoulders, so a resample can be badly
  unbalanced. Some replicates will contain **only one condition** — those must be detected
  and discarded, with the discard count reported, not silently dropped.
- `IDnum` must be **rebuilt inside each replicate**: a study drawn twice contributes two
  copies of its shoulders, and they must be distinct subjects to lcmm or the fit is wrong.
- The spline basis must use **the knots and boundary knots of the original fit**, not
  recomputed per replicate — otherwise each replicate estimates a slightly different
  quantity and the interval is not for anything.

## Suggested procedure

1. Fit the reference model once on the real data: diagonal random spline, trimmed to the
   x-overlap, explicit knots — i.e. iteration 04's retained specification plus the trim.
2. For b = 1…B: resample studies with replacement, rebuild subject ids, refit with the
   **original** basis, compute ES (and the extras above). Record convergence.
3. Report the **percentile interval** of the retained replicates, alongside the point
   estimate from step 1, plus how many replicates failed to converge or lacked both
   conditions.
4. Cross-check with **leave-one-study-out**: refit 17 times dropping each study, to see
   whether any single study drives the result. Cheap, and it answers a question a reviewer
   will certainly ask.

**Budget.** The retained model fits in ~17 s, so B = 500 is roughly 2½ hours single-threaded
and well under an hour across cores. Do a B = 20 dry run first to measure the convergence
failure rate before committing.

## What would make this iteration a success

- An interval for the ~10° effect that **does not depend on residual normality** — closing
  the impasse rather than arguing about it.
- An interval whose width reflects **17 studies**, not 3,724 rows. Expect it to be
  substantially wider than any model-based CI so far. That is the point, not a problem.
- A leave-one-study-out panel showing the estimate is not one study's artefact.
- A defensible sentence for the manuscript, of the form: *"in vivo differed from ex vivo by
  X° on average over 15–120° of elevation (study-level bootstrap 95% CI [L, U]); the
  comparison is exploratory because condition is confounded with source study."*

## What this will not fix

The bootstrap widens the interval honestly; it **cannot un-confound the design**. No
resampling scheme recovers a causal contrast when no study measured both conditions. The
conclusion stays descriptive, and iteration 05 should say so as plainly as its predecessors
do.
