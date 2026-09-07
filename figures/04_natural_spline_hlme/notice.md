# Notice — how to read the exported values

Column-by-column reference for the seven CSVs in this folder. For what the iteration *is* and
what it found, see [`README.md`](README.md). The equivalent documents for the earlier
iterations are [`../02_spline_mixed_model/notice.md`](../02_spline_mixed_model/notice.md)
and [`../03_natural_spline_hlme/notice.md`](../03_natural_spline_hlme/notice.md).

All files are written by [`../../analysis_random_spline.R`](../../analysis_random_spline.R).
The retained model:

```r
hlme(fixed  = Y ~ ns(TIME, knots=c(40,70,100,130), Boundary.knots=c(0,160)) * cond,
     random = ~ ns(TIME, knots=c(40,70,100,130), Boundary.knots=c(0,160)),
     subject = "IDnum", ng = 1)
```

---

## `00_variants.csv` — the three random-effect specifications

The decisive table of this iteration: one row per specification tried.

| column | what it is |
| --- | --- |
| `variant` | short id — `R1_intercept`, `R2_spline_diag`, `R3_spline_unstr` |
| `spec` | the `random =` argument in words |
| `converged` | `TRUE` when lcmm reported `conv == 1` **and** every parameter is finite |
| `n_param` | total estimated parameters: 12 fixed + the variance block + 1 residual sd |
| `loglik`, `BIC` | maximised log-likelihood and lcmm's BIC (**uses the number of subjects**, 44) |
| `resid_acf1` | **the metric this iteration exists to move** — lag-1 residual autocorrelation *within* shoulder, never crossing a shoulder boundary |
| `resid_sd` | sd of the subject-specific residuals |
| `skew`, `kurtosis` | of those residuals; **excess** kurtosis, so 0 is normal |
| `wald_pass` / `wald_total` | how many fixed effects cleared \|coef/se\| ≥ 1.96 |
| `seconds` | wall-clock fit time |

> **`resid_acf1` alone understates the improvement.** It compares only neighbouring points.
> Iteration 03's autocorrelation was a slow *drift* still near 1 many lags out; here the ACF
> decays to zero within about six lags and then oscillates. Read
> `05_autocorrelation.png`, not just this column.

> **Beware the trade-off in this table.** Moving from `R1` to `R2`/`R3` improves BIC, ACF and
> residual sd, and *worsens* kurtosis from 1.18 to ~13. The two requirements — normal
> residuals and modelled correlation — pull in opposite directions here. See the README.

---

## `00_coefficients.csv` — all 19 parameters of the retained model

| column | what it is |
| --- | --- |
| `parameter` | lcmm's internal name, `names(m$best)` |
| `estimate` | the fitted value |
| `se` | standard error, `sqrt(diag(VarCov(m)))` |
| `block` | `fixed_effect` (12) or `variance_component` (7) |
| `wald` | `estimate / se`, fixed effects only |

### The 12 fixed effects

Identical in meaning to iteration 03 — `intercept`, `ns1…ns5` for the ex-vivo reference
curve, `condin vivo` for the level shift, `ns1:condin vivo…ns5:condin vivo` for the shape
difference.

One difference matters. Because the **boundary knot is now at 0°** rather than at a
data-driven −3.46°, and every `ns()` basis function is zero at the lower boundary knot,
`intercept` and `condin vivo` are the ex-vivo level and the in-vivo minus ex-vivo difference
**at 0° of elevation — the arm at the side**. In iteration 03 those two coefficients were
anchored at an extrapolated point and `condin vivo` even changed sign with `df`. Here they
are interpretable anatomical quantities.

### The 7 variance components

`varcov 1 … varcov 6` are the **variances of the 6 random effects** (intercept + 5 basis
terms) and `stderr` is the residual sd. Because the retained model uses `idiag = TRUE`,
$D$ is diagonal and **no covariances are estimated** — there is one entry per random effect,
in order, and nothing off-diagonal.

For a labelled version with tests, call `VarCovRE(model)` in R.

> Had the unstructured variant been kept, `varcov` would instead hold the 21 lower-triangle
> entries of a full 6×6 covariance, read row by row, with variances at positions 1, 3, 6, 10,
> 15, 21. It was dropped — see the README — so that layout does not arise here.

---

## `00_df_sweep.csv` — the diagonal random spline at every df

One row per `K`. Same idea as iteration 03's `00_model_selection.csv`, but the random part
is `~ns(...)` rather than `~1`.

| column | what it is |
| --- | --- |
| `K` | spline df — the number of `ns()` basis columns |
| `n_random` | random effects, always `K + 1` (intercept + one per basis column) |
| `knots` | the interior knots used, `/`-separated — quantiles of the fitted x rounded to the nearest 10° |
| `converged` | `TRUE` when `conv == 1` and every parameter is finite |
| `n_param` | total estimated parameters: `2K+2` fixed + `K+1` variances + 1 residual sd |
| `BIC` | lcmm's BIC (**number of subjects**, 44, not observations) |
| `resid_acf1` | lag-1 within-shoulder residual autocorrelation |
| `resid_sd` | sd of the subject-specific residuals — **the driver of everything else** |
| `skew`, `kurtosis` | of those residuals; **excess** kurtosis, so 0 is normal |
| `wald_pass` / `wald_total` | fixed effects clearing \|coef/se\| ≥ 1.96 |
| `wald_pass_ref` / `wald_total_ref` | the same, over the **reference-curve** terms only |
| `wald_pass_diff` / `wald_total_diff` | the same, over the **difference-curve** terms only |
| `seconds` | wall-clock fit time |

> **Read `resid_sd` and `kurtosis` together.** They move in lockstep and in opposite
> directions: as `K` rises the random curve absorbs more of each shoulder's trajectory,
> `resid_sd` falls 1.108 → 0.437 and `kurtosis` climbs 5.70 → 12.55, while `resid_acf1`
> improves. There is no `K` satisfying both criteria — see the README.

> **`n_random` is what matters, not `K` alone.** At K = 5 each of the 44 shoulders carries 6
> random effects — 264 predicted quantities from 3,724 rows.

---

## `00_wald_tests.csv` — the 12 fixed effects, grouped by curve

Same schema as iteration 03:

| column | what it is |
| --- | --- |
| `curve` | `reference (ex vivo)` = `intercept` + `ns1…ns5`; `difference (in vivo - ex vivo)` = `condin vivo` + the interactions |
| `term`, `estimate`, `se` | as in the coefficients file |
| `wald` | `estimate / se`, signed |
| `p_value` | lcmm's two-sided p |
| `passed` | `TRUE` when \|wald\| ≥ 1.96 |

> **These standard errors are the honest ones**, and they are *larger* than iteration 03's
> for the difference terms — which is the point. Allowing each shoulder its own curve removes
> the pseudo-replication that made iteration 03's Wald tests optimistic, and the difference
> curve consequently falls from 4/6 to 2/6.
>
> The residuals of this model are **not** normal (excess kurtosis ≈ 13), so the Wald
> approximation rests on an assumption the diagnostics reject. Treat these as descriptive.

---

## `00_wald_tests_by_df.csv` — the same Wald table for every df

36 rows: 6 + 8 + 10 + 12, one block per `K`. Same columns as `00_wald_tests.csv`, with three
added at the front:

| column | what it is |
| --- | --- |
| `K` | the spline df of the model this row belongs to |
| `n_fixed` | that model's fixed effects, `2K + 2` |
| `retained_df` | `TRUE` on the `K = 5` block, the df used for the main figures |

These come from the **diagonal** random-spline sweep, so the `retained_df` rows are *not*
identical to `00_wald_tests.csv`, which is the **unstructured** retained model. The two
differ slightly in standard errors; compare rather than assume.

Use this file to check whether a conclusion survives the choice of `df` — most individual
coefficients do not. Only `ns1:condin vivo` passes at every df, and `condin vivo` fails at
every df.

---

## `00_overlap_trim.csv` — full range vs trimmed to the x-overlap

Two rows, **the same model** — K = 5, interior knots 40/70/100/130, boundary knots 0/160 —
fitted on the full thinned data and on the subset where **both** conditions have data.
Only the data differs; the basis is identical, so the comparison is single-factor.

| column | what it is |
| --- | --- |
| `setting` | `full range` or `trimmed to overlap` |
| `x_lo`, `x_hi` | the elevation range fitted |
| `n_rows` | rows in that fit |
| `resid_sd`, `kurtosis`, `resid_acf1` | the residual diagnostics to compare |

> **BIC is deliberately absent.** The two rows are *different datasets*, so their likelihoods
> are not comparable and a BIC column would invite exactly the wrong comparison. Only the
> residual diagnostics can be read across the rows.

---

## `00_thinning.csv`

Not written by this script — the thinning is identical to iteration 03, so see
[`../03_natural_spline_hlme/00_thinning.csv`](../03_natural_spline_hlme/00_thinning.csv)
and its notice. 22,102 rows become 3,724 at ≤100 points per shoulder.

## What is *not* in these files

- **The ε<sub>ij</sub>** — one per fitted row; residuals, not parameters. `m$pred$resid_ss`.
- **The b<sub>i</sub>** — 6 per shoulder, *predictions* rather than free parameters; the
  model estimates their covariance. Plotted in `04_random_effects.png`, available as
  `m$predRE`.
- **The other 71 joint × motion × DoF combinations** — worked example only.
