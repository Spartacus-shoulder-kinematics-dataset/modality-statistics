# Notice — how to read the exported values

Column-by-column reference for the four CSVs in this folder. For what the iteration *is*
and what it found, see [`README.md`](README.md). For the equivalent document covering the
iteration-02 GAMM, see [`../02_spline_mixed_model/notice.md`](../02_spline_mixed_model/notice.md).

All files are written by [`../../analysis_hlme.R`](../../analysis_hlme.R).

The model behind them:

```r
hlme(fixed = Y ~ ns(TIME, df = 5) * cond, random = ~1, subject = "IDnum", ng = 1)
```

---

## `00_coefficients.csv` — 14 rows, every fitted parameter

| column | what it is |
| --- | --- |
| `parameter` | lcmm's internal name, i.e. `names(m$best)` |
| `estimate` | the fitted value |
| `se` | standard error, `sqrt(diag(VarCov(m)))` |
| `block` | `fixed_effect` (the 12 regression coefficients) or `variance_component` (the 2 variances) |
| `wald` | `estimate / se`, filled in for fixed effects only — see below |

### The rows

| rows | parameter | meaning |
| --- | --- | --- |
| 1 | `intercept` | level of the **ex-vivo** reference curve (ex vivo is the reference level of `cond`) |
| 2–6 | `ns1` … `ns5` | the natural-spline basis weights that shape the **ex-vivo reference curve** |
| 7 | `condin vivo` | the constant **level shift** applied to in-vivo rows |
| 8–12 | `ns1:condin vivo` … `ns5:condin vivo` | the **shape difference** — how the in-vivo curve departs from the ex-vivo one. All zero would mean identical shapes |
| 13 | `varcov 1` | **variance** of the per-shoulder random intercept. Its square root, 6.37, is the typical shoulder offset in degrees |
| 14 | `stderr` | **residual standard deviation** — the sd of ε<sub>ij</sub> |

> `varcov 1` is a **variance** while `stderr` is a **standard deviation**. Take the square
> root of the first before comparing them: √40.58 = 6.37 versus 3.19, i.e. between-shoulder
> scatter is about twice the residual scatter.

### Why 5 spline terms, not 6

`ns(TIME, df = 5)` produces exactly 5 basis columns; the sixth degree of freedom is the
intercept, which the model fits separately. Unlike the GAMM in iteration 02, **nothing is
penalised here** — every coefficient costs a full degree of freedom, which is precisely why
each one can carry an ordinary Wald test.

---

## `00_wald_tests.csv` — the 12 fixed effects, grouped by curve

| column | what it is |
| --- | --- |
| `curve` | which curve the term belongs to — `reference (ex vivo)` or `difference (in vivo - ex vivo)` |
| `term` | the coefficient name |
| `estimate`, `se` | as above |
| `wald` | **`estimate / se`** — signed, so the direction is visible |
| `p_value` | lcmm's two-sided p for that Wald statistic |
| `passed` | `TRUE` when **\|wald\| ≥ 1.96** |

### What the Wald statistic is

The estimate divided by its own standard error: *how many standard errors away from zero is
this coefficient?* Under the usual normal approximation, |Wald| ≥ 1.96 corresponds to a
two-sided p < 0.05, which is the rule used for `passed`.

It answers one question per coefficient — *is this term distinguishable from zero?* — and
that is the point of this iteration: the GAMM of iteration 02 could not answer it per
coefficient, only per smooth.

### The `curve` split, and how the counts are formed

The 12 terms divide into two groups of 6, and each is counted separately:

- **reference (ex vivo)** = `intercept` + `ns1…ns5` — the ex-vivo trajectory itself.
- **difference (in vivo − ex vivo)** = `condin vivo` + `ns1:condin vivo…ns5:condin vivo` —
  everything that switches on only for in-vivo rows.

This is what "how many Wald tests passed per curve" means, and it is the count annotated on
`01_population_curves_by_condition.png` and printed in the digest.

> **Read the counts as descriptive, not as calibrated tests.** The within-shoulder residual
> autocorrelation of 0.970 is unmodelled, so these standard errors are optimistic. See the
> caveat section of [`README.md`](README.md).

> **`intercept` and `condin vivo` are anchors, not quantities.** Every `ns()` basis function
> is exactly zero at the lower boundary knot, so those two coefficients are the ex-vivo level
> and the in-vivo minus ex-vivo difference **at x = −3.46°** — outside the observed range,
> and moving with the knots. `condin vivo` even changes sign between df = 3 and df = 4. Their
> Wald tests are legitimate tests of those particular numbers; they are *not* tests of "is
> there a level difference between conditions". Use the difference curve for that.

---

## `00_wald_tests_by_df.csv` — the same table for every df

36 rows: 6 + 8 + 10 + 12, one block per converged `K`. Same columns as
`00_wald_tests.csv`, with three added at the front:

| column | what it is |
| --- | --- |
| `K` | the `df` of the model this row belongs to |
| `n_fixed` | that model's number of fixed effects, `2K + 2` |
| `retained` | `TRUE` on the rows of the model retained by BIC (here `K = 5`) |

Filtering to `retained == TRUE` reproduces `00_wald_tests.csv` exactly. Use this file to
check whether a conclusion survives the choice of `df` — several do not, which is the point
of keeping it.

---

## `00_model_selection.csv` — one row per natural-spline df tried

| column | what it is |
| --- | --- |
| `K` | the `df` passed to `ns(TIME, df = K)` |
| `n_fixed` | number of fixed effects, always **2K + 2** (intercept, K basis terms, level shift, K interactions) |
| `converged` | `TRUE` when lcmm reported `conv == 1` |
| `loglik` | maximised log-likelihood |
| `AIC` | −2·loglik + 2·npm |
| `BIC` | −2·loglik + npm·log(**number of subjects**) — see the warning below |
| `wald_pass` / `wald_total` | how many of the `n_fixed` coefficients cleared 1.96 |
| `wald_pass_ref` / `wald_total_ref` | the same, counted over the **reference-curve** terms only |
| `wald_pass_diff` / `wald_total_diff` | the same, over the **difference-curve** terms only |
| `sigma_b` | sd of the shoulder random intercept, `sqrt(varcov 1)` — barely moves with `K` |
| `sigma_resid` | residual sd for that `K` |
| `skew`, `kurtosis` | of the residuals; **excess** kurtosis, so 0 is normal |
| `shapiro_p` | Shapiro–Wilk on a random 5,000-row subsample (its hard cap) |
| `seconds` | wall-clock fit time |

> **lcmm's BIC uses the number of SUBJECTS (44), not observations.** `mgcv` uses
> observations. The two iterations' BIC values are therefore **not comparable** — compare
> `07_vs_iteration02_gamm.png` instead.

> **Do not act on `shapiro_p` alone.** At these sample sizes it rejects almost any real
> data; every K here returns p < 1e-17 while showing excess kurtosis of only ~1.1. The Q-Q
> plot in `03_diagnostics.png` is the primary evidence, and it is straight.

---

## `00_thinning.csv` — one row per shoulder

| column | what it is |
| --- | --- |
| `shoulder` | the shoulder ID, `article_shoulderid` |
| `condition` | `ex vivo` or `in vivo` |
| `n_before` | points in the full dataset |
| `n_after` | points actually fitted |
| `thinned` | `TRUE` when `n_after < n_before` |

Shoulders with more than 100 points were reduced to 100, evenly spaced in x-rank, so each
curve's x-coverage and both its endpoints survive intact. 29 of 44 shoulders were thinned;
22,102 rows became 3,724. See the README for why this is defensible and what it changes.

---

## What is *not* in these files

- **The ε<sub>ij</sub> themselves.** One per fitted row (3,724 of them); they are residuals,
  not parameters. What is identified about them is `stderr`. Get them from
  `m$pred$resid_ss` if you need them.
- **The b<sub>i</sub> themselves.** The 44 shoulder intercepts are *predictions*, not free
  parameters — the model estimates only their variance, `varcov 1`. They are plotted in
  `04_random_effects.png` and available as `m$predRE`.
- **The fitted curves.** Those are the PNGs; the coefficients here are the ingredients.
- **The other 71 joint × motion × DoF combinations.** This iteration covers the worked
  example only.
