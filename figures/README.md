# Figures

Two things live here:

1. **The worked example** (`00_data_exploration/`, `01_sigmoid_nlme/`,
   `02_spline_mixed_model/`, `03_natural_spline_hlme/`, `04_natural_spline_hlme/`, `05_bootstrap/`) — the deep,
   pedagogical single-subset analysis, kept as iteration sets (see below).
2. **`generalized/`** — the same model applied to **every** joint x movement x
   DoF, with one **planche per movement**. See
   [`generalized/00_SUMMARY.md`](generalized/00_SUMMARY.md).

---

## Part A — the worked example (iteration sets)

These figures are produced by [`../analysis_shoulder.R`](../analysis_shoulder.R)
for the feature-example subset (scapulothoracic, frontal plane elevation,
DoF = 2, rad). We deliberately **keep every iteration** so the reasoning — *why
we moved from one model to the next* — stays visible. Read the folders in order.

> NOTE: `analysis_shoulder.R` regenerates the PNG folders but **not** this file.
> If you `rm -rf figures`, keep (or restore) this README separately.

Regenerate the figures with:

```bash
python3 prepare_monolix_data.py   # writes monolix_st_frontal_dof2.csv
Rscript analysis_shoulder.R       # writes every figure below
```

---

## `00_data_exploration/` — look before you model

Model-agnostic views of the data. They establish the three facts that drive
every later decision:

| figure | what it shows | consequence |
| --- | --- | --- |
| `01_raw_scatter.png` | all 22,102 points, coloured by condition | y rises monotonically with elevation; **no plateau** |
| `02_spaghetti_by_unit.png` | one curve per `article × shoulder` unit | big **between-unit** differences ⇒ need random effects |
| `03_x_sampling_density.png` | where x is sampled, per condition | sampling is **irregular / study-specific** ⇒ wide-format LGCM is inappropriate |
| `04_points_per_unit.png` | points per unit (sorted) | **11 to 2,987** points/unit — extreme imbalance |
| `05_binned_mean_sd.png` | binned mean ± SD per condition | first look at the average shapes |

---

## `01_sigmoid_nlme/` — attempt 1: parametric sigmoid (why we moved on)

Following METHODOLOGY.md's first suggestion, we tried a **4-parameter logistic
(sigmoid) nonlinear mixed-effects model** — the "Monolix-style" parametric
approach.

**Why we moved on:**

1. **The data never saturates.** Over the observed elevation range the ST angle
   keeps falling; it does not flatten to an upper/lower plateau. A logistic's
   asymptotes therefore sit *outside* the data (fitted `xmid` lands at the data
   edge, `B` far below any observation). See `01_sigmoid_pooled_fit.png`.
2. **Random effects were non-identifiable.** Putting random effects on all four
   logistic parameters (with only 44 units, each covering a different, partial
   slice of the curve) made the optimiser blow up: random-effect SDs exploded to
   ~5e4–5e5 and *every* fixed-effect p-value collapsed to ≈ 1 — a degenerate
   model, not a result. The only stable version keeps a **random intercept
   only**, and even then it **fits worst** (BIC ≈ 125,700).
3. **It misses the shape.** `02_parametric_vs_spline.png` overlays the sigmoid
   (dashed) on the spline (solid): the sigmoid is forced to plateau at both ends
   while the real curves keep descending past ~130°.

Lesson: **match the structural model to the data.** A sigmoid is the wrong shape
here, so we switched to a flexible spline.

> Note for the Monolix path: the same lesson applies there. If you use
> MonolixSuite later, prefer a structural model that does **not** impose
> saturation (e.g. a low-order polynomial / spline-like form), or expect the
> same identifiability trouble with a sigmoid.

---

## `02_spline_mixed_model/` — the fix: penalised spline mixed model + AR(1)

The model the data actually supports — still a (penalised) **nonlinear mixed
model**, fit with `mgcv`:

```
Y ~ condO + s(TIME) + s(TIME, by = condO) + s(ID, bs = "re")
```

- `s(TIME)` — the ex-vivo (reference) nonlinear trajectory.
- `s(TIME, by = condO)` — the **in-vivo difference** smooth (the hypothesis).
- `s(ID, bs = "re")` — a per-unit random intercept (the mixed part).

**One more correction that matters:** each unit is a dense, smooth curve, so
residuals are strongly autocorrelated (**lag-1 ≈ 0.996**, see
`05_autocorrelation.png`). Ignoring this makes every p-value look absurdly
significant. We refit with an **AR(1)** residual model and report the corrected
inference.

| figure | what it shows |
| --- | --- |
| `01_population_curves_by_condition.png` | **headline** — fitted population trajectory per condition (± 95% CI) over the raw data |
| `02_difference_invivo_minus_exvivo.png` | **the test** — in-vivo minus ex-vivo difference across elevation, with 95% CI; where the band excludes 0, the conditions differ |
| `03_diagnostics.png` | residuals vs fitted, and normal Q-Q |
| `04_random_effects.png` | spread of the per-unit random intercepts |
| `05_autocorrelation.png` | residual autocorrelation (why AR(1) is needed) |
| `06_model_comparison.png` | BIC: spline+AR1 vs spline-naive vs sigmoid NLME |
| `00_results_summary.txt` | the numeric digest + how to read it |

**Result (AR(1)-corrected, honest inference).** The in-vivo and ex-vivo
trajectories differ significantly — both a constant **level shift**
(`condO.L`: t = −4.3, p < 0.001) and a **shape difference**
(`s(TIME):condOin vivo`: p < 0.001). The difference figure shows in-vivo sitting
~10–15° below ex-vivo through the low-to-mid elevation range, converging toward
zero near end-range.

**Caveats** (see `00_results_summary.txt`): values labelled `rad` look like
degrees (verify provenance); only 44 heavily imbalanced units, so between-study
variance is modest to estimate; AR(1) handles within-curve autocorrelation only
approximately.

---

## `03_natural_spline_hlme/` — how few parameters can we get away with?

Produced by [`../analysis_hlme.R`](../analysis_hlme.R), not `analysis_shoulder.R`.
It has its own [`README.md`](03_natural_spline_hlme/README.md) and
[`notice.md`](03_natural_spline_hlme/notice.md).

Iteration 02 works but spends **64 coefficients**, none of which can be read
individually: the penalised basis weights are correlated, meaningless on their own,
and only testable as a whole smooth. So iteration 03 asks the opposite question —
replace the penalised thin-plate basis with a plain **natural spline**
(`splines::ns`), fit by **`lcmm::hlme`**, and test **every coefficient** with a Wald
statistic.

```bash
Rscript analysis_hlme.R           # ~1 min 35 s
```

| | iteration 02 | iteration 03 |
| --- | --- | --- |
| parameters | 64 coefficients + 3 λ + 2 variances | **12 fixed effects + 2 variances** |
| inference | F-test per smooth | **Wald test per coefficient** |
| autocorrelation | AR(1), ρ fixed at 0.996 | **not modelled** — see below |
| residual normality | not checked | **checked — the acceptance criterion** |

**Result.** `ns(TIME, df=5)` retained by BIC. Wald tests (`|coef/se| ≥ 1.96`) pass
**5/6** on the ex-vivo reference curve and **4/6** on the difference curve — the
shape difference is well determined (three interaction terms pass, one at Wald 6.7),
while the constant level shift lands at 1.94, just under threshold. Residuals are
near-normal (skew −0.25, excess kurtosis 1.18, straight Q-Q). The between-shoulder
sd comes out at **6.37** against the GAMM's **6.31** — an independent confirmation
that the two iterations agree on the parts they share.

**The caveat.** Within-shoulder residual autocorrelation is **0.970 and unmodelled**,
so the Wald SEs are optimistic — the counts are descriptive, not calibrated tests.
`cor = AR(TIME)` cannot fix it here: its decay estimates to ~0, which makes it
*aliased with the random intercept*, and it then swallows both the random-effect and
residual variances. That failed fit is kept and documented in
`00_results_summary.txt`, the same way iteration 01's sigmoid is kept.

**Iteration 02 remains the reference for inference.** Iteration 03 is the answer to
"can each parameter be tested", not a replacement.

---

## `04_natural_spline_hlme/` — a random CURVE per shoulder

Produced by [`../analysis_random_spline.R`](../analysis_random_spline.R). Its own
[`README.md`](04_natural_spline_hlme/README.md) and
[`notice.md`](04_natural_spline_hlme/notice.md).

Iteration 03 left its residual autocorrelation at **0.978, unmodelled**, and said the honest
fix was a random *smooth* per shoulder. This iteration does that — the canonical `lcmm`
pattern, putting the whole spline basis in the random effects so each shoulder gets its own
curve rather than just its own height:

```r
hlme(fixed = Y ~ ns(TIME, knots) * cond, random = ~ ns(TIME, knots), ...)
```

It also switches to **explicit knots** (interior 40/70/100/130°, boundary 0/160°), which
makes the intercept and level-shift terms interpretable at 0° of elevation instead of at an
extrapolated point.

**It fixes 03's defect and breaks 03's guarantee.** BIC falls 19,510 → 6,348, residual sd
3.17° → 0.43°, and the ACF changes from a slow drift to a short-range oscillation that dies
within ~6 lags. But excess kurtosis goes **1.18 → 12.98**: once the random curve absorbs each
shoulder's trajectory, what is left is crumbs, and the normality criterion fails.

The two requirements pull against each other:

| | autocorrelation handled | residuals normal |
| --- | :---: | :---: |
| iteration 03 (`random = ~1`) | ✗ | ✓ |
| iteration 04 (`random = ~ns`) | mostly ✓ | ✗ |

The scientific consequence is real: with shoulders allowed their own shapes, the difference
curve drops from **4/6 to 2/6** passing Wald terms — iteration 03's standard errors were
optimistic, as it had warned. What survives in both is a shape difference at low-to-mid
elevation; the level difference at 0° is −4.06° ± 2.31 and not distinguishable from zero.

Two follow-ups came out of it. **The unstructured covariance was dropped** — 40× the runtime
(659 s vs 17 s) for 101 BIC units and an identical ACF — so the retained model has 19
parameters with a diagonal *D*. And **trimming to the x-overlap is the best single lever on
the kurtosis found anywhere**: restricting to the range where both conditions have data
(here [14°, 150°], computed from the data so it carries to the other 71 cells) keeps 90.7% of
rows and takes excess kurtosis from **12.55 to 5.25** with the autocorrelation unchanged —
at an unchanged model (K = 5, same knots), so only the data differs. Narrowing the
*boundary knots* instead makes it worse (12.55 → 14.84): that changes the basis, not which
points are fitted. It halves the problem; it does not solve it.

---

## `05_bootstrap/` — planned, not implemented

A design note only — no script, no outputs. See
[`05_bootstrap/README.md`](05_bootstrap/README.md).

Iterations 03 and 04 reach an impasse that no further modelling choice resolves: model-based
standard errors need assumptions the diagnostics reject, and every test so far is
pseudo-replicated because condition is confounded with source study. A **study-level
bootstrap** — resampling the 17 studies, not the 44 shoulders and never the rows — fixes both
with one mechanism, and gives the stable ~10° effect size the interval it currently lacks.

---

## Part B — `generalized/` (every joint x movement x DoF)

Produced by [`../analyze_all.R`](../analyze_all.R), which runs the same spline
mixed model (+AR1) over **all 72** angular combinations and assembles the results
as **planches** — one plate per movement, in the shoulder-kinematics house style
(`shoulder-kinematics/spartacus/plots`):

```bash
python3 prepare_monolix_data.py   # writes spartacus_angles_long.csv
Rscript analyze_all.R             # writes everything under generalized/
```

- **`planche_<movement>.png`** — a 4 x 3 plate: **rows = joints**
  (sternoclavicular → acromioclavicular → scapulothoracic → glenohumeral,
  proximal → distal), **cols = DoF 1/2/3** (anatomical angle names as panel
  titles). Each panel overlays the in-vivo (green) / ex-vivo (orange) spline
  fits (±95% CI) on the raw study curves. **All three panels in a row share one
  y-axis** so the rotations are comparable within a joint. Each panel is
  annotated:
  - **top-right stars** — the in-vivo vs ex-vivo **trajectory-difference** test
    (`*** p<0.001, ** p<0.01, * p<0.05, ns`). Note these are almost always
    significant given the sample sizes — read the effect size, below.
  - **bottom-left text** — n shoulders (ex/in), mean & max `|in-vivo − ex-vivo|`
    in degrees, and the **level shift** (deg) with its own significance stars
    (the level-shift test is the discriminating one).
  - Cells with only one condition show a single pooled **black** fit labelled
    "single group (…)"; cells with too few units (<4 shoulders) show raw data
    only.
- **`00_master_summary.csv`** — one row per combo: mode, effect sizes, p-values.
- **`00_SUMMARY.md`** — index + the compared combos sorted by **effect size**
  (`mean |in-vivo − ex-vivo|`, degrees, measured only on the overlapping x-range).
- **`<joint>/<movement>/dof<k>/`** — per-combo drill-down curves.

Of the 72 combos, **33** have both conditions (in-vivo vs ex-vivo tested), 30
have a single condition, 9 are skipped (too few units). Read **magnitude, not
just p**: with these sample sizes the difference is almost always "significant",
so the effect size in degrees is the meaningful quantity.
