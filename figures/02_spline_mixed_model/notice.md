# Notice — how to read the exported model values

The two CSVs in this folder hold **every fitted quantity** of the worked-example model
(scapulothoracic joint · frontal-plane elevation · DoF 2 · 22,102 observations · 44
shoulders). They are written by [`../../analysis_shoulder.R`](../../analysis_shoulder.R);
regenerate with:

```bash
python3 prepare_monolix_data.py
Rscript analysis_shoulder.R
```

The model they describe is

> y<sub>ij</sub> = β₀ + β₁·𝟙(c = in vivo) + f(x<sub>ij</sub>) + g(x<sub>ij</sub>)·𝟙(c = in vivo) + b<sub>i</sub> + ε<sub>ij</sub>

fitted as `bam(Y ~ condO + s(TIME) + s(TIME, by = condO) + s(ID, bs = "re"), method = "fREML", discrete = TRUE, rho = ..., AR.start = ...)`.

There are two files because the model holds two different kinds of number: things that are
**coefficients** (one column of the design matrix each) and things that are not
(smoothing parameters, variances, the AR(1) correlation).

---

## `00_coefficients.csv` — 64 rows, one per coefficient

### The columns

| column | what it is |
| --- | --- |
| `block` | which term of the equation this coefficient belongs to. Five values — see below. This column is not from `mgcv`; it is added by the export script so the file can be grouped and filtered. |
| `term` | the internal name `mgcv` gives the coefficient, i.e. `names(coef(model))`. Useful for matching back to `summary(model)` output. |
| `label` | a readable name. Identical to `term` except for the random-effect rows, where it is replaced by the **shoulder ID** (`article_shoulderid`) so you can tell which specimen each b<sub>i</sub> belongs to. |
| `estimate` | the fitted value of the coefficient — the number the fit actually chose. In the units of y (degrees), for β₀ and the b<sub>i</sub>; for spline coefficients it is a weight on a basis function and is **not** directly a number of degrees. |
| `se` | **standard error** — how well determined that estimate is. |
| `edf` | **effective degrees of freedom** — how much of a free parameter that coefficient really is. |

### `se`, in more detail

The standard error is the square root of the corresponding diagonal element of the model's
Bayesian posterior covariance matrix (`sqrt(diag(model$Vp))`). Read it as the spread of
values that would have been about equally compatible with the data: roughly, the plausible
range is `estimate ± 1.96 × se`.

It is the denominator of the usual test — for the level shift,
`estimate / se = −6.4197 / 1.4891 = −4.311`, which is exactly the *t* value printed in
`00_results_summary.txt`.

> **Do not read the spline coefficients' SEs one at a time.** The nine basis coefficients of
> f are strongly correlated with each other, and an individual B<sub>k</sub> has no
> anatomical meaning on its own. A single β<sub>k</sub> being "non-significant" says nothing
> about whether the curve matters. Judge f and g through the fitted curve and its confidence
> band (`01_population_curves_by_condition.png`,
> `02_difference_invivo_minus_exvivo.png`) and through the smooth-level F-test in
> `00_results_summary.txt`.

### `edf`, in more detail

A plain regression coefficient costs exactly one degree of freedom. A **penalised** one
costs less, because the penalty pulls it toward zero and it is not free to take any value it
likes. The effective degrees of freedom measures what it actually cost — in practice a
number between 0 and 1 per coefficient:

- **edf ≈ 1** — essentially unpenalised, free to be whatever the data says.
- **edf ≈ 0.5** — the penalty is doing half the work of deciding this coefficient.
- **edf ≈ 0** — squeezed out; the coefficient contributes almost nothing.

Summed over a block, it gives the edf reported for that term in `00_results_summary.txt`:

| block | Σ edf | matches the summary line |
| --- | --- | --- |
| `f_spline` | 8.063 | `s(TIME)` |
| `g_difference_spline` | 8.497 | `s(TIME):condOin vivo` |
| `b_i_shoulder_re` | 40.933 | `s(ID)` |

The 44 shoulder offsets together cost only 40.9 degrees of freedom rather than 44 — that
gap is the **shrinkage** the random effect applies, pulling extreme shoulders toward the
population mean.

### The rows

| `block` | rows | symbol | what it is |
| --- | --- | --- | --- |
| `beta0_intercept` | 1 | β₀ | the ex-vivo baseline level, in degrees |
| `beta1_level_shift` | 1 | β₁ | the constant in-vivo level shift (`condO.L`) |
| `f_spline` | 9 | β<sub>k</sub> of f | basis weights for the ex-vivo reference trajectory |
| `g_difference_spline` | 9 | β<sub>k</sub> of g | basis weights for the in-vivo *departure* from f |
| `b_i_shoulder_re` | 44 | b<sub>i</sub> | one intercept offset per shoulder, in degrees |
| | **64** | | = `n_coefficients` in the other file |

**Why 9 and not 10.** The basis is set up with `k = 10`, but a smooth is only identifiable
up to a constant — otherwise it could trade level with the intercept β₀. `mgcv` therefore
applies a sum-to-zero constraint that removes one column, leaving 9. This is also why
`Ref.df` in the summary is not a round number.

**Why 44.** One coefficient per level of `ID`, and there are 44 shoulders in this subset.
Note that one label is `Ludewig et al._` with an empty shoulder id: that study's
`shoulder_id` field is blank in the raw data, so all of its rows collapse into a single
unit.

### Reading one row

```
block             term        label             estimate   se       edf
beta1_level_shift condO.L     condO.L           -6.41969   1.48906  1.000
```

"The in-vivo curve sits about 6.4 units below the ex-vivo one; that offset is determined to
within roughly ±2.9 (1.96 × 1.489); it cost a full degree of freedom because parametric
terms are not penalised."

⚠️ Condition is coded as an **ordered factor**, so `condO.L` is on a polynomial-contrast
scale — read its **sign and t value**, and take the size of the difference in degrees off
the difference curve rather than from this number directly.

```
block           term      label              estimate   se       edf
b_i_shoulder_re s(ID).1   Begon et al._1.0   -7.57606   1.51388  0.942
```

"This particular shoulder sits 7.6° below the population curve; the offset cost 0.94 of a
degree of freedom, i.e. it was shrunk slightly toward zero."

---

## `00_variance_components.csv` — 10 rows

Everything the model estimated that is **not** a coefficient.

### The columns

| column | what it is |
| --- | --- |
| `parameter` | the quantity's name |
| `value` | its fitted value |
| `note` | a one-line reminder of what it is |

### The rows

| parameter | value | how to read it |
| --- | --- | --- |
| `lambda[s(TIME)]` | 0.008656 | **smoothing parameter** of f. Sets how much curvature is penalised: large λ → straighter, small λ → wigglier. Estimated by fREML, not chosen. |
| `lambda[s(TIME):condOin vivo]` | 0.004402 | the same, for g. It is **separate** from f's, so g can be smoothed toward zero while f stays wiggly — which is what lets "the shapes agree" be a reachable result. |
| `lambda[s(ID)]` | 0.038632 | for a random effect the smoothing parameter *is* a variance ratio: λ = σ²/σ<sub>b</sub>². |
| `sigma_b` | 6.3072 | **sd of the shoulder random intercept** — the typical size of a shoulder's offset, in degrees. Here it is comparable in magnitude to the β₁ level shift being tested, which is worth keeping in view. Satisfies `sigma_b = sqrt(scale / lambda[s(ID)])`. |
| `sigma_residual` | 1.2397 | residual sd of ε<sub>ij</sub> — the scatter left after the curves and the shoulder offsets. This is `sqrt(model$sig2)`. |
| `rho` | 0.995554 | the **AR(1) correlation** of successive residuals within a shoulder. It is **fixed**, not estimated jointly: it is read off the lag-one residual autocorrelation of a first fit and passed back in. (The docs round it to 0.996.) |
| `edf_total` | 59.4934 | effective degrees of freedom used by the whole model — the sum of the `edf` column in the other file. Compare to 64 coefficients: the penalties bought back about 4.5. |
| `n_coefficients` | 64 | `length(coef(model))` — the row count of the other file. |
| `n_observations` | 22102 | rows fitted. |
| `BIC` | −31453.55 | model comparison score, lower is better. The same model without the AR(1) correction scores 111,037 and the rejected sigmoid NLME 125,717. |

---

## What is *not* in these files

- **ε<sub>ij</sub>, the residuals.** There is one per observation — 22,102 of them — and they
  are not parameters. What is identified about them is `sigma_residual` and `rho` above.
  Get them with `resid(model)` if you need them.
- **The fitted curves and their confidence bands.** Those are in the PNGs in this folder;
  the coefficients above are the ingredients, not the result.
- **The other 71 joint × motion × DoF combinations.** Those are summarised — one row each,
  no coefficients — in [`../generalized/00_master_summary.csv`](../generalized/00_master_summary.csv).

## One caution about all of it

These are the values of a model whose condition contrast is **perfectly confounded with
source study** — no study measured both in vivo and ex vivo. β₁ and g are therefore
descriptive of the difference between two groups of studies, not causal estimates of an
in-vivo effect. See
[`../../docs/STATISTICAL_METHODS_REVIEW.md`](../../docs/STATISTICAL_METHODS_REVIEW.md).
