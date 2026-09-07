# Iteration 03 — natural spline + `lcmm::hlme`

Produced by [`../../analysis_hlme.R`](../../analysis_hlme.R). Same subset and the same
figure set as [iteration 02](../02_spline_mixed_model/), but a different model, chosen so
that **every parameter can be tested individually**.

```bash
python3 prepare_monolix_data.py
Rscript analysis_hlme.R          # ~1 min 35 s
```

For what each CSV column means, see [`notice.md`](notice.md).

---

## Why this iteration exists

Iteration 02's GAMM works, but it spends **64 coefficients** on this one cell and none of
them can be read on its own: the penalised thin-plate basis weights are correlated, have no
anatomical meaning, and are only judged through a whole-smooth F-test. This iteration asks
the opposite question — *how few parameters can we get away with, and can we test each one?*

| | iteration 02 | iteration 03 |
| --- | --- | --- |
| basis | penalised thin-plate, k = 10 | natural spline, `ns(TIME, df = 5)` |
| fitter | `mgcv::bam`, fREML | `lcmm::hlme`, full ML |
| parameters | 64 coefficients + 3 λ + 2 variances | **12 fixed effects + 2 variances** |
| inference | F-test per smooth | **Wald test per coefficient** |
| correlation | AR(1), ρ fixed at 0.996 | none — see the caveat below |
| residual normality | not checked | **checked; it is the acceptance criterion** |
| data | all 22,102 rows | thinned to 3,724 |
| runtime | ~4 s | ~95 s |

## The model

```r
hlme(fixed   = Y ~ ns(TIME, df = 5) * cond,
     random  = ~ 1,
     subject = "IDnum",
     ng      = 1)
```

For shoulder $i$, observation $j$, measured in condition $c_i \in \{\text{ex vivo},\ \text{in vivo}\}$,
at thoracohumeral elevation $x_{ij}$:

$$
\begin{aligned}
y_{ij}
= \underbrace{\beta_0 + \sum_{k=1}^{5}\beta_k B_k(x_{ij})}_{\text{ex-vivo reference curve}}
+ \underbrace{\left[\gamma_0 + \sum_{k=1}^{5}\gamma_k B_k(x_{ij})\right]\mathbb{1}(c_i = \text{in vivo})}_{\text{in-vivo departure}}
+ b_i + \varepsilon_{ij}
\qquad b_i \sim \mathcal{N}(0,\sigma_b^{2}),
\quad \varepsilon_{ij} \stackrel{\text{iid}}{\sim} \mathcal{N}(0,\sigma_\varepsilon^{2}).
\end{aligned}
$$

$\mathbb{1}(c_i = \text{in vivo})$ is 1 on in-vivo rows and 0 on ex-vivo rows, so the two
conditions are two different curves sharing one equation:

$$
\text{ex vivo:}\quad y_{ij} = \beta_0 + \sum_{k=1}^{5}\beta_k B_k(x_{ij}) + b_i + \varepsilon_{ij}
$$

$$
\text{in vivo:}\quad y_{ij} = (\beta_0+\gamma_0) + \sum_{k=1}^{5}(\beta_k+\gamma_k) B_k(x_{ij}) + b_i + \varepsilon_{ij}
$$

So $\gamma_0$ shifts the curve vertically and $\gamma_1,\dots,\gamma_5$ bend it. All six
$\gamma$ equal to zero would mean the conditions are indistinguishable — which is exactly
what the Wald tests below are testing, one coefficient at a time.

### What the $B_k$ are

$B_1,\dots,B_5$ are the **natural cubic spline basis** returned by `splines::ns(TIME, df = 5)`.
They are fixed, known functions of $x$ — computed before any fitting, and identical for both
conditions. Only the twelve $\beta$ and $\gamma$ are estimated.

`df = 5` places **4 interior knots** at the 20/40/60/80th percentiles of the fitted $x$,
with the boundary knots at its range:

$$
\xi = (38.15,\ 71.40,\ 102.69,\ 126.85), \qquad
[\xi_{\min}, \xi_{\max}] = [-3.46,\ 187.02]
$$

The spline is a cubic polynomial between consecutive knots, joined so that the function and
its first two derivatives are continuous. **"Natural"** adds the boundary constraint that
makes it usable here: $S''(x) = S'''(x) = 0$ outside $[\xi_{\min}, \xi_{\max}]$, i.e. the
curve is forced to be **linear beyond the outermost knots** instead of letting a cubic fly
off. That constraint is what buys back the 2 degrees of freedom a free cubic spline on 4
knots would need, so 5 basis columns suffice.

The space these five functions span can be written explicitly in truncated-power form, with
$m = 6$ knots (4 interior plus the 2 boundaries) and $(u)_+ = \max(u,0)$:

$$
N_1(x) = 1,\qquad N_2(x) = x,\qquad N_{k+2}(x) = d_k(x) - d_{m-1}(x),
$$

$$
\text{where}\quad d_k(x) = \frac{(x-\xi_k)^3_+ - (x-\xi_m)^3_+}{\xi_m - \xi_k}.
$$

> `ns()` does **not** return these columns. It returns an equivalent, far better-conditioned
> B-spline parametrisation of the *same* function space, so the fitted $\beta_k$ are not
> truncated-power coefficients and should not be read as slopes or curvatures. Any individual
> $B_k$ is meaningless on its own; the curve is the sum. What each $\beta_k$ *does* support —
> unlike the penalised basis of iteration 02 — is an honest one-degree-of-freedom Wald test.

### The twelve fixed effects, fitted

| symbol | term in `00_coefficients.csv` | estimate | se | Wald |
| --- | --- | ---: | ---: | ---: |
| $\beta_0$ | `intercept` | 4.084 | 2.861 | 1.43 |
| $\beta_1$ | `ns1` | −10.480 | 1.995 | −5.25 |
| $\beta_2$ | `ns2` | −22.911 | 2.456 | −9.33 |
| $\beta_3$ | `ns3` | −52.024 | 1.367 | **−38.05** |
| $\beta_4$ | `ns4` | −53.581 | 5.954 | −9.00 |
| $\beta_5$ | `ns5` | −30.175 | 5.765 | −5.23 |
| $\gamma_0$ | `condin vivo` | −6.064 | 3.123 | −1.94 |
| $\gamma_1$ | `ns1:condin vivo` | −13.716 | 2.049 | −6.69 |
| $\gamma_2$ | `ns2:condin vivo` | −9.518 | 2.540 | −3.75 |
| $\gamma_3$ | `ns3:condin vivo` | 8.533 | 1.427 | 5.98 |
| $\gamma_4$ | `ns4:condin vivo` | −5.543 | 6.087 | −0.91 |
| $\gamma_5$ | `ns5:condin vivo` | −18.342 | 5.807 | −3.16 |

plus two variance components: $\sigma_b^2 = 40.581$ (so $\sigma_b = 6.370$) and
$\sigma_\varepsilon = 3.188$.

`df` was selected over $K \in \{2,3,4,5\}$ by BIC — all four converged, all four gave
near-normal residuals, and $K = 5$ won (see `00_model_selection.csv` and figure 06).

`08_by_df_population_curves.png` shows all four side by side, one panel per df on its own axes,
with the interior knots marked as rug ticks and the retained model boxed:

| df | params | BIC | Wald |
| ---: | ---: | ---: | :--- |
| 2 | 6 | 19713 | 6/6 |
| 3 | 8 | 19674 | 8/8 |
| 4 | 10 | 19517 | 6/10 |
| **5** | **12** | **19508** | **9/12** |

Two things are worth reading off that figure. First, **df = 2 and 3 are visibly too rigid** —
both conditions come out as near-straight lines, and although *every* coefficient passes its
Wald test (6/6, 8/8), that is a symptom of under-fitting rather than a good sign: with too
few basis functions each one is forced to carry a lot of signal. Second, the extra degrees of
freedom buy shape almost entirely for the **ex-vivo** curve, which only develops its
characteristic S-bend at df = 4–5; the in-vivo curve stays close to linear throughout. That
asymmetry is consistent with the ex-vivo data being far denser per shoulder.

Note also that BIC barely separates df = 4 from df = 5 (19517 vs 19508, a difference of 9),
so the retained model is not a decisive win — df = 4 would be a defensible alternative with
two fewer parameters.

`09_by_df_difference.png` repeats the comparison for the **difference curve**
(in vivo − ex vivo). These four panels deliberately **share one y-axis**, because the point
is to compare the size of the difference across df, which separate scales would hide. Marks
along the bottom show where the 95% band excludes zero.

It carries the clearest warning in this iteration:

| df | mean \|difference\| | behaviour above ~130° |
| ---: | ---: | :--- |
| 2 | 10.05° | **crosses zero, rises to +7°** |
| 3 | 9.88° | **crosses zero, rises to +5°** |
| 4 | 10.43° | stays negative, flattens to ≈ −4° |
| **5** | **10.48°** | stays negative, flattens to ≈ −6° |

**df = 2 and 3 predict a sign reversal at high elevation that df = 4 and 5 do not.** Under
the rigid models the in-vivo curve ends up *above* ex-vivo past ~135°; under the flexible
ones the difference simply shrinks toward zero without ever changing sign. That is a
different biomechanical claim, produced entirely by how many basis functions the spline was
allowed — and it is exactly the kind of artefact a near-parabola generates when forced onto
an S-shaped trajectory. It is also the region where the ex-vivo data thins out, so neither
version is well supported there.

The reassuring half: the **mean absolute difference is stable at 9.9–10.5° across every
df**. The overall effect size does not depend on the basis; only its shape at the extremes
does. Read the ~10° as the robust number and treat the end-range behaviour as unresolved.

`10_by_df_diagnostics.png` and `11_by_df_random_effects.png` complete the sweep with the two
distributional assumptions, checked at every df rather than only for the winner.

`10_by_df_diagnostics.png` gives each df a row: the normal Q-Q on the left, the histogram
against a fitted normal on the right — the same pair `03_diagnostics.png` shows for the
retained model.

**By the summary statistics, residual normality does not depend on the df at all.** Skew
stays in −0.37…−0.25 and excess kurtosis in 1.03…1.19, and all four Q-Q plots are straight
through the body with the same mild S at the tails (sample quantiles slightly inside the
line at both ends, i.e. tails a little *shorter* than normal). Nothing threatens the
acceptance criterion at any df.

**But the histograms disagree with the moments, and they are the more useful diagnostic
here.** At df = 2 and 3 the residual distribution is visibly lumpy — a dip at zero with mass
piled to either side, most obvious at df = 3 — while at df = 4 and 5 it settles into a
single smooth mode. That is unmodelled curvature leaking into the residuals: the rigid
splines cannot follow the trajectory, so points sit systematically above the fit in some
elevation bands and below it in others. Skew and kurtosis are almost blind to it, because a
symmetric two-humped distribution has roughly the moments of a normal. It is a good argument
for reading the graph rather than the two numbers — and an independent reason to prefer
df = 4–5, quite apart from BIC.

**σ<sub>b</sub> is essentially invariant across df** — 6.368, 6.385, 6.349, 6.370 at
df = 2, 3, 4, 5 (now exported as the `sigma_b` column of `00_model_selection.csv`, along
with `sigma_resid`). Adding basis functions changes the *population curve*, not the spread
between shoulders, which is what you would hope and another sign the between-shoulder term
is well identified — compare the mgcv GAMM's 6.31. The b<sub>i</sub> Q-Q plots in
`11_by_df_random_effects.png` are close to straight with a visible flattening in the upper
tail at every df: the two or three highest shoulders sit lower than a normal would put them,
so normality on b<sub>i</sub> is adequate but not perfect. With n = 44 that is not worth
chasing.

### Manuscript form

Ready to paste into `article/`:

```latex
\begin{equation}
y_{ij} = \beta_{0} + \sum_{k=1}^{5}\beta_{k}B_{k}(x_{ij})
       + \left[\gamma_{0} + \sum_{k=1}^{5}\gamma_{k}B_{k}(x_{ij})\right]
         \mathbb{1}(c_{i}=\text{in vivo})
       + b_{i} + \varepsilon_{ij},
\label{eq:hlme}
\end{equation}
```

where $B_{1},\dots,B_{5}$ is a natural cubic spline basis with interior knots at the
quintiles of the thoracohumeral angle, $b_{i}\sim\mathcal{N}(0,\sigma_b^{2})$ is a
shoulder-specific random intercept and
$\varepsilon_{ij}\sim\mathcal{N}(0,\sigma_\varepsilon^{2})$.

## Thinning

`hlme` costs roughly **cubically** in points per shoulder, so the dense ex-vivo curves were
capped at 100 points each, evenly spaced in x-rank (both endpoints always kept, so each
curve's x-coverage is preserved exactly).

```
22,102 rows -> 3,724 rows   (29 of 44 shoulders thinned)
ex vivo: 18,302 -> 1,101    in vivo: 3,800 -> 2,623
```

This is defensible because at ρ ≈ 0.996 a 2,987-point curve carries only a handful of
independent points — but it *is* a subsample, and `00_thinning.csv` records exactly what
was dropped. Note it also reverses the row imbalance: ex vivo has 13 shoulders and in vivo
31, so capping leaves in-vivo rows in the majority. The statistical unit is the shoulder,
so this does not bias the fit, but it is worth knowing when reading the raw scatter.

---

## Results

**Wald tests** — `|coef/se| ≥ 1.96`, in `00_wald_tests.csv`:

| curve | passed | which ones fail |
| --- | --- | --- |
| reference (ex vivo) | **5 / 6** | the intercept (1.43) |
| difference (in vivo − ex vivo) | **4 / 6** | the level shift (1.94, just under) and `ns4:cond` (0.91) |

Read that as: the ex-vivo trajectory is very well determined (Wald up to 38 on `ns3`), and
the two conditions clearly differ in **shape** — three of the five interaction terms pass,
one at Wald 6.7. What is *not* resolved is the constant level shift, which lands at 1.94,
fractionally below threshold. So the evidence is for a difference in the *form* of the
rhythm rather than a clean vertical offset.

**The same table for every df** is in `00_wald_tests_by_df.csv` (36 rows: 6 + 8 + 10 + 12,
with `K` and `retained` columns), and the per-curve counts are now columns of
`00_model_selection.csv`:

| df | reference curve | difference curve | total |
| ---: | :--- | :--- | :--- |
| 2 | 3/3 | 3/3 | 6/6 |
| 3 | 4/4 | 4/4 | 8/8 |
| 4 | 4/5 | **2/5** | 6/10 |
| **5** | **5/6** | **4/6** | **9/12** |

Note df = 4 resolves *fewer* difference terms than df = 5 (2/5 against 4/6) despite being
the simpler model — pass counts are not monotone in df, so "more tests passed" is not a
model-selection criterion.

### Do not read $\gamma_0$ as "the level shift"

The per-df table exposes a trap. $\gamma_0$ (`condin vivo`) **changes sign** with the basis:

| df | $\gamma_0$ | Wald | passes |
| ---: | ---: | ---: | :--- |
| 2 | **+5.913** | +2.69 | yes |
| 3 | **+8.163** | +3.48 | yes |
| 4 | −3.743 | −1.44 | no |
| 5 | −6.064 | −1.94 | no |

That is not instability in the data — it is what $\gamma_0$ *means*. All five `ns()` basis
functions are **exactly zero at the lower boundary knot** (verified: `predict(B, -3.46)`
returns zeros at every df), so $\gamma_0$ is the in-vivo minus ex-vivo difference **at
x = −3.46°** — a point outside the observed range, reached by extrapolation, whose location
and surroundings move as the knots move.

So $\gamma_0$ is an anchor for the parametrisation, not a biomechanical quantity. It is
*not* the mean level shift, and neither its sign nor its Wald test should be reported as
"the difference between conditions". The interpretable quantity is the **difference curve**
$\gamma_0 + \sum_k \gamma_k B_k(x)$ evaluated across the observed range — figures 02 and 09
— whose mean magnitude is a stable ~10° at every df. The same caution applies to the
intercept $\beta_0$, which fails its Wald test for exactly the same reason.

**Residual normality** — the acceptance criterion, and it is met:

```
skewness  -0.25      excess kurtosis  1.18      within +/-1.96 sd  94.2%  (normal 95.0%)
```

The Q-Q plot in `03_diagnostics.png` is straight and the histogram sits on the normal curve.
Kurtosis stays in 1.03–1.19 across every K, so this is a stable property of the model, not
an artefact of the retained df. (Shapiro–Wilk returns p ≈ 3e-18, but at n = 3,724 it rejects
almost any real data — the Q-Q plot is the evidence.)

**Cross-check against iteration 02:** the between-shoulder sd comes out at **6.37**, against
**6.31** from the mgcv GAMM. Two different bases, two different fitters, two different row
counts, essentially the same answer — see `07_vs_iteration02_gamm.png`.

---

## The caveat that matters

**Within-shoulder residual lag-1 autocorrelation is 0.970, and it is not modelled.**

`04_random_effects.png` and especially the *residuals vs x* panel of `03_diagnostics.png`
show why: each shoulder's residuals trace a smooth arc. A shoulder's deviation from the
population curve is **itself a smooth function of elevation**, not noise, and a random
intercept absorbs only its *level*, not its *shape*.

Consequence: the Wald standard errors are **optimistic**. Treat the counts above as a
descriptive summary of which terms are well determined, not as calibrated significance
tests. Iteration 02, which corrects for this with a fixed ρ, remains the reference for
inference.

### Why `cor = AR(TIME)` does not fix it

It was tried, and it is documented in `00_results_summary.txt` rather than deleted:

- The estimated decay was `cor1 = 0.000566`, so the correlation `exp(-cor1·Δx) ≈ 1` even
  100° apart. **A correlation that never decays is a random intercept**, so the AR process
  is aliased with `random = ~1`.
- It won that competition and collapsed the others: between-shoulder variance
  40.58 → 0.0007, residual sd 3.19 → 0.10.
- Its reported `resid_ss` then has kurtosis 175 — an artefact, because `resid_ss` subtracts
  the *predicted AR path*, which tracks the data almost exactly. It is not the ε<sub>ij</sub>
  whose normality we require.

Fitting `AR(TIME)` with no random intercept gives an identical fit (`cor1 = 0.0006`,
`cor2 = 7.139`), confirming the aliasing. `BM(TIME)` diverges outright. Random slopes on the
first two spline terms only moved the residual ACF from 0.977 to 0.961 while worsening
kurtosis to 7.2.

The honest fix is a **random smooth per shoulder** with enough flexibility to absorb shape
deviation — outstanding item 4 in
[`../../docs/STATISTICAL_METHODS_RESPONSE.md`](../../docs/STATISTICAL_METHODS_RESPONSE.md).

### And the standing caveat

Condition is perfectly confounded with source study — no study measured both in vivo and ex
vivo. Everything here is **descriptive, not causal**.

---

## The files

| file | what it is |
| --- | --- |
| `00_results_summary.txt` | human-readable digest of everything below |
| `00_model_selection.csv` | the K grid: BIC, Wald counts, normality, runtime per df |
| `00_coefficients.csv` | all 14 parameters with SE and Wald |
| `00_wald_tests.csv` | the retained model's 12 fixed effects, grouped by curve, with pass/fail |
| `00_wald_tests_by_df.csv` | the same, for **every** df — 36 rows with `K` and `retained` columns |
| `00_thinning.csv` | per shoulder, points before and after thinning |
| `01_population_curves_by_condition.png` | fitted curves with 95% bands, Wald counts annotated |
| `02_difference_invivo_minus_exvivo.png` | in-vivo − ex-vivo difference on the x-overlap only |
| `03_diagnostics.png` | **the normality evidence** — Q-Q, histogram, residuals vs fitted and vs x |
| `04_random_effects.png` | the 44 shoulder intercepts, histogram and Q-Q |
| `05_autocorrelation.png` | the unmodelled residual ACF — the caveat, drawn |
| `06_model_comparison.png` | BIC and Wald pass-rate across K |
| `07_vs_iteration02_gamm.png` | this model overlaid on the iteration-02 GAMM |

The four `by_df` figures repeat the retained-model views across every df, so each is named
after the single-model figure it mirrors:

| retained model (df = 5) | same view, all four df |
| --- | --- |
| `01_population_curves_by_condition.png` | `08_by_df_population_curves.png` |
| `02_difference_invivo_minus_exvivo.png` | `09_by_df_difference.png` |
| `03_diagnostics.png` | `10_by_df_diagnostics.png` |
| `04_random_effects.png` | `11_by_df_random_effects.png` |

| `08_by_df_population_curves.png` | **2×2: the fitted curves at df = 2, 3, 4, 5**, each on its own axes, knots rugged, retained model boxed |
| `09_by_df_difference.png` | **2×2: the in-vivo − ex-vivo difference at each df**, shared y-axis, with the elevations where the band excludes zero marked |
| `10_by_df_diagnostics.png` | **one row per df: normal Q-Q + histogram against a fitted normal**, the same pair `03_diagnostics.png` shows for the retained model |
| `11_by_df_random_effects.png` | **2×2: normal Q-Q of the 44 shoulder intercepts at each df**, shared axes, σ<sub>b</sub> annotated |
