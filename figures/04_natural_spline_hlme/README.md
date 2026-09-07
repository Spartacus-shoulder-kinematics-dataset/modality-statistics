# Iteration 04 — random spline + `lcmm::hlme`

Produced by [`../../analysis_random_spline.R`](../../analysis_random_spline.R).

```bash
python3 prepare_monolix_data.py
Rscript analysis_random_spline.R      # ~13 min (R3 dominates the runtime)
```

For what each CSV column means, see [`notice.md`](notice.md).

> **Read this one with iteration 03 open.** It fixes 03's main defect and breaks 03's main
> guarantee. Neither iteration satisfies both requirements at once, and that tension is the
> real result of this iteration.

---

## Why this iteration exists

[Iteration 03](../03_natural_spline_hlme/) left one problem open. Its within-shoulder
residual autocorrelation was **0.978 and unmodelled**, because a shoulder's deviation from
the population curve is itself a smooth function of elevation, and a random *intercept*
absorbs only its level. Its own README said the honest fix was a random *smooth* per
shoulder.

That is the canonical `lcmm` remedy, and the pattern used in the reference script by the
package's own authors (Prague & Proust-Lima): put the whole spline basis in the random
effects.

Two changes from iteration 03, both adopted from that reference:

| | iteration 03 | iteration 04 |
| --- | --- | --- |
| random part | `~ 1` | **`~ ns(TIME, knots)`** — a random curve per shoulder |
| knots | `df = 5`, data-driven quantiles | **explicit**: interior 40/70/100/130°, boundary 0/160° |
| random effects | 1 | **6** (intercept + 5 basis terms) |

The explicit knots matter for interpretation, not just reproducibility. Iteration 03's
boundary knots landed at −3.46° and 187.02°, so its $\gamma_0$ was the difference at an
*extrapolated* point and flipped sign with `df`. Anchoring the lower boundary at 0° makes
$\gamma_0$ the in-vivo minus ex-vivo difference **with the arm at the side** — a real
quantity. 76 of 3,724 rows fall outside [0°, 160°] and get the linear extrapolation a
natural spline is designed for.

## The model

```r
hlme(fixed  = Y ~ ns(TIME, knots = c(40,70,100,130), Boundary.knots = c(0,160)) * cond,
     random = ~ ns(TIME, knots = c(40,70,100,130), Boundary.knots = c(0,160)),
     subject = "IDnum", ng = 1)
```

$$
\begin{aligned}
y_{ij} =\;& \underbrace{\beta_0 + \sum_{k=1}^{5}\beta_k B_k(x_{ij})}_{\text{ex-vivo reference}}
 + \underbrace{\Big[\gamma_0 + \sum_{k=1}^{5}\gamma_k B_k(x_{ij})\Big]\mathbb{1}(c_i=\text{in vivo})}_{\text{in-vivo departure}} \\[4pt]
&+ \underbrace{b_{i0} + \sum_{k=1}^{5} b_{ik} B_k(x_{ij})}_{\text{this shoulder's own curve}}
 + \varepsilon_{ij},
\qquad \mathbf{b}_i \sim \mathcal{N}_6(\mathbf{0}, D),
\quad \varepsilon_{ij}\sim\mathcal{N}(0,\sigma_\varepsilon^2).
\end{aligned}
$$

The fixed part is **unchanged** from iteration 03 — same 12 coefficients, same meaning. The
only difference is $\mathbf{b}_i$: six random effects instead of one, so each shoulder gets
its own *shape*, not merely its own height.

---

## The comparison

`00_variants.csv`. `resid_acf1` is the metric this iteration exists to move.

| variant | random part | params | BIC | lag-1 ACF | resid sd | kurtosis | Wald | time |
| --- | --- | ---: | ---: | ---: | ---: | ---: | :--- | ---: |
| `R1_intercept` | `~1` (iteration 03) | 14 | 19510 | **0.978** | 3.170 | **1.18** | 8/12 | 11 s |
| `R2_spline_diag` | `~ns(…)`, diagonal | 19 | 6449 | 0.836 | 0.437 | 12.55 | 7/12 | 16 s |
| `R3_spline_unstr` | `~ns(…)`, unstructured | 34 | **6348** | **0.835** | 0.433 | **12.98** | 7/12 | 659 s |

### What improved — a great deal

The random curve is unambiguously the better description of the data. **BIC drops from
19,510 to 6,348**, and the residual sd falls from 3.17° to 0.43°: most of what iteration 03
called "residual" was in fact each shoulder's own trajectory shape.

The autocorrelation headline number, 0.978 → 0.836, badly *undersells* the improvement.
Look at `05_autocorrelation.png`: under iteration 03 the ACF was a slow persistent drift,
still near 1 many lags out. Now it **decays to zero within about six lags**, swings negative
to −0.38 around lag 11, and damps out — a short-range oscillation instead of a long-range
drift. The integrated autocorrelation, which is what actually inflates effective sample
size, has collapsed even though lag-1 is still high.

### What broke — the acceptance criterion

**Residual excess kurtosis goes from 1.18 to 12.98.** `03_diagnostics.png` shows it plainly:
the histogram is a needle at zero, and the Q-Q plot is severely S-shaped with residuals
reaching ±10 standard deviations.

The mechanism is the same one that discredited `cor = AR(TIME)` in iteration 03. As the
random structure absorbs more of each shoulder's trajectory, $\varepsilon_{ij}$ stops being
measurement noise and becomes *the crumbs left after a nearly perfect fit* — sd 0.43° on
data spanning 60°. What remains is dominated by the handful of points the random curve
cannot reach, which is exactly a heavy-tailed distribution.

So the requirement that $\varepsilon_{ij}$ be normal and the requirement that the
within-curve correlation be modelled **pull against each other** here:

| | autocorrelation handled | residuals normal |
| --- | :---: | :---: |
| iteration 03 (`~1`) | ✗ (ACF 0.978) | ✓ (kurtosis 1.18) |
| iteration 04 (`~ns`) | mostly ✓ (short-range) | ✗ (kurtosis 12.98) |

### Does a smaller basis rescue the normality? No.

The obvious hope is that the kurtosis is a side-effect of using `df = 5`, and that a
smaller basis — fewer random effects, less absorbed — would leave a larger, better-behaved
residual. `00_df_sweep.csv` tests it directly, refitting the **diagonal** random spline at
every df with knots from one rule (quantiles of the fitted x rounded to the nearest 10°,
which reproduces 40/70/100/130 at K = 5):

| K | random effects | knots | BIC | lag-1 ACF | resid sd | excess kurtosis |
| ---: | ---: | :--- | ---: | ---: | ---: | ---: |
| 2 | 3 | 90 | 12361 | 0.938 | 1.108 | **5.70** |
| 3 | 4 | 60/110 | 9248 | 0.899 | 0.694 | 6.73 |
| 4 | 5 | 50/90/120 | 7607 | 0.871 | 0.532 | 9.44 |
| 5 | 6 | 40/70/100/130 | **6449** | **0.836** | 0.437 | **12.55** |

The mechanism is confirmed exactly, and it is **monotone in both directions at once**
(`12_by_df_tradeoff.png`):

- more basis functions → more of each shoulder's trajectory absorbed → **residual sd falls**
  1.108 → 0.437,
- and as the residual shrinks toward nothing, **kurtosis climbs** 5.70 → 12.55,
- while the thing we wanted, **autocorrelation, improves** 0.938 → 0.836.

So every step that buys autocorrelation costs normality, in lockstep. There is **no K at
which both are acceptable**. Even the most restrictive model tried, K = 2 with three random
effects, still has excess kurtosis 5.70 — five times worse than *any* model in iteration 03,
which sat at 1.03–1.19 for every K — and its ACF of 0.938 is barely better than iteration
03's 0.978. Meanwhile BIC punishes the retreat heavily: 12,361 at K = 2 against 6,449 at
K = 5.

The contrast with iteration 03 is what settles it. Under a random **intercept**, df is
irrelevant to normality — residual sd stays ≈ 3.2 and kurtosis ≈ 1.1 at every K. Under a
random **curve**, df merely slides you along the trade-off. The kurtosis is therefore not a
df problem to be tuned away; it is intrinsic to letting each shoulder have its own curve on
data this dense, and the answer has to come from somewhere else — a residual model that
expects heavy tails, or fewer points per shoulder so that ε<sub>ij</sub> is genuinely
measurement noise again.

### The curves and the difference across df

`07_by_df_population_curves.png` and `08_by_df_difference.png` repeat the retained model's
two headline views at every df, mirroring iteration 03's pair.

**The population curves show an asymmetry worth noting.** The ex-vivo curve develops its
S-bend as df rises, exactly as in iteration 03 — but the **in-vivo curve stays close to a
straight line at every df**, from K = 2 through K = 5. Extra basis functions buy shape almost
entirely for the ex-vivo condition. That is consistent with the ex-vivo data being far
denser per shoulder even after thinning (1,101 rows across 13 shoulders versus 2,623 across
31), so the in-vivo curve is being estimated from many shoulders each contributing few, more
scattered points.

**The difference curves are the honest version of iteration 03's.** The mean absolute
difference is essentially invariant — **9.95, 10.36, 10.01, 9.89°** at K = 2…5 — matching
iteration 03's stable ~10° and confirming that this quantity does not depend on the basis or
on the random-effect structure. It is the number to report.

What *has* changed is the uncertainty. The bands are roughly twice as wide as iteration
03's, because allowing each shoulder its own curve removes the pseudo-replication. The
consequence is visible along the bottom of each panel:

| K | mean \|diff\| | band excludes 0 over | high-elevation behaviour |
| ---: | ---: | :--- | :--- |
| 2 | 9.95° | 27–150° (68%) | crosses 0 at 135°, reaches +8.3° |
| 3 | 10.36° | 25–150° (77%) | crosses 0 at 130°, reaches +13.9° |
| 4 | 10.01° | 16–120° (76%) | crosses 0 at 136°, reaches +4.9° |
| **5** | **9.89°** | **14–122° (79%)** | reaches 0 only at 150° |

So the difference is firmly established across the **low-to-mid range** and becomes
indistinguishable from zero above roughly 120°, at every df. The sign-reversal artefact that
iteration 03 showed at low df is still visible here at K = 2 and 3 (+8.3° and +13.9° at the
top end) and still absent at K = 5 — but now the confidence band covers zero there anyway,
so it no longer supports any claim in either direction. That is the cleaner reading: **the
end-range behaviour is simply unresolved**, rather than being a difference whose sign depends
on the basis.

### Which coefficients survive, at every df

`00_wald_tests_by_df.csv` gives the full per-coefficient table for all four df (36 rows:
6 + 8 + 10 + 12, with `K` and `retained_df` columns), and `11_by_df_wald.png` draws it as a
grid — rows are coefficients, columns are df, the rule separates reference-curve terms
(below) from difference-curve terms (above). Per-curve counts are also columns of
`00_df_sweep.csv`:

| K | reference | difference |
| ---: | :--- | :--- |
| 2 | 3/3 | 2/3 |
| 3 | 4/4 | 3/4 |
| 4 | 4/5 | 2/5 |
| 5 | 5/6 | 2/6 |

The grid says considerably more than the counts, and two things in it are stable enough to
report:

- **γ<sub>0</sub> (`condin vivo`) never passes, at any df** — |Wald| = 1.5, 1.2, 1.0, 1.9.
  With the boundary knot at 0°, that is a clean and consistent statement: **the in-vivo and
  ex-vivo curves are not distinguishable in level with the arm at the side.** In iteration 03
  this coefficient flipped sign with df and could not be read at all; here it is stable, and
  stably null.
- **`ns1:condin vivo` passes at every df** — 5.4, 5.3, 9.1, 7.6, always the strongest
  difference term. This is the low-elevation shape departure, and it is the single most
  robust piece of evidence for a condition effect in this iteration.

Everything else is unstable. `ns2:condin vivo` passes at K = 2, 3 and 5 but collapses to 0.3
at K = 4; `ns3:condin vivo` passes only at K = 3; `ns4:condin vivo` only at K = 4. Even on
the reference curve the intercept passes at K = 2 and 3 (2.9, 2.6) then fails at K = 4 and 5
(0.1, 0.5).

That instability is not noise in the data — the higher-order basis coefficients are strongly
correlated with each other, and which one carries a given feature of the curve shifts as the
knots move. It is a good reason **not** to report individual γ<sub>k</sub> as findings, and
to lean instead on the quantity that does not move: the mean |difference| of ≈ 10°.

### Are the b<sub>i</sub> normal — is that tested, or assumed?

**Assumed.** $\mathbf{b}_i \sim \mathcal{N}_6(\mathbf{0}, D)$ is part of the likelihood
`hlme` maximises; nothing in the fit checks it, and there is no test statistic for it in the
output. `04_random_effects.png` is a *diagnostic* we chose to draw, not something the model
validates.

That figure is also **weaker evidence than it looks**, for a reason worth knowing.
`predRE` returns empirical Bayes predictions (BLUPs), and those are **shrunk toward zero**
— by a factor that depends on how much data each shoulder contributes. Shrinkage pulls them
toward normality, so a straight Q-Q of predicted random effects is optimistically biased:
it can look normal even when the underlying $\mathbf{b}_i$ are not. Use it to spot outlying
shoulders, not to confirm the assumption.

Two things follow:

- **Mis-specifying the b<sub>i</sub> distribution is usually less damaging than
  mis-specifying $\varepsilon_{ij}$.** The fixed effects $\beta,\gamma$ stay consistent under
  a fairly wide range of random-effect distributions; what degrades is $D$ itself and the
  per-shoulder predictions. So the kurtosis problem on $\varepsilon_{ij}$ documented above is
  the more serious of the two.
- **If the b<sub>i</sub> really are non-normal, `lcmm`'s own answer is `ng > 1`.** Latent
  classes are precisely how this package represents a random-effect distribution that is not
  a single Gaussian — a mixture over sub-populations of shoulders instead. That is a
  different scientific claim (it says there are *kinds* of shoulder) and is out of scope
  here, but it is the natural next step if figure 04 ever shows real structure rather than
  mild tail flattening.

### What it says about iteration 03's Wald tests

Iteration 03 warned its standard errors were optimistic. They were. Once shoulders are
allowed their own shapes, the **difference curve drops from 4/6 to 2/6** passing terms —
only $\gamma_1$ (Wald −6.99) and $\gamma_2$ (−4.69) survive. $\gamma_3$ falls to 1.75,
$\gamma_4$ to −0.14, $\gamma_5$ to 0.55.

Read plainly: **the evidence for a difference in trajectory *shape* is weaker than iteration
03 suggested**, and concentrated in the low-to-mid elevation range. The reference curve is
untouched at 5/6.

$\gamma_0$ is now interpretable thanks to the boundary knot at 0°: the in-vivo curve sits
**−4.06° ± 2.31** below ex-vivo with the arm at the side (Wald −1.76, p = 0.079) — a real
anatomical statement, though not distinguishable from zero.

### What "diagonal" and "unstructured" mean

Both variants assume $\mathbf{b}_i = (b_{i0},\dots,b_{i5})^{\top} \sim \mathcal{N}_6(\mathbf{0}, D)$.
They differ only in what $D$ is allowed to be.

**Unstructured** ($D$ free): any symmetric positive-definite $6\times6$ matrix, so
$6\cdot7/2 = 21$ free parameters — 6 variances and 15 covariances.

$$
D=\begin{pmatrix}
d_{00} & d_{01} & \cdots & d_{05}\\
d_{01} & d_{11} & \cdots & d_{15}\\
\vdots & & \ddots & \vdots\\
d_{05} & d_{15} & \cdots & d_{55}
\end{pmatrix}
$$

The off-diagonals carry real content: $d_{01}\neq 0$ says a shoulder sitting high overall
*also* tends to have a particular mid-range curvature. Shoulders vary along correlated
directions rather than independently.

**Diagonal** (`idiag = TRUE`): $D = \operatorname{diag}(d_{00},\dots,d_{55})$, just 6
parameters. Each random coefficient still varies between shoulders, but knowing one tells
you nothing about the others — the $b_{ik}$ are mutually independent. It is a strict
restriction of the unstructured model (15 covariances set to 0), so the two are nested and
comparable by likelihood ratio.

> **One subtlety worth knowing.** "Diagonal" is not a basis-free statement. If the basis is
> reparametrised, $B' = BA$ for an invertible $A$, then $\mathbf{b}' = A^{-1}\mathbf{b}$ has
> covariance $A^{-1}DA^{-\top}$ — generally *not* diagonal. So `idiag = TRUE` imposes
> independence on **the particular coefficients `ns()` happens to return**, and a different
> though mathematically equivalent spline parametrisation would impose a different
> assumption. The unstructured model has no such dependence: it is invariant to how the
> basis is written. That is the real argument for unstructured — not fit, but that the
> assumption means the same thing regardless of parametrisation.

### Diagonal or unstructured?

BIC selects `R3_spline_unstr`, and the script follows BIC, so every other output here is R3.
But that is a thin win: **101 BIC units for 15 extra parameters and 40× the runtime**
(659 s vs 16 s), for a lag-1 ACF identical to three decimals (0.8352 vs 0.8356) and a
*slightly worse* kurtosis. On 44 shoulders, an unstructured 6×6 covariance is 21 parameters
estimated from 44 units.

**For the 72-cell sweep, use `R2_spline_diag`.** R3 is not worth 40× the compute for a
difference you cannot see in any figure.

---

## Verdict

Iteration 04 is the better *model* and iteration 03 is the better *inference*, and that is
not a satisfying place to stop.

- The random curve is clearly right about the data: BIC, residual sd and the ACF structure
  all say so, and it corrects an over-optimism in iteration 03 that we had flagged but not
  quantified.
- But its $\varepsilon_{ij}$ no longer satisfy the normality requirement, so its own Wald
  standard errors rest on an assumption its diagnostics reject.

Neither is the finished answer. What both agree on, and what should be reported:

- the ex-vivo reference trajectory (5/6 Wald in both),
- a **shape** difference concentrated at low-to-mid elevation ($\gamma_1$, $\gamma_2$ pass
  in both, with large margins),
- a level difference at 0° that is **not** distinguishable from zero.

And the standing caveat is unchanged: condition is perfectly confounded with source study,
so all of this is descriptive, not causal — see
[`../../docs/STATISTICAL_METHODS_REVIEW.md`](../../docs/STATISTICAL_METHODS_REVIEW.md).

---

## The files

| file | what it is |
| --- | --- |
| `00_results_summary.txt` | human-readable digest |
| `00_variants.csv` | **the decisive table** — the three random-effect specifications compared |
| `00_df_sweep.csv` | the df sweep: K = 2…5 with the diagonal random spline, ACF / residual sd / kurtosis per K |
| `00_coefficients.csv` | all 34 parameters of the retained model with SE and Wald |
| `00_wald_tests.csv` | the retained model's 12 fixed effects, grouped by curve, with pass/fail |
| `00_wald_tests_by_df.csv` | the same, for **every** df — with `K` and `retained_df` columns |
| `01_population_curves_by_condition.png` | fitted curves with 95% bands, knots rugged, Wald counts |
| `02_difference_invivo_minus_exvivo.png` | in-vivo − ex-vivo difference on the x-overlap |
| `03_diagnostics.png` | **where the normality fails** — Q-Q, histogram, residuals vs fitted and vs x |
| `04_random_effects.png` | Q-Q of each of the 6 random effects |
| `05_autocorrelation.png` | the residual ACF against iteration 03's, as a reference line |
| `06_variant_comparison.png` | ACF and BIC across the three specifications |
| `07_by_df_population_curves.png` | **2×2: fitted curves at K = 2…5**, own axes, knots rugged, retained df boxed |
| `08_by_df_difference.png` | **2×2: the in-vivo − ex-vivo difference at each K**, shared y-axis, band-excludes-zero marked |
| `09_by_df_diagnostics.png` | one row per df: Q-Q + histogram — **the evidence that no df rescues normality** |
| `10_by_df_random_effects.png` | **2×2: Q-Q of the random intercept at each K**, shared axes |
| `11_by_df_wald.png` | **grid of \|Wald\| per coefficient × df** — which terms survive, and at which df |
| `12_by_df_tradeoff.png` | kurtosis, ACF and residual sd against df — the trade-off in three panels |
