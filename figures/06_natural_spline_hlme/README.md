# Iteration 06 — Melanie Prague's plan, and what it found

Produced by [`../../analysis_lrt.R`](../../analysis_lrt.R) in ~5 minutes.

```bash
python3 prepare_monolix_data.py
Rscript analysis_lrt.R
```

> **The plan below was executed, then step 1 was re-run over K = 1…20.**
> [Results](#results) are at the bottom. On the K = 1…10 grid originally specified, BIC never
> turned and the retained K was the edge of the grid rather than a minimum. **Extending the
> grid to 20 fixed that**: BIC has a genuine interior minimum at **K = 10**, and that single
> criterion is now the whole selection rule — see
> [how K is chosen](#how-k-is-chosen--one-method-one-number). A separate experiment shows
> what that K is really a statement about:
> [hold the random part fixed](#the-experiment-that-proves-it-hold-the-random-part-fixed).

Same cell as [iteration 04](../04_natural_spline_hlme/): **scapulothoracic / frontal plane
elevation / DoF 2**, 22,102 observations from 44 shoulders.

---

## Why this iteration exists

Iterations 03, 04 and [`../generalized_v2/`](../generalized_v2/) selected the spline df by
BIC on the *full* model and reported **Wald tests per coefficient**. Those iterations
established, and then documented at length, that the per-coefficient route is unstable:
which coefficient carries the signal moves as the knots move, so a cell can read `○○○○○`
and still be strongly significant jointly. Melanie's verdict is blunter than ours was:

> **4) PAS de test de wald, on abandonne !** ce n'est pas assez informatif et trop dur à
> interpréter ! Regarder les graphs de diff, c'est plus informatif, c'est ce sur quoi
> j'étais parti, c'est parfait !

So the inference changes shape. The **likelihood-ratio test** replaces the Wald battery, the
df is chosen on a **reduced** model rather than the full one, and the **difference curve**
becomes the thing we report.

Three further changes:

| | iterations 03–05 | **iteration 06** |
| --- | --- | --- |
| reference condition | ex vivo | **in vivo** — so the departure is the *ex-vivo* departure |
| df selection | BIC on the full model, all data, K = 2…5 | **BIC on the IN-VIVO data alone, K = 1…20** |
| global test | joint Wald χ² on the γ block | **likelihood-ratio test** |
| per-coefficient tests | Wald + ●/○ strip | **none** |
| knots on figures | dashed verticals + dots | **not drawn** |

---

## The model

With **in vivo as the reference**, for shoulder *i*, observation *j*, condition
$c_i \in \{\text{in vivo},\ \text{ex vivo}\}$:

$$
y_{ij} = \underbrace{\beta_0 + \sum_{k=1}^{K}\beta_k B_k(x_{ij})}_{\textbf{in-vivo reference curve}}
+ \underbrace{\Big[\gamma_0 + \sum_{k=1}^{K}\gamma_k B_k(x_{ij})\Big]\mathbb{1}(c_i=\text{ex vivo})}_{\textbf{ex-vivo departure}}
+ \underbrace{b_{i0} + \sum_{k=1}^{K} b_{ik}B_k(x_{ij})}_{\text{this shoulder's own curve}}
+ \varepsilon_{ij}
$$

with $\mathbf{b}_i \sim \mathcal{N}_{K+1}(\mathbf{0}, \Omega)$, $\Omega$ **diagonal**, and
$\varepsilon_{ij}\sim\mathcal{N}(0,\sigma_\varepsilon^2)$.

$\Omega$ is the covariance matrix of the random effects — written $\Omega$ rather than $D$,
following the usual convention in the mixed-model literature. Diagonal here means
$\Omega = \operatorname{diag}(\omega_0^2,\dots,\omega_K^2)$: each shoulder's basis
coefficients vary independently, so $K+1$ variances are estimated and no covariances.

> ### ⚠ The difference curve changes sign
> Every earlier iteration plotted **in vivo − ex vivo**. Because in vivo is now the
> reference, this one plots **ex vivo − in vivo** — the mirror image. A reader comparing
> figure 03 here against iteration 04's without noticing will read every conclusion
> backwards. The axis label says so explicitly.

---

## The four steps

### 1 — Choose K by BIC on the **in-vivo data alone**, K = 1…20

**in vivo is the reference condition, so the shape of the reference trajectory is chosen
from the reference's own data.** The ex-vivo shoulders get no vote on how many degrees of
freedom the common curve needs — they are what will be tested against it.

```r
din <- subset(data, cond == "in vivo")
hlme(fixed = Y ~ ns1 + … + nsK,          # no condition term: there is only one condition
     random = ~ ns1 + … + nsK,
     subject = "IDnum", ng = 1, idiag = TRUE, data = din)
```

31 of the 44 shoulders are in vivo. The basis is still built on the **full trimmed range**,
so the K chosen here transfers unchanged to the H₀ and H₁ models of steps 2–3, which are
fitted on every shoulder.

> **H₀ is refitted on all the data.** Step 1's in-vivo model *chooses K*; it is not the null
> hypothesis. The likelihood-ratio test needs H₀ and H₁ on identical rows, so H₀ —
> `Y ~ ns(df = K)` with no condition terms — is fitted again on all 44 shoulders at the
> retained K. Differencing step 1's in-vivo log-likelihood against a full-data H₁ would be
> meaningless.

The random part **mirrors K**, as in iteration 04 — so K = 10 means 11 random effects per
shoulder. Each K records `converged`, and selection runs over the converged rows only.
**This mirroring turns out to be the whole story**; see
[what "K = 10" is really a statement about](#what-k--10-is-really-a-statement-about).

> *Anticipated, and half wrong:* the plan expected high K to fail outright. **All twenty
> reported convergence**, the slowest in 88 s, at K = 20 asking 21 diagonal variances from
> 31 shoulders — nothing errored, which is precisely why BIC kept going on the short grid.
> But convergence is not the same as finding the global optimum, and above K ≈ 12 the
> likelihood itself degrades, so the prediction was right in spirit: the model does run out
> of road, it just does so quietly instead of crashing. See
> [the tail caveat](#read-the-tail-with-care--the-fit-itself-gets-worse-past-k--10).

`ns(df = K)` gives K basis columns and K − 1 interior knots; K = 1 is a straight line,
K = 10 has 9 interior knots.

#### Per-shoulder RMSE at every K, not just the retained one

BIC is a single number about the whole fit; it says nothing about whether some individual
shoulder is being left behind. So the sweep also records, at **every** K, the per-shoulder
RMSE across the in-vivo shoulders — its **maximum** and its **median** — and
`01_choosing_k.png` plots them beside the BIC curve on a log scale.

The two panels answer different questions:

| panel | question |
| --- | --- |
| BIC vs K | how much structure does the *likelihood* want, net of the parameter penalty |
| max / median RMSE vs K | is the *worst-fitted* shoulder still being missed, and by how much |

### 1a — Individual-fit check, before the condition enters

> *trouver un indicateur pour vérifier que le max d'erreur chez un individu n'est pas trop
> grand par rapport aux autres*

From the subject-specific (BLUP-based) predictions `m$pred$pred_ss` of **step 1's in-vivo
model** at the retained K, compute per shoulder *i*: `n_i`, `rmse_i`, `max_abs_i`. Two
summaries answer the question:

- **Headline ratio** — $\max_i(\text{rmse}_i) / \operatorname{median}_i(\text{rmse}_i)$.
  One number: *the worst-fitted shoulder is X times the typical one.*
- **Robust flag** — RMSE is right-skewed, so a mean ± sd rule would be dragged by the very
  outliers it is meant to catch. Use a MAD-based z on the log scale:

$$
z_i = \frac{\log \text{rmse}_i - \operatorname{median}_j(\log \text{rmse}_j)}
{1.4826 \cdot \operatorname{MAD}_j(\log \text{rmse}_j)}, \qquad \text{flag } z_i > 3
$$

This runs on the **in-vivo** model deliberately: a badly fitted shoulder found here is a
data or model-adequacy problem inside the reference condition itself, which no amount of
condition modelling could explain away.

#### Look at the dots: $y_{ij}$ against $\hat{y}_{ij}$

Summary statistics can hide a shoulder that is fitted with the right *average* error but the
wrong *shape*, so figure 04 plots every observation against its prediction, with the 1:1
line. `hlme` gives two predictions per row, and **the contrast between them is the whole
point**:

| prediction | what it uses | what a tight plot means |
| --- | --- | --- |
| $\hat{y}^{\,\text{marg}}_{ij}$ = `pred_m` | the **population** curve only, $\mathbf{b}_i = \mathbf{0}$ | the shared trajectory alone explains the data |
| $\hat{y}^{\,\text{ss}}_{ij}$ = `pred_ss` | population curve **+ that shoulder's BLUP** $\mathbf{b}_i$ | each shoulder is individually captured |

Plotted side by side, the marginal panel should be visibly looser than the subject-specific
one — the gap between them *is* the work the random curves are doing, and it is exactly what
Melanie's note (a) asks to see. If the subject-specific panel is still scattered, the random
structure is not capturing individuals and no amount of condition modelling will fix that.

Each panel gets the 1:1 line (not a regression line — the question is agreement, not
correlation), points coloured by condition, and the concordance
$R^2 = 1 - \frac{\sum (y-\hat{y})^2}{\sum (y-\bar{y})^2}$ annotated. Reporting $R^2$ this
way rather than as a squared correlation matters here: a shoulder fitted with a constant
offset can correlate at 0.99 while sitting well off the 1:1 line, and only the residual form
catches it.

### 2 — Full model at the **same** K, in vivo as reference

```r
d$cond <- factor(..., levels = c("in vivo", "ex vivo"))   # in vivo FIRST = reference
hlme(fixed = Y ~ (ns1 + … + nsK) * cond, random = ~ ns1 + … + nsK, …)
```

The K is inherited from step 1 and **not re-selected**. The script asserts the coefficient
is named `condex vivo` — if it comes back `condin vivo` the reference never flipped and
every sign downstream is wrong.

### 3 — Likelihood-ratio test

```r
LR <- 2 * (m1$loglik - m0$loglik)
df <- K + 1                       # gamma_0 plus K shape terms
p  <- pchisq(LR, df = df, lower.tail = FALSE)
```

- `p < 0.05` → reject H₀: the conditions differ.
- `p > 0.05` → cannot reject H₀.

df = K + 1 matches Melanie's `df = 6` worked example at K = 5. Both models must be fitted on
**exactly the same rows** — the script asserts equal row and subject counts before
differencing the log-likelihoods, since a silent row drop would make the statistic
meaningless rather than merely wrong.

**No multiple-comparison adjustment here**: iteration 06 is one cell, one test. Bonferroni
across the whole plate belongs to the generalized rewrite (see below).

### 4 — No Wald tests

Dropped entirely: no `00_wald_tests.csv`, no ●/○ strip, no per-coefficient stars. The
coefficient table is still exported for the record, but carries **estimates and standard
errors only** — no `wald`, no `passed` column, so nobody is tempted to reintroduce them.

The reporting object is **figure 03, the difference curve with its confidence band**.

---

## Figures — no knots drawn anywhere

| figure | content |
| --- | --- |
| `01_choosing_k.png` | **BIC against K = 1…20** beside per-shoulder max/median RMSE — the df selection, shown rather than asserted. The chosen K is a filled dot; non-converged K are marked distinctly. If the test ever runs at a K below BIC's minimum (a cost cap), *both* are drawn and labelled, so an off-minimum dot never reads as an error |
| `02_population_curves.png` | both conditions with 95% bands over the raw scatter, stacked vertically for **every** K |
| `03_difference.png` | **ex vivo − in vivo** with 95% CI and the zero line, at every K; the headline |
| `04_observed_vs_predicted.png` | **2 panels: $y_{ij}$ vs $\hat{y}^{\,\text{marg}}_{ij}$ and $y_{ij}$ vs $\hat{y}^{\,\text{ss}}_{ij}$**, 1:1 line, coloured by condition, concordance $R^2$ annotated — the gap between the panels is what the random curves buy |
| `05_individual_fit.png` | per-shoulder RMSE sorted with the flag threshold drawn, and observed-vs-predicted **per shoulder** for the worst few — one small panel each, so a wrong *shape* is visible and not just a large average error |
| `06_diagnostics.png` | residuals vs fitted, Q-Q, histogram, residuals vs x |
| `08_bic_two_penalties.png` | BIC computed on **subjects** (as `lcmm` does) against BIC on **observations**, with the max/median RMSE panel beside it — shows the choice of K does not hinge on what *n* means |
| `09_what_the_df_buys.png` | deviance and BIC against fixed K for three random-part sizes — the experiment showing the fixed df buy ~0.2% of what the random dimension does |

`08_*` and `09_*` come from [`../../analysis_random_dim.R`](../../analysis_random_dim.R),
which reads `00_bic_by_k.csv` and so must be run **after** `analysis_lrt.R`:

```bash
Rscript analysis_lrt.R          # ~19 min — steps 1-4, figures 01-06
Rscript analysis_random_dim.R   # ~13 min — figures 08-09
```

> `07_*` was the L-method knee analysis. **Dropped** — one method decides K, and it is BIC.
> Those files are deleted. `utils_lmethod.R` and `analysis_k_grid.R` are left on disk but are
> now unreferenced dead code; `analysis_k_grid.R` in particular opens by asserting that BIC
> never turns, which the K = 1…20 grid disproves. Delete both when convenient.

Every figure as high-resolution PNG **and** vector PDF (`cairo_pdf`, so Unicode glyphs
survive), reusing the renderer from `../../analysis_random_spline.R`.

## Outputs

| file | contents |
| --- | --- |
| `00_bic_by_k.csv` | **the selection table** — K = 1…20 with n_fixed, n_random, converged, loglik, BIC, per-shoulder rmse max/median/ratio, marginal and subject-specific $R^2$, seconds |
| `00_lrt.csv` | loglik0, loglik1, LR, df, p, and the row/subject counts the assertion checked |
| `00_lrt_by_k.csv` | the same test repeated at **every** K in `KS_FULL` — the evidence that the verdict is not a df artefact |
| `00_coefficients.csv` | both models' estimates and SEs — **no Wald column** |
| `00_individual_fit.csv` | per shoulder: n, rmse, max_abs, ratio to median, robust z, flagged, and the concordance $R^2$ under both the marginal and the subject-specific prediction |
| `00_results_summary.txt` | digest |
| `08_bic_by_k_extended.csv` | `00_bic_by_k.csv` plus the recovered parameter count and BIC recomputed on the **observation** count |
| `09_decoupled_random.csv` | the decoupling experiment: fixed K = 1…20 crossed with a random part held at 2 and at 4 random effects |

---

## Settings and open items

**`TRIM_TO_OVERLAP <- TRUE`.** Melanie's TODO — *garder toute la plage des abscisses même
si que deux individus* — refers to a degenerate case in `analyze_all_v2`, where the trim cut
to the common portion of just two curves. That does not apply to this cell (13 ex-vivo and
31 in-vivo shoulders), so iteration 06 keeps iteration 04's overlap trim and stays directly
comparable to it. Flipping the constant to `FALSE` fits the full union range instead.

**Carried to the generalized rewrite, not done here:**

- **Bonferroni across the plate** — Melanie asked for Bonferroni rather than BH, being the
  more conservative choice, applied to the LRT p-values of all compared cells at once.
- **Per-marker df** — each cell selects its own K by BIC on its own reduced model, so the
  spline degree varies by marker instead of being fixed at 4.
- **The full-range question above**, which is where it actually bites.
- **Knot placement at quantiles** — marked optional in the notes; already what the current
  knot rule does.

## What carries over unchanged

The data preparation, the thinning to ≤100 points per shoulder, the diagonal random-curve
structure and the `cairo_pdf` renderer all come from
[`../../analysis_random_spline.R`](../../analysis_random_spline.R). The standing caveat also
carries over: **condition is perfectly confounded with source study**, so whatever the LRT
says, the result is descriptive and not causal — see
[`../../docs/STATISTICAL_METHODS_REVIEW.md`](../../docs/STATISTICAL_METHODS_REVIEW.md).

---

## Results

**The likelihood-ratio test is decisive.** At the selected **K = 10**, on 3,378 rows from all
44 shoulders:

| | |
| --- | ---: |
| log-likelihood, reduced (H₀) | −726.66 |
| log-likelihood, full (H₁) | −674.21 |
| **LR = 2(LL₁ − LL₀)** | **104.90** |
| df = K + 1 | 11 |
| **p** | **1.9 × 10⁻¹⁷** |
| mean \|ex vivo − in vivo\| | **10.11°** over 14–150° |

→ **reject H₀: the ex-vivo trajectory departs from the in-vivo reference.**

The test now runs at the K the method selected, with no cap in between. An earlier run
reported K = 8 here: `KS_FULL` had been capped at 8 on a measurement of **3,425 s** for the
K = 10 pair of fits. That measurement was an outlier — **the same pair now takes 91 s**, and
the whole K = 1…10 sweep 426 s. The cap was removed.

The mean |ex vivo − in vivo| difference over 14–150° is in
[`00_lrt.csv`](00_lrt.csv) and drawn in [`03_difference.png`](03_difference.png), which is
the reporting object Melanie asked for. Note the sign: **ex vivo minus in vivo**, mirrored
against every earlier iteration.

### The test does not depend on the K argument at all

The whole K debate below turns out not to touch the conclusion. H₀ and H₁ were refitted on
all 44 shoulders at **every** K from 1 to 10 ([`00_lrt_by_k.csv`](00_lrt_by_k.csv)), and the
test rejects everywhere:

| K | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | **10** |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| LR | 24.7 | 42.8 | 54.2 | 68.1 | 83.1 | 100.6 | 132.7 | **144.4** | 129.5 | 104.9 |
| df | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
| p | 4e−06 | 3e−09 | 5e−11 | 2e−13 | 8e−16 | 8e−19 | 8e−25 | **1e−26** | 6e−23 | 2e−17 |

Even K = 1 — a straight line for each condition, the crudest model on the grid — rejects at
p = 4 × 10⁻⁶. **K changes the curve you draw, not the verdict you reach.** So the sections
that follow are about getting an honest *description* of the trajectory; they are not
load-bearing for the finding.

**The LR is not monotone in K**: it climbs to 144.4 at K = 8, then falls back to 129.5 and
104.9. Extending the sweep from 8 to 10 is what exposed that — on K ≤ 8 alone the statistic
looks like it rises without limit. The reason is that H₀ and H₁ both gain flexibility as K
grows, and past K ≈ 8 the *shared* curve starts absorbing departure that the condition block
had been carrying. The selected K = 10 therefore sits on the **conservative** side of the
grid, which is the right direction to err.

### How K is chosen — one method, one number

> **The method.** Fit the reduced model (no condition terms) to the **in-vivo data alone**
> over **K = 1…20**, with the random part mirroring K, and **take the BIC minimum.**
> → **K = 10.**
>
> That is the whole rule. It has no tuning constant and no judgement call, it is Melanie's
> step 1 as written, and it scores the model that is actually fitted downstream — steps 2–3
> mirror the random part the same way, so BIC is grading the real thing rather than a proxy.

Earlier drafts of this file reported three candidate K side by side — a BIC minimum, an RMSE
plateau, and an L-method knee — without saying which one decided. That was the confusing part,
and it is resolved: **BIC decides.** The L-method has been dropped entirely and its outputs
deleted. The per-shoulder RMSE stays below, as a description of what K = 10 buys, not as a
rival answer.

The one objection to BIC was real and is now gone. On the original K = 1…10 grid BIC fell
monotonically and its minimum sat on the last value tried — it reported the edge rather than
making a choice. **Extending to K = 20 removes that objection: BIC turns.** All twenty
reported convergence — though see [the tail caveat](#read-the-tail-with-care--the-fit-itself-gets-worse-past-k--10)
for why that is weaker than it sounds at high K.

| K | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | **10** | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BIC | 9407 | 6280 | 4322 | 3323 | 2419 | 2144 | 1871 | 1710 | 1578 | **1536** | 1621 | 1699 | 1803 | 1940 | 1969 | 1993 | 1991 | 1967 | 1958 | 1947 |
| ΔBIC | — | −3127 | −1958 | −1000 | −904 | −274 | −273 | −161 | −131 | **−43** | +86 | +78 | +103 | +137 | +29 | +25 | −2 | −24 | −9 | −11 |

**Minimum at K = 10**, and it is a real interior minimum: BIC climbs 458 points over the next
six df and never returns — the best value in the whole tail is 1947 at K = 20, still 411
worse than K = 10.

#### Read the tail with care — the fit itself gets worse past K = 10

The BIC rise above K = 10 is **not** the parameter penalty doing its job. The raw deviance
worsens too, and by far more than the penalty could explain:

| K | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| deviance | **1457** | 1536 | 1607 | 1703 | 1834 | 1855 | 1873 | 1864 | 1833 | 1817 | 1799 |

**A model with 20 df fits worse than one with 10** — 417 deviance points worse at K = 16 —
and then partially recovers. Two things are going on, and only one of them is benign:

- **`ns(df = K)` bases are not nested across K.** The K − 1 interior knots sit at quantiles,
  so changing K moves *every* knot; the K = 10 basis is not a subspace of the K = 11 basis.
  Deviance is therefore not *required* to fall monotonically, unlike in a properly nested
  sequence.
- **But 417 points is far too large for that to be the whole story.** At K = 16 the model is
  estimating 17 diagonal random-effect variances from 31 shoulders. `converged = TRUE` means
  lcmm's own criteria were satisfied, **not** that the global maximum was found, and these
  fits are almost certainly sitting in worse local optima.

So: **do not over-read the shape of the curve above K ≈ 12.** What survives is the
comparison that matters — every K from 1 to 20 was fitted identically, and K = 10 gave both
the best likelihood and the best BIC of all of them. That is enough to select it. The claim
"BIC rises monotonically after 10" is not, and is not made here.

**The choice does not hinge on what *n* means.** `lcmm` computes BIC on the number of
*subjects*, so an extra spline df costs only 2 log(31) = 6.87. Recomputed on the 2,277
in-vivo observations the penalty is 2 log(2277) = 15.46 — more than twice as steep — and the
minimum is **still K = 10** ([`08_bic_by_k_extended.csv`](08_bic_by_k_extended.csv),
left panel of [`08_bic_two_penalties.png`](08_bic_two_penalties.png)). The two curves are
almost indistinguishable at this scale, which is itself the point: the deviance differences
dwarf either penalty.

#### What K = 10 actually buys — the per-shoulder view

This is **not** a second criterion competing for the answer. K = 10 is the answer. But BIC is
one number about the whole fit, and it is worth knowing how the flexibility it bought gets
distributed.

| K | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 | 15 | 20 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **max** RMSE (°) | 4.25 | 1.52 | 1.11 | 0.97 | 0.77 | 0.71 | 0.62 | 0.60 | 0.61 | 0.59 | 0.45 |
| **median** RMSE (°) | 1.06 | 0.59 | 0.36 | 0.29 | 0.22 | 0.17 | 0.15 | 0.12 | 0.11 | 0.10 | 0.075 |
| **max / median** | 4.0 | 2.6 | 3.1 | 3.3 | 3.5 | 4.3 | 4.2 | 4.9 | 5.4 | 5.7 | 6.0 |

Max RMSE plateaus at ~0.60° from K = 8 all the way to K = 15 while the median keeps falling.
The **ratio rises monotonically, 4.0 → 6.0**: every df past 6 makes the worst shoulder
*relatively* worse off. The added flexibility is going to shoulders that were already fitted
well — which is exactly what BIC, a single whole-fit number, cannot distinguish from "fits
everyone better".

#### What "K = 10" is really a statement about

Because the random part **mirrors K**, every added degree of freedom buys flexibility
*twice* — once in the shared curve, once again in each shoulder's own. BIC charges for K+1
fixed parameters and K+1 variances, but the likelihood gain from 31 shoulders each fitting
an 11-parameter curve dwarfs that penalty.

The consequence is in [`04_observed_vs_predicted.png`](04_observed_vs_predicted.png):

| prediction | concordance R² |
| --- | ---: |
| marginal — population curve only, **b**ᵢ = 0 | **0.796** |
| subject-specific — + that shoulder's BLUP | **1.000** |

Residual sd is a fraction of a degree on data spanning 73.7°. **The random curves are
reproducing the in-vivo data almost exactly.** Same over-absorption that gave
[iteration 04](../04_natural_spline_hlme/) excess kurtosis 12.55, seen from the other side.

And the tell is the marginal R², which across **all twenty K spans 0.7950 to 0.7969** — it
does not move. Nineteen added degrees of freedom change the population curve's explanatory
power by two parts in a thousand. **The population curve is not what any of them are
buying.**

### Step 1a — the individual-fit check

At the selected K = 10:

- **worst shoulder / median RMSE = 4.93×**
- **1 of 31 flagged** at robust z > 3: **`Henninger et al._19.0`** (z = 4.01, RMSE 0.60°
  against a median of 0.12°, max absolute error 1.31°, 73 points)
- marginal → subject-specific R² gap of **0.204**

Restricting step 1 to in vivo is what exposed that shoulder: pooled with the 13 ex-vivo
shoulders it sat inside the spread and nothing was flagged. Its **marginal** R² is 0.686
against 0.999 subject-specific — the population curve fits it worst of any in-vivo shoulder,
and only its own BLUP rescues it.

**The flag is K-dependent, and that is worth knowing.** The max/median ratio climbs steadily
across the grid — 4.2× at K = 8, 4.9× at K = 10, 6.0× at K = 20
([`00_bic_by_k.csv`](00_bic_by_k.csv), `rmse_ratio`). At K = 8 nothing is flagged; at K = 10
this shoulder is. Whether a shoulder looks like an outlier is therefore partly a statement
about how unevenly the added flexibility gets distributed, not a fixed property of the data.

### The experiment that proves it: hold the random part fixed

The claim above — that K is really choosing the random dimension — was tested rather than
asserted. Step 1 was re-run over K = 1…20 with the random part **decoupled** from K and held
at a fixed size ([`../../analysis_random_dim.R`](../../analysis_random_dim.R) →
[`09_decoupled_random.csv`](09_decoupled_random.csv),
[`09_what_the_df_buys.png`](09_what_the_df_buys.png)):

| random part | deviance at low K → K = 20 | gain from 19 extra fixed df | BIC minimum | max per-shoulder RMSE |
| --- | ---: | ---: | ---: | ---: |
| **mirrors K** (K+1 effects) | 9390 → 1799 | **−7591** | sharp, K = 10 | 4.25° → 0.45° |
| fixed, **4** random effects | 4291 → 4055 | −236 | flat basin K ≈ 8–20 | 1.06° → 1.04° |
| fixed, **2** random effects | 9390 → 9372 | **−18** | K = 1 | 4.25° → 4.20° |

Two rows of that table are the same models as the mirrored sweep and reproduce it exactly —
R = 1 at K = 1 gives BIC 9407.5 against the mirrored 9407.45, and R = 3 at K = 3 gives
deviance 4291.3 against 4291.29 — so the comparison is wired correctly.

**With the random part held at 2 effects, nineteen extra fixed degrees of freedom buy 18
deviance points: 0.24% of what the mirrored sweep gains.** With 4 effects they buy 236, and
BIC's basin is flat to within 10 points from K = 8 to K = 20 while max per-shoulder RMSE does
not move at all (1.06° → 1.04°).

So the steep BIC curve is almost entirely the **random** dimension. Changing the number of
random effects moves the deviance by thousands (9372 at 2 effects, 4055 at 4, 1799 at 11);
changing the fixed df at a fixed random dimension moves it by tens.

**This does not invalidate K = 10.** The model's random part mirrors K by construction — that
is iteration 04's model, kept deliberately — so K = 10 is the correct answer *for the model
being fitted*, and BIC is scoring it honestly. What the experiment settles is the
**interpretation**: K is a joint statement about how flexible the population curve *and* each
shoulder's own curve are allowed to be, and the second half dominates. Anyone reading "the
spline df is 10" as "the mean trajectory needs 10 df" would be wrong — the mean trajectory's
explanatory power is 0.795 at K = 1 and 0.796 at K = 20.

If a future iteration wants K to mean *the df of the common trajectory* — which is closer to
the literal reading of Melanie's step 1 — the random part must be fixed independently. That
is a different model, and it belongs to its own iteration, not to a second K here.

### None of this weakens the test

Steps 2 and 3 compare two models at the *same* K on the *same* rows, so the LRT is valid at
whatever K. And empirically it does not care: it rejects at every K from 1 to 10, at
p ≤ 4 × 10⁻⁶ throughout. **The K question is about the picture, not the p-value.**

### Cross-check against iteration 04

Iteration 04's joint Wald χ² on the difference block gave p ≈ 10⁻¹² for this same cell at
K = 5. The LRT here gives 1.9 × 10⁻¹⁷ at K = 10 — different statistic, different df,
different reference level, and a df chosen from the in-vivo data alone — but the same
verdict, comfortably. The two approaches agree that the conditions differ; they disagree
only about how to say it.

### Residuals

At K = 10: skew **0.210**, excess kurtosis **12.70**, Shapiro p = 2.7 × 10⁻⁴⁸ (subsampled).

The normality target Melanie set is **not met**, and this is the same trade-off documented
in [iteration 04](../04_natural_spline_hlme/): the random curves absorb so much that what is
left is a narrow core with heavy tails. Note the direction — kurtosis at K = 10 (12.70) is
*worse* than the K = 8 run's 7.15. **More df make the residuals less normal, not more**,
because each added df pulls more structure into the fit and leaves a spikier remainder. This
is unresolved and is what [`../05_bootstrap/`](../05_bootstrap/) exists for; it does not
affect the LRT, which is a likelihood comparison rather than a Wald statistic leaning on
normal standard errors.
