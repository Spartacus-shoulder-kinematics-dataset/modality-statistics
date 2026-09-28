# Generalized sweep v3 — per-cell spline df, likelihood-ratio test, Bonferroni

[Iteration 06](../06_natural_spline_hlme/)'s inference applied to all 72
joint × humeral motion × degree-of-freedom cells.

```bash
python3 prepare_monolix_data.py

# the sweep, split by step. Step 1 is the long pole — shard it.
for s in 0 1 2 3 4; do V3_NSHARD=5 V3_SHARD=$s V3_NPROC=4 Rscript analyze_all_v3_step1.R & done; wait
for s in 0 1 2 3 4; do V3_NSHARD=5 V3_SHARD=$s V3_NPROC=4 Rscript analyze_all_v3_step2.R & done; wait
Rscript analyze_all_v3_step3.R        # LRT + Bonferroni + the master CSVs
Rscript analyze_all_v3_figures.R      # seconds, as often as you like
```

`00_SUMMARY.md` carries the generated results table and is rewritten by every
figures run. This file is the method.

---

## Why this iteration exists

v2 fixed the spline df at 4 for every cell and reported **Wald tests per
coefficient**. Mélanie Prague's plan
([docs/meeting-questions.md](../../docs/meeting-questions.md)) replaces that:

> 1) Best BIC donc le degré de spline change **en fonction du marqueur**, dans un
>    modèle réduit sans in vivo departure. (a) vérification model predict vs data …
>    BLUP. (b) choix des knots au quantile.
> 2) Même degré + in vivo departure.
> 3) Test du rapport de vraisemblance, **pour chaque marqueur** — ajustement
>    multiple ici, sur la planche au complet. **Prendre Bonferroni**.
> 4) **PAS de test de Wald, on abandonne !**

[Iteration 06](../06_natural_spline_hlme/) proved that recipe on one cell. v3 is
the generalization it deferred: Bonferroni across the plate, a df per marker, and
a figure per cell.

| | v2 | **v3** |
| --- | --- | --- |
| spline df | fixed `ns(df = 4)` everywhere | **chosen per cell by BIC over K = 1…20** |
| what BIC scores | the full model | **the reduced model, on the reference condition's rows** |
| reference condition | ex vivo | **in vivo** — so the departure is the *ex-vivo* departure |
| difference curve | in vivo − ex vivo | **ex vivo − in vivo**, mirrored |
| global test | joint Wald χ² on the γ block | **likelihood-ratio test**, df = K + 1 |
| per-coefficient tests | Wald + ●/○ strip | **none** |
| multiplicity | BH-FDR | **Bonferroni** (BH kept as a secondary column) |
| per-cell figures | none | **six, one folder per cell** |
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

with $\mathbf{b}_i \sim \mathcal{N}_{K+1}(\mathbf{0}, \Omega)$, $\Omega$
**diagonal**, and $\varepsilon_{ij}\sim\mathcal{N}(0,\sigma_\varepsilon^2)$.
**K is now a property of the cell**, not a constant.

> ### ⚠ The difference curve changes sign
> v2 and iterations 03–05 plotted **in vivo − ex vivo**. Because in vivo is now the
> reference, v3 plots **ex vivo − in vivo** — the mirror image. Every axis label
> says so explicitly.

---

## The steps, and why they are separate scripts

Step 1 fits 20 models per cell; steps 2 and 3 fit one each. Splitting them means
the expensive part is sharded and cached on its own, and changing the LRT or the
multiplicity adjustment costs seconds rather than a re-sweep.

### 1 — Choose K by BIC on the reduced model, K = 1…20

Fitted on the **reference condition's own rows**. In vivo is the reference, so
the shape of the reference trajectory is settled from the reference's own data —
the ex-vivo shoulders get no vote on how many degrees of freedom the common curve
needs, since they are what will be tested against it.

```r
hlme(fixed = Y ~ ns1 + … + nsK, random = ~ ns1 + … + nsK,
     subject = "IDnum", ng = 1, idiag = TRUE, data = dref)
```

The random part **mirrors K**, as in iterations 04 and 06. Selection runs over the
**converged** rows only. Three things are recorded besides the winner, because
"BIC chose K" and "BIC ran out of grid" must not look alike:

| column | meaning |
| --- | --- |
| `n_converged` | how many of the 20 K produced a usable fit. A cell with 3–4 reference shoulders will fail at high K, and that is information |
| `bic_edge` | the minimum sat on K = 1 or K = 20 — BIC reported the **edge** rather than turning |
| `dBIC_next` | what the runner-up costs, i.e. how sharp the choice was |
| `var_per_shoulder`, `overparam` | (K + 1) variances ÷ shoulders. **At 1 or more the retained K is not identifiable** |

### The grid is 1…20 for every cell, including cells that cannot support it

This is deliberate — one rule, applied everywhere, no per-cell judgement call. But it
has a consequence that has to be read off the flags rather than assumed away.

The random part mirrors K, so a fit estimates **K + 1 diagonal variances** from however
many shoulders the cell has. And `lcmm` penalises BIC on the **subject** count, so in a
thin cell an extra parameter costs only 2 log(n) — with 9 shoulders that is 4.4, which
is no resistance at all. BIC will happily keep buying random-effect flexibility.

Acromioclavicular / horizontal flexion / DoF 1 is the worked example: 9 shoulders,
**K\* = 19 — twenty variances from nine curves**, BIC still falling at the ceiling
(K = 18 → 19 gains 24 points), and 90 % of the cell's 2333 s spent on K ≥ 12.

That K is not wrong arithmetic; it is BIC answering the question it was asked. It just
does not mean "this trajectory needs 19 degrees of freedom". Iteration 06 established
the same point on a data-rich cell — the marginal R² moved by two parts in a thousand
across the whole grid, so **K is overwhelmingly a statement about the random
dimension**, and in thin cells it is *only* that.

So every such cell is marked, in three places: `overparam` in
`00_master_summary.csv`, a red `N var/M sh` on its `00_k_map` tile, and a red line on
its own `01_choosing_k`. **A flagged K should be read as "as much flexibility as the
optimiser would bear", not as a selected df.** The LRT is unaffected either way — it
compares two models at the same K on the same rows, and iteration 06 showed the verdict
holds at every K from 1 to 8.

A cell with `bic_edge` is marked in red on `01_choosing_k`, carries a `*` on the
forest and the significance plates, and is hatched on `00_k_map`.

For single-condition cells there is no reference to speak of, so K is chosen on
whatever condition is present — recorded in `k_basis`, since several of them are
ex-vivo-only and "in vivo" would be a lie there.

### 1a — Individual-fit check, before the condition enters

> *trouver un indicateur pour vérifier que le max d'erreur chez un individu n'est
> pas trop grand par rapport aux autres*

From the subject-specific (BLUP-based) predictions at the retained K: per shoulder
`n`, `rmse`, `max_abs`. Two summaries answer it — the **headline ratio**
$\max_i(\text{rmse}_i) / \operatorname{median}_i(\text{rmse}_i)$, and a **robust
flag**. RMSE is right-skewed, so a mean ± sd rule would be dragged by the very
outliers it is meant to catch; the flag is a MAD-based z on the log scale:

$$
z_i = \frac{\log \text{rmse}_i - \operatorname{median}_j(\log \text{rmse}_j)}
{1.4826 \cdot \operatorname{MAD}_j(\log \text{rmse}_j)}, \qquad \text{flag } z_i > 3
$$

### 2 — Full model at the **same** K

```r
d$cond <- factor(..., levels = c("in vivo", "ex vivo"))   # in vivo FIRST = reference
hlme(fixed = Y ~ (ns1 + … + nsK) * cond, random = ~ ns1 + … + nsK, …)
```

K is inherited from step 1 and **not re-selected**. The script asserts the
coefficient came back named `condex vivo` — if the reference ever failed to flip,
every sign downstream would invert silently.

### 3 — Likelihood-ratio test, then Bonferroni over the plate

H₀ is the reduced model at the same K, **refitted on all shoulders**: step 1's
model chose K but was fitted on the reference condition alone, so differencing its
log-likelihood against a full-data H₁ would be meaningless. Row and subject counts
are asserted equal before the log-likelihoods are differenced.

```r
LR <- 2 * (m1$loglik - m0$loglik)
df <- K + 1                       # gamma_0 plus K shape terms
p  <- pchisq(LR, df = df, lower.tail = FALSE)
```

Then, over **every compared cell at once**:

```r
p_bonf <- p.adjust(p_lrt, "bonferroni")   # the headline
p_bh   <- p.adjust(p_lrt, "BH")           # secondary, so v3 reads against v2
```

Bonferroni is what Mélanie asked for — the conservative choice. Every star on
every plate is drawn from `p_bonf`.

### 4 — No Wald tests

Dropped entirely: no `00_wald_tests.csv`, no ●/○ strip, no per-coefficient stars.
`00_coefficients.csv` carries **estimates and standard errors only**, so nobody is
tempted to reintroduce them. The reporting object is the **difference curve**.

---

## What carries over from v2 unchanged

The per-cell data recipe is v2's, and it is shared by all three steps — a given
cell and K must yield identical rows and an identical basis in every step, or the
LRT is not a statistic:

- thin to ≤ 100 points per shoulder;
- interior knots at the cell's **own quantiles** (Mélanie's 1b);
- the `cairo_pdf` renderer, so every figure is a high-resolution PNG **and** a
  vector PDF.

### What v3 changed: boundary knots at the **fitted** range

v2 placed the boundary knots at the *untrimmed* range, on the grounds that tying
them to the overlap made the residual tails worse. **v3 reverses that**, and the
reason is the reporting object itself.

A natural spline is constrained to be **linear beyond its boundary knots**. When
the fitted data occupies only part of the basis's support — untrimmed
−3.5…187° against a fitted 14…150°, say — the outer basis columns are
near-collinear over the range actually observed. Their random-effect variances
then blow up, and because the fixed-effect covariance is conditioned on the
random structure, the difference band inflates with them.

Measured on scapulothoracic / frontal plane elevation / DoF2, **same K = 7, same
3378 rows**, changing nothing but the boundary knots:

| | boundary (−3.5, 187) | boundary (14, 150) |
| --- | --- | --- |
| random-effect params | 63 40 46 39 43 59 **2·10²** **5·10²** | 48 13 20 23 35 46 58 1.1·10² |
| SE of the difference at x = 14 | 3.62 | **2.33** |
| 95 % CI half-width at the left edge | 7.1° | **4.6°** |

A 1.6× inflated edge band on the figure the whole analysis reports is not a
detail, and it gets worse in thin cells. v2's original measurement was made at a
**fixed df = 4**, where the mismatch is mild; it does not survive per-cell K up
to 20.

The overlap trim — not the boundary knots — is what keeps the curve off empty
space, so nothing is lost by anchoring them. It also makes v3's basis identical
to iteration 06's `ns(df = K)` up to the 0.1° rounding `knots_for()` applies to
interior knots for readable labels.

### What v3 changed: `MIN_UNITS_X = 1`

This is Mélanie's standing TODO — *garder toute la plage des abscisses même si que
deux individus*.

`support_range()` returns the widest x-interval spanned by at least
`MIN_UNITS_X` shoulders. v2 set it to **2**, so a stretch only one shoulder
reached was trimmed away. v3 sets it to **1**, which degenerates to the absolute
min/max: **the whole abscissa is kept and fitted**, and the per-condition overlap
becomes the plain intersection of the two conditions' full ranges — the same rule
iteration 06 used.

It is a real trade-off, not a free win. v2 chose 2 precisely because a lone unit
reaching far past the others drags the fitted curve across empty space —
sternoclavicular horizontal flexion has one unit (17 of 13,064 rows) reaching
100° while the other nine stop near 0°. Those stretches are now fitted, and **the
95 % band is what shows how little is holding them up**: where one shoulder
carries the curve, the band is wide, and the difference is correspondingly far
from significant.

The setting matters to K, not only to the picture. The sparse tails are exactly
where extra spline df get spent, so a wider abscissa generally wants a larger K.
On the shared cell (scapulothoracic / frontal plane elevation / DoF 2) the two
settings gave:

| | `MIN_UNITS_X = 2` | **`MIN_UNITS_X = 1`** (current) |
| --- | --- | --- |
| x-range fitted | 20.1 – 138.9 | **14.0 – 150.0** |
| rows | 3002 | **3378** |
| K chosen | 7 | see `00_SUMMARY.md` |

Both K figures above were measured under v2's old boundary-knot rule, so the `= 2`
column is history rather than a clean isolation of this setting.

At `MIN_UNITS_X = 1` the rows and range match iteration 06's exactly, so that
cell is a direct check on the sweep — see `cells/scapulothoracic_frontal_plane_elevation_dof2/`.

**`TRIM_TO_OVERLAP` is still `TRUE`.** The abscissa is kept in full *within* each
condition, but a cell is still fitted on the intersection of the two conditions'
ranges, because the difference curve is only defined where both exist.

---

## Caching, resuming and sharding

The sweep is long, so interruption is treated as normal.

- **One RDS per cell per step**, under `cache/step{1,2,3}/`, written the moment
  that cell's step finishes. A killed run loses at most the cell it was inside.
- Each step **skips a cell whose RDS already exists** and whose fingerprint
  matches. `V3_REFIT=1` ignores the caches. Changing a constant in
  `analyze_all_v3_common.R` makes every cache stale automatically — the
  fingerprint is checked, not the timestamp.
- **Sharding across processes, not `mclapply`**: `hlme(nproc =)` already forks
  internally and nesting the two makes a mess. `V3_NSHARD` / `V3_SHARD` split the
  plate; one writer per file, so no locking. Bring `V3_NPROC` down when sharding
  or the shards fight for cores.
- `V3_ONLY='<joint>|<motion>|<dof>'` runs exactly one cell.
- `V3_KS='1:3'` truncates the grid for a wiring smoke test. It is part of the
  fingerprint, so a cache built under a short grid reads as stale rather than
  passing as a real selection.

`cache/` is gitignored — it is derived from `spartacus_angles_long.csv`, which is
itself untracked.

---

## Figures

### Per cell — `cells/<joint>_<motion>_dof<n>/`

| file | content |
| --- | --- |
| `01_choosing_k` | **BIC against K = 1…20** beside per-shoulder max/median RMSE — the df selection, shown rather than asserted. Log y, because at K = 1 the reduced model is a straight line whose BIC dwarfs the rest of the grid and flattens the minimum into the axis |
| `02_population_curves` | both conditions at K\*, 95 % bands over the raw scatter |
| `03_difference` | **ex vivo − in vivo** with 95 % CI, the zero line, the significant stretches shaded, and LR / df / p / Bonferroni annotated — **the headline** |
| `04_observed_vs_predicted` | marginal vs subject-specific, 1:1 line, concordance $R^2$ — the gap between the panels is what the random curves buy |
| `05_individual_fit` | per-shoulder RMSE with the robust-z threshold drawn, plus the worst few shoulders individually, so a wrong *shape* shows and not just a large average error |
| `06_diagnostics` | residuals vs fitted, Q-Q, histogram, residuals vs elevation |

Concordance $R^2 = 1 - \frac{\sum (y-\hat{y})^2}{\sum (y-\bar{y})^2}$ against the
1:1 line, **not** a squared correlation: a shoulder fitted with a constant offset
can correlate at 0.99 while sitting well off the diagonal, and only the residual
form catches it.

### Plate level

| file | content |
| --- | --- |
| `00_k_map` | **how K was chosen, for the whole plate at once** — one tile per cell carrying its K, its own BIC-vs-K curve as a sparkline with the minimum marked, and a red flag when the minimum sat on a grid endpoint or K was lost to non-convergence |
| `00_forest` | every compared cell on one difference axis: pale band = the range the difference covers across elevation, dark bar = the 95 % CI of the *mean*, plus K and the Bonferroni stars |
| `00_significance_map` | **where** in the movement each cell differs |
| `00_sigmap_elevation`, `00_sigmap_rotation` | the same, with the nesting inverted — rows are joint × DoF and the motions are the inner axis, so the planes sit adjacent on one shared abscissa |
| `planche_<motion>` | population curves, rows = joints, cols = DoF |
| `planche_diff_<motion>` | the ex vivo − in vivo difference, same layout |

## Outputs

| file | contents |
| --- | --- |
| `00_master_summary.csv` | one row per cell: counts, x-ranges, `K_star`, `k_basis`, `bic_edge`, `n_converged`, `dBIC_next`, knots, `LR`, `df`, `p_lrt`, `p_bonf`, `p_bh`, the difference block, `sig_frac_x`, $R^2$, RMSE, `n_flagged`, residual moments, seconds |
| `00_bic_by_k.csv` | **the selection evidence** — every cell × every K |
| `00_individual_fit.csv` | every cell × every shoulder |
| `00_coefficients.csv` | both models' estimates and SEs — **no Wald column** |
| `00_SUMMARY.md` | generated index and results table |

---

## Open item

**`TRIM_TO_OVERLAP <- TRUE`.** With `MIN_UNITS_X = 1` the abscissa is now kept in
full within each condition, but a cell is still fitted on the **intersection** of
the two conditions' ranges. Fitting the union instead would mean extrapolating the
difference into elevations where one condition has no data at all, which is why it
is left on. It is a single constant in `analyze_all_v3_common.R`; flipping it means
another full sweep, so it is a deliberate decision rather than a silent one.

## Standing caveat

**Condition is perfectly confounded with source study** — no study measured both.
Whatever the LRT says, the result is descriptive and not causal. See
[`../../docs/STATISTICAL_METHODS_REVIEW.md`](../../docs/STATISTICAL_METHODS_REVIEW.md).
