# Generalized sweep v2 — iteration 04's model over all 72 cells

Produced by [`../../analyze_all_v2.R`](../../analyze_all_v2.R) in ~10 minutes.

```bash
python3 prepare_monolix_data.py
Rscript analyze_all_v2.R
```

[`00_SUMMARY.md`](00_SUMMARY.md) carries the model equation, the full glyph key and the
per-cell results table. This file is the orientation: what changed from v1, what the three
summary plates are for, and what the sweep actually found.

---

## What changed from v1

[`../generalized/`](../generalized/) applied the mgcv GAMM of
[iteration 02](../02_spline_mixed_model/). This applies the random-spline `hlme` model of
[iteration 04](../04_natural_spline_hlme/).

| | v1 (`analyze_all.R`) | v2 (this) |
| --- | --- | --- |
| fitter | `mgcv::bam`, fREML | `lcmm::hlme`, full ML |
| basis | penalised thin-plate, k = 10 | natural spline, `ns(df = 4)` — 3 interior knots |
| per cell | 64 coefficients | **10 fixed effects** |
| shoulder effect | random **intercept** | random **curve** (`~ns(...)`), D **diagonal** |
| correlation | AR(1), ρ fixed | absorbed by the random curve |
| x-range | full | **trimmed to the x-overlap**, per cell |
| inference | F-test per smooth | **Wald per coefficient + a joint χ² on the difference block** |
| runtime | ~4 s/cell | ~10 s/cell |

Two rules are applied per cell and both matter:

- **Knots** come from that cell's own quantiles; **boundary knots from its untrimmed range**.
  Tying boundary knots to the overlap was tested and makes the residual tails worse
  (kurtosis 12.55 → 14.84 in iteration 04), so it is deliberately not done.
- **Curves are drawn only where at least 2 shoulders have data.** Sternoclavicular
  horizontal flexion is the case that forced this: one unit (`Oki et al._`, 17 of 13,064
  rows) reaches 100° while the other nine stop near 0°, and the absolute range was dragging
  the fitted curve across ~90° of empty space. This is outstanding item 3 of the
  [review](../../docs/STATISTICAL_METHODS_REVIEW.md) ("per-bin support rule"), and applying
  it moved the significant count from 27 to **25 of 33** — signal that lived in
  one-shoulder territory, correctly removed.

## Results

**33 compare · 30 single · 9 skipped · 0 convergence failures.**
**25 of 33** reach BH-FDR p < 0.05 on the joint difference test.

### The three summary plates

| plate | x-axis | answers |
| --- | --- | --- |
| [`00_forest.png`](00_forest.png) | degrees of difference | **how much** does each cell differ, and how certain |
| [`00_significance_map.png`](00_significance_map.png) | thoracohumeral elevation | **where in the movement** it differs, and over what fraction of the range |
| `planche_<motion>` / `planche_diff_<motion>` | elevation | the curves and the difference themselves, 4 joints × 3 DoF |

On the forest, the **wide pale band** is the range the difference covers across elevation
and the **narrow dark bar** is the 95% CI of the mean — computed properly as
$\bar{x}^\top V \bar{x}$ on the averaged design row, not by averaging pointwise SEs, which
would badly understate it.

### What the significance map shows that the forest does not

**A large mean difference does not imply a difference you can locate.** The `% of range`
column (`sig_frac_x` in the master CSV) is the fraction of the fitted range where the
95% band excludes zero, and it
dissociates sharply from effect size:

| cell | mean \|Δ\| | % of range significant | joint p (FDR) |
| --- | ---: | ---: | ---: |
| Glenohumeral · DoF1 — plane of elevation (IER 0°) | 48.9° | **0%** | 0.62 |
| Glenohumeral · DoF3 — int/ext rotation (IER 0°) | 80.3° | **12%** | 0.40 |
| Glenohumeral · DoF3 — int/ext rotation (frontal) | 60.5° | **100%** | 6.5e-06 |
| Acromioclavicular · DoF2 — medial/lateral rot (sagittal) | 20.9° | **100%** | 2.5e-33 |
| Scapulothoracic · DoF3 — post/ant tilt (sagittal) | 2.3° | 9% | 2.4e-14 |

Nine cells show **no significant region at all** — a grey bar with no blue — including two
of the largest effects on the whole sweep. Those two are also the ones the joint test
declines to call significant, so the plates agree; but a reader scanning the forest for big
numbers would land on exactly the cells with nothing behind them. **Read the two plates
together.**

The reverse pattern also occurs and is worth naming: **scapulothoracic DoF3 in sagittal
elevation is highly significant (joint p = 2.4e-14) with a mean difference of 2.3° that is
separable from zero over only 9% of the range.** Statistically solid, biomechanically
negligible — the joint test is answering "is it exactly zero", which with this much data it
almost never is.

### The replication problem, quantified

Of the 9 cells with a mean difference above 15°, **8 have exactly one ex-vivo study**
(`st_ex = 1`); the ninth has two. Across all 33 compared cells, 18 have `st_ex = 1`, and the
Spearman correlation between effect size and total study count is **−0.42**.

**The bigger the reported difference, the thinner the evidence behind it.** This is the
confounding the review's first point is about, made visible: condition is a property of the
study, so a cell with one ex-vivo study has one ex-vivo *observation* in the only sense that
matters, however many shoulders it contributes.

### Residual normality varies by cell

Excess kurtosis across the 33 compared cells: median **6.9**, quartiles 2.7–13.5, range
0.6–51.5. Only 5 cells are close to normal (< 2); 12 exceed 10. Iteration 04's residual
pathology is therefore the common case across the dataset, not a quirk of the worked
example — so the Wald standard errors here rest on an assumption most cells reject, and the
counts are descriptive. See [`../05_bootstrap/`](../05_bootstrap/) for the intended fix.

## How to read the plates

The full glyph key is in [`00_SUMMARY.md`](00_SUMMARY.md) and repeated in the footer of every
planche. The two things most easily misread:

- **The ●/○ strip is indexed.** It has 5 dots read left to right as γ₀, γ₁, …, γ₄, with the
  index printed under each dot on the figure. **γ₀ is the constant level shift**
  (`condin vivo`); γ₁…γ₄ are the shape terms `ns1:condin vivo` … `ns4:condin vivo`. Any
  strip can be checked coefficient by coefficient in `00_wald_tests.csv`.
- **Stars are the joint test, not the strip.** Because the basis coefficients are
  correlated, a cell can read `○○○○○` and still be strongly significant jointly. Trust the
  stars; the strip is texture.

## Files

| file | what it is |
| --- | --- |
| `00_SUMMARY.md` | the model equation, glyph key, and per-cell results table |
| `00_master_summary.csv` | one row per cell — mode, ranges, knots, effect sizes, `sig_frac_x` (the fraction of the range that is significant), joint and shape tests, Wald counts, kurtosis, runtime |
| `00_wald_tests.csv` | one row per **coefficient**, every cell |
| `00_forest.png/.pdf` | all compared cells on one difference axis |
| `00_significance_map.png/.pdf` | where along elevation each cell differs |
| `planche_<motion>.png/.pdf` | population curves, 4 joints × 3 DoF, shared abscissa |
| `planche_diff_<motion>.png/.pdf` | the in-vivo − ex-vivo difference, same layout |

Every figure is written as both a high-resolution PNG (2.2× scale) and a vector PDF.

## Caveat

Condition is perfectly confounded with source study — no study measured both in vivo and ex
vivo. Everything here is **descriptive, not causal**, and the effective replication is the
number of studies, not shoulders or rows.
