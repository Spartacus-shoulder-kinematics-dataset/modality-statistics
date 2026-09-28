# Generalized sweep v4 — K chosen on all the data, K ≤ 10, likelihood-ratio test, Bonferroni

[Iteration 06](../06_natural_spline_hlme/)'s inference applied to all 72
joint × humeral motion × degree-of-freedom cells.

```
1  K by BIC on the REDUCED model (no condition terms), K = 1..10,
   fitted on EVERY row (ex + in vivo).  Random part mirrors K.
   Candidates: BIC minimum, first K whose BIC gain < 5 % of the total drop (default), < 1 %.
   The K in use is selected_K in K_selection.csv, editable by hand.
1a per-shoulder RMSE, robust flag, marginal vs subject-specific R2.
2  full model at the SELECTED K, all shoulders, in vivo as the reference.
3  LR = 2(LL1 - LL0) ~ chi2(K+1), Bonferroni across every compared cell.
4  no Wald tests anywhere — the difference curve is the reporting object.
```

## What came out

| cells | 72 |
| --- | ---: |
| compared | 33 |
| single condition | 30 |
| skipped | 9 |
| tested | 0 |
| significant, Bonferroni | **0** |
| significant, BH | 0 |
| K in use, range | 3 – 10 |
| K set by hand | 0 |

**Bonferroni is the headline** — Melanie asked for the conservative adjustment,
applied to the LRT p-values of every compared cell at once. BH is reported
beside it for comparison with the earlier sweeps.

## The difference curve is the reporting object

In vivo is the reference, so every difference is **ex vivo − in vivo** — mirrored
against v2 and iterations 03–05. `03_difference` in each cell folder is the figure
that carries the finding.

## Figures

| file | what it shows |
| --- | --- |
| `00_k_map` | the spline df every cell chose, with its own BIC-vs-K curve inside each tile |
| `00_forest` | every compared cell on one difference axis |
| `00_significance_map` | where in the movement each cell differs |
| `00_sigmap_elevation`, `00_sigmap_rotation` | the same, grouped by DoF so the planes sit adjacent |
| `planche_<motion>` | population curves, rows = joints, cols = DoF |
| `planche_diff_<motion>` | the ex vivo − in vivo difference, same layout |
| `cells/<cell>/01_choosing_k` | that cell's df selection, shown rather than asserted |
| `cells/<cell>/01b_fit_cloud` | step 1's fit at the K in use: data, each shoulder's own curve, the population curve |
| `cells/<cell>/02_population_curves` … `06_diagnostics` | the per-cell drill-down |

## Compared cells, largest mean |Δ| first

| joint | motion | DoF | K | shoulders ex/in | mean \|Δ\| | LR | df | p | Bonferroni | |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |

† selected_K set by hand in K_selection.csv.

## Standing caveat

**Condition is perfectly confounded with source study** — no study measured both.
Whatever the LRT says, the result is descriptive and not causal. See
[`../../docs/STATISTICAL_METHODS_REVIEW.md`](../../docs/STATISTICAL_METHODS_REVIEW.md).
