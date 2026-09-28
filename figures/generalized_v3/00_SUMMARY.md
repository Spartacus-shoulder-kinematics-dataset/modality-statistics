# Generalized sweep v3 — per-cell spline df, likelihood-ratio test, Bonferroni

[Iteration 06](../06_natural_spline_hlme/)'s inference applied to all 72
joint × humeral motion × degree-of-freedom cells.

```
1  K by BIC on the REDUCED model (no condition terms), K = 1..20,
   fitted on the REFERENCE condition's own rows.  Random part mirrors K.
1a per-shoulder RMSE, robust flag, marginal vs subject-specific R2.
2  full model at the SAME K, all shoulders, in vivo as the reference.
3  LR = 2(LL1 - LL0) ~ chi2(K+1), Bonferroni across every compared cell.
4  no Wald tests anywhere — the difference curve is the reporting object.
```

## What came out

| cells | 72 |
| --- | ---: |
| compared | 33 |
| single condition | 30 |
| skipped | 9 |
| tested | 33 |
| significant, Bonferroni | **23** |
| significant, BH | 27 |
| K retained, range | 3 – 19 |
| cells whose BIC minimum sat on a grid edge | 0 |

**Bonferroni is the headline** — Melanie asked for the conservative adjustment,
applied to the LRT p-values of every compared cell at once. BH is reported
beside it only so v3 can be read against v2.

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
| `cells/<cell>/02_population_curves` … `06_diagnostics` | the per-cell drill-down |

## Compared cells, largest mean |Δ| first

| joint | motion | DoF | K | shoulders ex/in | mean \|Δ\| | LR | df | p | Bonferroni | |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Glenohumeral | internal-external rotation 0 degree-abducted | 3 | 17 | 10/21 | 66.8° | 24.5 | 18 | 1.4e-01 | 1.0e+00 | ns |
| Glenohumeral | frontal plane elevation | 3 | 9 | 10/23 | 58.0° | 46.4 | 10 | 1.2e-06 | 4.0e-05 | *** |
| Glenohumeral | internal-external rotation 0 degree-abducted | 1 | 17 | 10/21 | 56.6° | 34.4 | 18 | 1.1e-02 | 3.7e-01 | ns |
| Glenohumeral | sagittal plane elevation | 3 | 9 | 10/24 | 54.9° | 416.7 | 10 | 2.6e-83 | 8.7e-82 | *** |
| Acromioclavicular | sagittal plane elevation | 2 | 11 | 10/3 | 20.8° | 299.6 | 12 | 5.6e-57 | 1.9e-55 | *** |
| Acromioclavicular | frontal plane elevation | 2 | 3 | 10/4 | 19.9° | 34.5 | 4 | 5.9e-07 | 1.9e-05 | *** |
| Scapulothoracic | horizontal flexion | 2 | 9 | 10/7 | 19.7° | 104.1 | 10 | 8.0e-18 | 2.7e-16 | *** |
| Glenohumeral | internal-external rotation 0 degree-abducted | 2 | 10 | 10/21 | 19.6° | 48.2 | 11 | 1.3e-06 | 4.3e-05 | *** |
| Glenohumeral | sagittal plane elevation | 1 | 9 | 10/24 | 18.8° | 406.4 | 10 | 4.2e-81 | 1.4e-79 | *** |
| Sternoclavicular | frontal plane elevation | 3 | 11 | 13/4 | 13.0° | 179.8 | 12 | 4.7e-32 | 1.6e-30 | *** |
| Scapulothoracic | horizontal flexion | 3 | 8 | 10/7 | 12.4° | 25.2 | 9 | 2.8e-03 | 9.1e-02 | ns |
| Sternoclavicular | sagittal plane elevation | 3 | 8 | 13/3 | 11.9° | 96.1 | 9 | 9.4e-17 | 3.1e-15 | *** |
| Scapulothoracic | horizontal flexion | 1 | 9 | 10/7 | 11.8° | 135.6 | 10 | 3.3e-24 | 1.1e-22 | *** |
| Acromioclavicular | frontal plane elevation | 3 | 9 | 10/4 | 10.8° | 23.2 | 10 | 1.0e-02 | 3.3e-01 | ns |
| Scapulothoracic | frontal plane elevation | 2 | 10 | 13/31 | 10.1° | 105.0 | 11 | 1.8e-17 | 6.0e-16 | *** |
| Glenohumeral | frontal plane elevation | 2 | 10 | 10/23 | 9.8° | 122.0 | 11 | 7.1e-21 | 2.3e-19 | *** |
| Scapulothoracic | internal-external rotation 0 degree-abducted | 2 | 15 | 10/21 | 9.6° | 30.0 | 16 | 1.8e-02 | 5.9e-01 | ns |
| Glenohumeral | frontal plane elevation | 1 | 9 | 10/23 | 9.4° | 42.2 | 10 | 6.9e-06 | 2.3e-04 | *** |
| Scapulothoracic | sagittal plane elevation | 2 | 9 | 13/25 | 9.1° | 119.1 | 10 | 7.7e-21 | 2.5e-19 | *** |
| Scapulothoracic | internal-external rotation 0 degree-abducted | 1 | 11 | 10/21 | 7.6° | 175.1 | 12 | 4.3e-31 | 1.4e-29 | *** |
| Glenohumeral | sagittal plane elevation | 2 | 13 | 10/24 | 7.4° | 385.9 | 14 | 1.2e-73 | 3.8e-72 | *** |
| Acromioclavicular | sagittal plane elevation | 1 | 8 | 10/3 | 7.2° | 29.5 | 9 | 5.4e-04 | 1.8e-02 | * |
| Scapulothoracic | sagittal plane elevation | 1 | 11 | 13/25 | 7.0° | 39.4 | 12 | 8.9e-05 | 2.9e-03 | ** |
| Scapulothoracic | frontal plane elevation | 3 | 10 | 13/31 | 7.0° | 63.0 | 11 | 2.6e-09 | 8.5e-08 | *** |
| Scapulothoracic | frontal plane elevation | 1 | 11 | 13/31 | 6.7° | 16.6 | 12 | 1.6e-01 | 1.0e+00 | ns |
| Sternoclavicular | sagittal plane elevation | 2 | 8 | 13/3 | 6.5° | 58.1 | 9 | 3.1e-09 | 1.0e-07 | *** |
| Sternoclavicular | frontal plane elevation | 1 | 8 | 13/4 | 5.9° | 17.4 | 9 | 4.3e-02 | 1.0e+00 | ns |
| Scapulothoracic | internal-external rotation 0 degree-abducted | 3 | 17 | 10/21 | 5.7° | 19.8 | 18 | 3.4e-01 | 1.0e+00 | ns |
| Sternoclavicular | frontal plane elevation | 2 | 11 | 13/4 | 5.2° | 138.9 | 12 | 1.0e-23 | 3.3e-22 | *** |
| Acromioclavicular | frontal plane elevation | 1 | 8 | 10/4 | 4.4° | 26.9 | 9 | 1.5e-03 | 4.8e-02 | * |
| Sternoclavicular | sagittal plane elevation | 1 | 8 | 13/3 | 3.9° | 10.3 | 9 | 3.3e-01 | 1.0e+00 | ns |
| Acromioclavicular | sagittal plane elevation | 3 | 8 | 10/3 | 3.0° | 2.0 | 9 | 9.9e-01 | 1.0e+00 | ns |
| Scapulothoracic | sagittal plane elevation | 3 | 11 | 13/25 | 2.7° | 80.3 | 12 | 3.7e-12 | 1.2e-10 | *** |

\* the BIC minimum sat on an endpoint of the K grid — the edge was reported, not a choice.

## Standing caveat

**Condition is perfectly confounded with source study** — no study measured both.
Whatever the LRT says, the result is descriptive and not causal. See
[`../../docs/STATISTICAL_METHODS_REVIEW.md`](../../docs/STATISTICAL_METHODS_REVIEW.md).
