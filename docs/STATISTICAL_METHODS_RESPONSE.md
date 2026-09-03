# Response to the statistical-methods review

We thank the reviewer. The review is largely correct; in particular the central
design point is not only valid but stronger than stated. Below we respond to each
section, mark the disposition, and note what was **implemented now** versus
**planned**. Evidence is from the exemplar re-analysis in
`review_response_analysis.R` and from the updated sweep in `analyze_all.R`.

**Legend:** ✅ accepted & implemented · ◐ accepted, partially implemented ·
📋 accepted, planned · 💬 clarification.

---

## Overall assessment — ✅ accepted

We agree the previous framing over-reached. The comparison is now presented as an
**exploratory, descriptive association between condition-labelled study
trajectories**, not a causal effect of muscle activation, throughout
`material_and_method.tex`, `figures/generalized/00_SUMMARY.md`, and the planches.

---

## 1. Condition is study-level and perfectly confounded — ✅ accepted (confirmed, and worse than described)

**Verified in the data.** No source article contains both conditions: **6 ex-vivo
studies vs 11 in-vivo studies**, disjoint. For the exemplar cell
(scapulothoracic / frontal / DoF 2) the 22,102 rows and 44 shoulders reduce to
**4 ex-vivo vs 5 in-vivo independent studies**. Across the 72-cell sweep the
support is often far thinner: **27 of 33 compared cells have ≤ 3 studies in a
condition, and 21 cells have a condition represented by a single study** — e.g.
the largest reported effects (glenohumeral rotations, |Δ|≈60–70°) come from cells
with **1 ex-vivo vs 2 in-vivo studies**. Those are study-to-study anecdotes, not
population contrasts.

**Implemented.**
- `prepare_monolix_data.py` now preserves `study` (article) as its own column.
- The model gains a **study-level random intercept** `u_{s(i)}` (Eq. 1); the text
  no longer claims the shoulder intercept captures between-study variability.
- Every cell now reports **distinct studies per condition** (`st_ex/st_in`) in the
  master table, `00_SUMMARY.md`, and the planche panels.
- The estimand is restated as a confounded association (see §2.3 of the M&M).

**Demonstrated (exemplar).** Adding the study random effect changes the level
shift only modestly (−6.42° → −6.73°) but **inflates its SE (1.49 → 1.65) and
weakens p (1.6e-5 → 4.6e-5)**; this is the best-supported cell, so the effect on
thinly-supported cells is larger.

**Planned.** Full study-level refit of all 72 cells with `u_{s(i)}`, plus the
leave-one-study-out and cluster-bootstrap intervals below, as the primary
inferential output (currently the sweep reports the shoulder-level model with
study counts and FDR; the exemplar is fully re-fit with the study term).

## 2. Hierarchical model too simple — ◐ accepted

Correct that a random *intercept* forces a common curve shape and that "study"
was absent. Eq. (1) now includes `u_{s(i)}`; the misleading sentence is removed.
We agree random *slope/shape* deviations (a factor-smooth per shoulder/study)
would be preferable and note that many shoulders are too sparse to identify them.

**Planned:** where estimable, a shoulder/study factor-smooth with a constrained
basis, checked against residual diagnostics; otherwise a **two-stage** approach
(fit each study curve on a common basis, then meta-analyse a small set of curve
summaries with study weights). We will show the chosen structure improves
residual dependence without materially changing conclusions or inflating them
spuriously.

## 3. AR(1) useful as a warning, inadequate as final inference — ◐ accepted

We agree on all four points and have fixed the estimator bug:
- **Estimator corrected.** ρ is now computed **within shoulder** (no cross-curve
  products) in both `analysis_shoulder.R` and `analyze_all.R`. For the exemplar
  this changes ρ from 0.9974 (pooled) to **0.9956** (within-shoulder) — small
  here, correct in principle.
- We accept that a single ρ over irregular x-spacing lacks a consistent physical
  meaning, that one ρ across all designs is an approximation, and that ρ's
  uncertainty is not propagated. The AR(1) is now described as a **sensitivity
  device**, not "honest inference".
- ρ≈0.996 is explicitly interpreted as evidence that thousands of rows are not
  thousands of independent measurements → inference must be study-level.

**Planned:** a continuous-distance residual correlation `exp(-|x_j-x_k|/range)`
where feasible; residual ACF and residual-vs-x plots **by condition and study**;
and reporting conclusions under {no AR, AR(1), study-balanced subsampling}.

## 4. Tests and intervals approximate and over-interpreted — ◐ accepted

Accepted. The M&M now states the F/Wald tests and `vcov`-based bands are
**approximate and conditional** on the smoothing parameters and plug-in ρ; that
the plotted bands are **pointwise, not simultaneous** (so "differs wherever the
band excludes 0" is not a calibrated claim); and that model-selection uncertainty
(post-sigmoid) is not reflected.

**Implemented:** bands relabelled pointwise; the "band excludes 0" reading is
retracted in text. **Demonstrated:** leave-one-study-out on the exemplar gives
β₁ ∈ [−8.5°, −6.1°] (full-data −6.7°), i.e. real study-level sensitivity.
**Planned:** study-level bootstrap / leave-study-out intervals for a few
pre-specified estimands, and a simultaneous band only if a valid method is
implemented.

## 5. Positivity/support is more than range overlap — ◐ accepted

Strongly agree, and the sweep now exposes the problem: many overlaps are covered
by one condition densely and the other by **one study**. Restricting the ES
integral to the range overlap is necessary but not sufficient.

**Implemented:** per-cell **study counts by condition** are reported (the crudest
support check) and the effect-size ranking is annotated as support-limited.
**Planned:** per-x-bin counts of shoulders *and studies* per condition; a
pre-declared common-support rule (retain bins with ≥ k independent studies each);
compute ES only on supported bins; present narrow/poorly-supported overlaps as
descriptive-only.

## 6. Effect size needs an estimand, uncertainty, context — ◐ accepted

- It is an average **over angle** with uniform weight; this is now stated as the
  chosen estimand, and `|𝒳|` is defined as the interval width (a per-degree
  average), not a norm.
- **Signed** mean difference is now reported alongside the absolute ES (the
  master table and planches show signed Δ), so crossings are visible.
- **Max |Δ|** is retained only as a secondary, flagged as endpoint/support
  sensitive.
- ES currently lacks a CI. **Planned:** study-level bootstrap interval for one or
  two pre-specified summaries (mean signed difference over a justified range;
  difference at selected elevations), plus a practical-relevance threshold if a
  biomechanically meaningful angular difference can be agreed.

## 7. Multiple comparisons and selective reporting — ◐ accepted

Accepted; "32 of 33 significant" is not 32 discoveries.
- **Implemented:** p-values are **Benjamini–Hochberg FDR-adjusted** across the
  compared family, stars on planches are labelled **exploratory**, and the
  headline emphasises effect size + study counts over stars.
- We chose the **exploratory path** for the sweep and say so explicitly.
- **Eligibility thresholds** (≥3 per condition, ≥4 total) are now described as
  pragmatic minima, not adequacy guarantees, and both shoulder and study counts
  are listed per cell.

**Planned (if any confirmatory claim is made):** a small pre-specified primary
family of joint–motion–DoF hypotheses with study-level FWER/FDR control, or a
hierarchical model partially pooling related DoF/motion effects.

## 8. Model checking and selection — 📋 accepted, planned

Agreed that rejecting the sigmoid and a low BIC do not validate the final model.
**Planned** for the final model: per-smooth edf and `k`-index basis checks;
random-effect variances and convergence warnings; residual distribution and
heteroscedasticity by condition/study/x; within-shoulder residual dependence
after fitting; sensitivity to basis dimension/type; results under dense-study
thinning / equal study weighting; and leave-one-study-out predictive checks. We
will triangulate rather than cite a single diagnostic as proof.

## 9. Data-processing and reproducibility — ◐ accepted

1. **Unit `rad` vs degrees — 💬 open.** Confirmed the discrepancy: values are
   filtered on `unit == "rad"` but ranges/magnitudes are degree-scaled. We have
   **not** resolved provenance and will not label results in degrees until the
   dataset metadata is confirmed; an automated unit/range check will be added.
   (Flagged to the data owner.)
2. **Article preserved — ✅ done** (`study` column; see §1).
3. **Row counts ≠ sample sizes — ✅ done**: shoulder and study counts per cell are
   now first-class outputs; the narrative leads with them.
4. **Exploratory vs confirmatory separation — ◐**: the representative-cell
   development is now explicitly labelled, and the 72-cell sweep is reported as
   exploratory rather than a universal confirmatory specification.

---

## On the revised analysis plan and suggested conclusion — ✅ adopted

We adopt the reviewer's framing. The M&M conclusion/estimand now reads, in
substance:

> Penalised mixed spline models provide **descriptive** estimates of differences
> between trajectories contributed by in-vivo- and ex-vivo-**labelled source
> studies** over regions of shared angular support. Because condition is
> completely confounded with source study and protocol, these contrasts **cannot
> be interpreted as causal effects of muscle activation**; their uncertainty and
> robustness are evaluated at the **study level**.

## Summary of changes already in the repository

| Area | Change | File(s) |
|---|---|---|
| Study identity | preserved as its own column | `prepare_monolix_data.py` |
| Study RE | added `u_{s(i)}` to Eq. (1); fixed b_i wording | `article/material_and_method.tex` |
| ρ estimator | within-shoulder (bug fix) | `analysis_shoulder.R`, `analyze_all.R` |
| Study counts | `st_ex/st_in` per cell (table, summary, planches) | `analyze_all.R` |
| Signed ES | `diff_mean_signed` reported with `|Δ|` | `analyze_all.R` |
| Multiplicity | BH-FDR across compared cells; stars marked exploratory | `analyze_all.R` |
| Framing | confounding / pointwise / exploratory caveats | `article/material_and_method.tex`, `figures/generalized/00_SUMMARY.md` |
| Evidence | exemplar study-RE + LOSO + within-ρ demo | `review_response_analysis.R` |

## Key outstanding items (planned, not yet done)

1. Refit **all 72 cells** with the study random effect as the primary model.
2. **Study-level bootstrap / leave-study-out** intervals for pre-specified
   estimands (the exemplar LOSO is a proof of concept).
3. **Per-bin study support** and a common-support rule; drop/flag single-study
   overlaps.
4. **Richer random structure** (shoulder/study shape deviations) or a two-stage
   meta-analysis.
5. **Continuous-distance** residual correlation; diagnostics by condition/study.
6. **Resolve the rad/degree provenance** before reporting any result in degrees.
