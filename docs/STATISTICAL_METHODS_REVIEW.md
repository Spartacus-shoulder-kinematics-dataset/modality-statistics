# Statistical-methods review

## Scope and overall assessment

This review addresses **statistical methodology only**. It evaluates the
analysis described in `article/material_and_method.tex` and implemented mainly
in `analysis_shoulder.R` and `analyze_all.R`; it does not assess the clinical or
biomechanical rationale.

The study has several strong foundations. It recognises that the observations
are repeated within shoulder, avoids forcing irregular trajectories onto an
artificial common grid, avoids extrapolating the reported contrast outside the
observed overlap, and distinguishes a trajectory-shape contrast from a constant
level shift. A spline mixed model is a sensible *descriptive* starting point for
these data, and rejecting the unconstrained sigmoid is statistically justified.

However, the current analysis does **not yet support the strength of its
in-vivo-versus-ex-vivo inferential claims**. The main issue is not the spline;
it is the design. Condition is completely confounded with source study, whereas
the fitted model treats it as though it were a shoulder-level covariate. The
AR(1), random-effects, uncertainty, and multiplicity choices further make the
displayed p-values and 95% bands too optimistic. I would recommend **major
statistical revision** before interpreting the comparisons as effects of muscle
activation or reporting widespread statistical significance.

## 1. The primary issue: condition is study-level and perfectly confounded

### What is happening

In the supplied angular data, every source article belongs to one condition:
six articles are ex vivo and eleven are in vivo. There is no article containing
both conditions. Thus `in_vivo` is a property of the **study/protocol**, not an
independently varying property of a shoulder or of an observation.

In symbols, if `study(i)` denotes the source article for shoulder `i`, then
`condition_i` is determined by `study(i)`. The proposed coefficient for
condition therefore estimates a mixture of:

* in-vivo/ex-vivo status;
* imaging and measurement protocol;
* coordinate-system implementation and data-processing choices;
* participant/cadaver characteristics and inclusion criteria;
* movement execution, range, and sampling density; and
* any other systematic difference between the two sets of studies.

It cannot isolate a causal effect of “muscle activity.” A random intercept for
shoulder, `s(ID, bs = "re")`, only models variation **among shoulders**. It
does not resolve the fact that all shoulders from an article share a protocol,
nor can it create within-study evidence for a condition contrast that the data
do not contain.

### Why this matters

Adding more rows from a dense cadaver curve does not add independent evidence
about condition. The relevant replication for the condition contrast is closer
to the number of independent studies (6 versus 11), not 809,241 rows, and in
some joint–motion–DoF cells it is smaller still. A very small p-value can
therefore mostly reflect pseudo-replication: treating many within-study points
as evidence for a study-level comparison.

### Required revision

1. Reframe the estimand as an **association between condition-labelled study
   trajectories**, not a causal muscle-activation effect. This is necessary
   unless a study measured both conditions under a comparable protocol or
   strong external adjustment data are available.
2. Preserve the article identifier in the analysis data. Add a study-level
   random effect (and, where identifiable, study-specific trajectory
   deviations) in addition to the shoulder-within-study effect. This absorbs
   some heterogeneity, but it cannot eliminate perfect confounding; say so
   explicitly.
3. Use study, not raw observations, as the unit for sensitivity analyses:
   leave-one-study-out fits; study-balanced/thinned fits; and a two-stage
   meta-analytic analysis of a pre-defined curve summary. Report how much the
   result changes when each source study is omitted.
4. If the scientific claim must be causal, the appropriate solution is a
   different design: within-study comparison, matched protocols, or a
   prospective experiment. No statistical model can identify the missing
   within-study contrast.

## 2. The hierarchical model is too simple for curve-to-curve heterogeneity

The current population model is

`Y ~ condition + s(x) + s(x, by = condition) + s(shoulder, bs = "re")`.

The random term permits each shoulder to move vertically, but assumes all
shoulders in a condition have the same curve shape. That is a very strong
assumption for pooled studies with different ranges and protocols. In practice,
a shoulder can differ in slope, curvature, onset, or range—not merely in its
overall intercept. When such variation is ignored, residual dependence remains
and standard errors for the population difference can be too small.

A better descriptive GAMM would allow curve-specific departures, for example a
factor-smooth interaction for shoulder (and preferably a study-level component)
with an appropriately constrained basis. Whether it is estimable must be
checked, because some shoulders are sparse. A pragmatic alternative is a
two-stage analysis: fit each shoulder/study curve with common, pre-specified
basis functions, then model a small set of curve summaries with study-aware
weights. Either way, show that the selected random structure improves residual
diagnostics and changes neither the substantive conclusion nor its uncertainty
materially.

Importantly, `s(ID, bs = "re")` should not be described as capturing
“between-shoulder and between-study variability”: it captures only a separate
intercept for each `ID`. Study is absent from the model.

## 3. The AR(1) correction is useful as a warning, but inadequate as final inference

The effort to account for serial dependence is commendable. Nevertheless, the
implemented AR(1) has limitations that need to be stated and investigated.

* An AR(1) process assumes correlation depends on the **number of adjacent
  observations**. Here the x values are irregularly spaced: two consecutive
  points may be 0.1° apart in one curve and many degrees apart in another. A
  single correlation `rho` consequently lacks a consistent physical meaning.
  A continuous-distance correlation, such as `exp(-|x_j-x_k|/range)`, is more
  natural when x-spacing is available.
* One common `rho` is imposed on all shoulders, motions, joints, and both
  conditions, despite radically different sampling designs. It may be a
  reasonable approximation for a sensitivity analysis, not an established
  residual model.
* `rho` is estimated from residuals of an initial independent-error fit and
  then held fixed. Its uncertainty is not propagated into the reported standard
  errors, p-values, or confidence bands. The estimator in the scripts also
  computes a pooled lag-one residual product before applying `AR.start`; it
  should calculate within-shoulder correlations only.
* Very large `rho` (about 0.996 in the exemplar) signals that hundreds or
  thousands of rows on a smooth reconstructed curve are not hundreds or
  thousands of independent measurements. It is a reason to base inference on
  shoulders/studies, not evidence of extraordinary precision.

Use residual ACFs and residuals-versus-x plots **by condition and by study**,
not only a pooled plot. Compare conclusions under plausible dependence models:
no AR term, the present AR(1), study-/curve-balanced subsampling, and a model
with a continuous residual correlation if computationally feasible. Report this
as a sensitivity analysis rather than calling the AR(1) results “honest
inference.”

## 4. The tests and intervals are approximate and currently over-interpreted

The smooth-term F-test and parametric Wald test produced by `mgcv` are
approximate, conditional on the estimated smoothing parameters and on the
plug-in `rho`. They are not exact confirmatory tests in a small, unbalanced,
study-confounded hierarchy. The confidence intervals calculated from
`vcov(model)` inherit the same conditional assumptions.

There are two further interpretation problems:

* The plotted bands are **pointwise** 95% intervals. They do not provide 95%
  simultaneous coverage across the full x-range. Reading “the curve differs
  wherever the band excludes zero” makes many location-by-location claims and
  is not justified by pointwise bands.
* The model and smoothing parameters were selected from the same data, after a
  failed sigmoid attempt. The final intervals do not reflect model-selection
  uncertainty.

Prefer uncertainty intervals for a small number of pre-specified, interpretable
estimands—for example mean difference over a clinically justified range and
difference at selected elevations. Obtain them with a study-level bootstrap or
leave-study-out resampling, which re-fits the complete procedure and naturally
captures clustering and smoothing selection better than a row-level bootstrap.
If plotting a whole difference curve, label pointwise bands as such and add a
simultaneous band only if a valid method is implemented.

## 5. Positivity/support is more than range overlap

Restricting the numerical effect-size integral to the intersection of the two
conditions’ x-ranges is an important improvement over extrapolation. Yet range
overlap alone is insufficient. There may be only a few in-vivo points in a
portion of that interval while ex-vivo contributes a dense curve there. The
estimated difference then relies primarily on smoothing assumptions rather
than a direct comparison.

For every reported comparison, show the **number of shoulders and studies with
data in x bins** for each condition, not just the overall x range and row count.
Define a common-support rule before analysis (for example, retain bins only
when each condition has observations from a minimum number of independent
studies). Calculate the effect size only on that supported region and report
the region itself. Very narrow or poorly supported overlaps should be presented
as descriptive only, or not compared.

## 6. The chosen effect size needs an estimand, uncertainty, and context

The reported quantity, mean absolute fitted difference over the overlap, is a
useful descriptive magnitude. It has advantages: it is in degrees and does not
cancel positive and negative differences. But it is not yet a complete effect
estimate.

* It is an average over **angle**, not an average over observations or people.
  This is fine only if uniform weighting over the chosen x interval is the
  intended scientific estimand. State that choice and justify the interval.
* It is unsigned, so it hides direction changes and curve crossings. Always
  show the signed difference curve and report a signed summary as well.
* It has no confidence interval, bootstrap interval, or study-level sensitivity
  range. A ranking of cells by an uncertain quantity should not be treated as a
  ranking of evidence.
* “Maximum absolute difference” is especially unstable because it is driven by
  endpoints, smoothing choices, and sparse support. It should not be a primary
  outcome without a pre-specified interior support rule.

Pre-specify one or a few primary summaries, report their uncertainty, and add a
threshold for practical relevance if one exists (e.g., a biomechanically
meaningful angular difference). A p-value cannot answer whether a 2° or 10°
difference matters.

## 7. Multiple comparisons and selective reporting must be addressed

Thirty-three condition comparisons are tested across joints, motions, and DoF.
With a 5% threshold, even independent null tests would yield false positives;
these tests are also correlated and the same studies contribute to many cells.
The result “32 of 33 significant” should therefore not be presented as 32
separate confirmatory discoveries.

Choose one of the following paths and state it in advance:

1. **Exploratory path:** describe all p-values as exploratory, emphasise effect
   estimates and uncertainty, and avoid binary star-based conclusions.
2. **Confirmatory path:** name a limited set of primary joint–motion–DoF
   hypotheses, then control the family-wise error rate or false discovery rate
   for the remaining planned family.
3. **Hierarchical path:** fit a multivariate/hierarchical model that partially
   pools related effects across DoFs and motions, while still retaining
   study-level clustering.

The arbitrary eligibility rules (at least 3 shoulders per condition; at least
4 total) also need justification. Three shoulders in a condition—potentially
from only three studies—is not enough to validate a flexible trajectory test.
At minimum, list both shoulders and distinct source studies for every cell.

## 8. Model checking and selection need stronger evidence

The manuscript says the spline was selected because the sigmoid was
non-identifiable and had worse fit. That is a reasonable rejection of that
specific sigmoid, but it does not validate the final model. BIC values from
different approximate fits (including a fixed-rho `bam` fit) should not be the
sole model-selection evidence, especially when residual correlation and random
effects differ.

For the final model, report at least:

* effective degrees of freedom and basis-dimension checks for every smooth;
* convergence/fit warnings and the estimated random-effect variance;
* residual distributions and heteroscedasticity by condition, study, and x;
* residual dependence within shoulder after fitting;
* sensitivity to basis dimension, basis type, and smoothing approach;
* whether results change when each dense study is thinned to a comparable x
  grid or each shoulder/study receives comparable weight; and
* out-of-sample or leave-one-study-out predictive checks.

Do not describe a low BIC or a non-significant diagnostic plot as proof that the
assumptions are true. They are pieces of evidence that should be triangulated.

## 9. Data-processing and reproducibility items that affect statistics

1. The input is filtered with `unit == "rad"`, but the methods and figures
   interpret the values as degrees; the x ranges and reported effects also look
   degree-scaled. This is not a cosmetic issue. Confirm the dataset metadata,
   convert explicitly if needed, and add automated unit/range checks. Do not
   label a result in degrees until provenance is resolved.
2. Preserve and model the article identifier separately from shoulder ID. The
   current derived `ID = article_shoulder` makes it impossible for the analysis
   script to fit or audit study-level effects.
3. Make the analysis input, exclusions, numbers of unique shoulders, and numbers
   of unique studies reproducible per cell. Row counts are not sample sizes.
4. Separate exploratory model development from final evaluation where possible.
   The representative-cell development process should not silently become a
   universal confirmatory specification for all 72 cells.

## A defensible revised analysis plan

The following plan would materially strengthen the work while respecting the
available data:

1. **Define the target carefully.** Primary wording: “difference between
   pooled in-vivo-labelled and ex-vivo-labelled source-study trajectories,” not
   an isolated causal effect of muscle activation.
2. **Retain study and shoulder hierarchy.** Fit a descriptive GAMM with study
   and shoulder structure plus constrained study/shoulder curve deviations, or
   use a two-stage study-aware approach if the richer GAMM is not stable.
3. **Pre-specify supported x regions and primary curve summaries.** Use study
   counts in x bins, not merely numeric range overlap. Report signed and
   absolute summaries with explicit angle weighting.
4. **Quantify uncertainty at the study level.** Refit under leave-one-study-out
   and cluster/study bootstrap resampling. Present the resulting intervals and
   influence diagnostics.
5. **Treat the 72-cell sweep as exploratory unless a primary family is set.**
   Replace or de-emphasise significance stars; report estimates, uncertainty,
   support, and study counts. Apply a stated multiplicity strategy if making
   confirmatory claims.
6. **Perform and report sensitivity analyses.** At minimum vary dense-curve
   weighting/thinning, residual-dependence assumptions, spline basis size, and
   inclusion of each source study.

## Suggested replacement for the main statistical conclusion

> Penalised mixed spline models provide descriptive estimates of differences
> between trajectories contributed by in-vivo-labelled and ex-vivo-labelled
> source studies over regions of shared angular support. Because condition is
> completely confounded with source study and protocol, these contrasts cannot
> be interpreted as causal effects of muscle activation. Their uncertainty and
> robustness should be evaluated at the study level.

This wording preserves the genuinely useful descriptive contribution while
avoiding an inference that the available design cannot sustain.
