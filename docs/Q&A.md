# Q&A — decisions & reasoning digest

Synthetic log of what was asked, what was decided, and the alternatives — kept so
the reasoning can be revisited. Generated/updated by the `/qa` skill.

---

## 1. Which statistical method should compare in-vivo vs ex-vivo shoulder rhythm in Spartacus, and is LGCM the right tool?   _(2026-07-15)_

**A.** Use a **mixed-effects model of the trajectory**, not Latent Growth Curve
Models. The data is unbalanced, irregularly sampled, and has few units per
group, so SEM growth curves (which need balanced, aligned time points and many
units) are inappropriate. Deliverable: `METHODOLOGY.md` separating the *method*
from the *software* (MonolixSuite and R).

**Decision — modelling family:**
- ★ **(Nonlinear) mixed-effects / GAMM** — uses every raw point at its true x, models between-study variance as random effects.
- LGCM / SEM growth curves — rejected: needs wide-format aligned grid + large N.
- Descriptive polynomial on group-mean curves — rejected: discards multi-study structure; keep only for visuals.

**Decision — software:**
- ★ **R-first, Monolix optional** — free, scriptable, reproducible for the intern.
- Monolix-centric (GUI + SAMBA) — good, but licensed and less reproducible.
- Framework-agnostic — clean but not runnable now.

## 2. Run the full analysis in R (no Monolix license), pedagogically, saving all figures.   _(2026-07-15)_

**A.** Built `analysis_shoulder.R`. The parametric **sigmoid was non-identifiable**
(the rhythm never plateaus, so its asymptotes fall outside the data), so the
workhorse is a **penalised-spline GAMM** (`mgcv`) with a per-shoulder random
intercept. Residuals are strongly autocorrelated (lag-1 ≈ 0.996) → added an
**AR(1)** term. Result: in-vivo differs from ex-vivo (level shift + shape).

**Decision — structural model:**
- ★ **Penalised spline (GAMM)** — flexible, identifiable on non-saturating data.
- 4-parameter logistic NLME — rejected: asymptotes unconstrained → blows up.
- Polynomial linear mixed model — usable baseline, but rigid and its random-slope covariance was degenerate.

## 3. Keep every modelling iteration, with a README on why we moved on.   _(2026-07-15)_

**A.** Split figures into iteration sets — `00_data_exploration/`,
`01_sigmoid_nlme/`, `02_spline_mixed_model/` — plus `figures/README.md` narrating
why each step superseded the last (match the model to the data; correct for
autocorrelation).

## 4. Generalise the analysis to every joint angle, every joint, every movement.   _(2026-07-15)_

**A.** `analyze_all.R` sweeps all **72** joint×movement×DoF combinations
(33 compared, 30 single-condition, 9 skipped for too few units). Key fix: measure
the in-vivo−ex-vivo difference **only on the x-overlap** of both conditions
(full-range gave impossible ~1500° extrapolations).

**Decision — cross-combo summary metric:**
- ★ **Effect size = mean |in-vivo − ex-vivo| (degrees)** on the overlap — biomechanically meaningful and discriminating.
- −log10(difference-smooth p) — rejected: p saturates (≈0) everywhere → uniform, uninformative.
- Raw level shift only — misses shape differences where curves cross.

## 5. Make the overviews *planches* — one per motion, rows = joints, cols = DoF.   _(2026-07-15)_

**A.** Six planches matching the `shoulder-kinematics/spartacus/plots`
convention: rows = joints proximal→distal (sternoclavicular → acromioclavicular →
scapulothoracic → glenohumeral), cols = DoF 1/2/3 with anatomical angle names;
each panel overlays the in-vivo/ex-vivo fits on raw curves.

## 6. Give each planche row the same y-axis.   _(2026-07-15)_

**A.** Per joint row, a common y-range is computed from the raw data plus fitted
CI across its three DoF panels, so the rotations are comparable within a joint
while different joints keep their own scale.

## 7. Show significance stars + per-box stats; explain the lone black curve.   _(2026-07-15)_

**A.** Each compared panel now shows top-right **trajectory-difference stars**
(`*** / ** / * / ns`) and a bottom-left block: n shoulders (ex/in), mean & max
|Δ| (deg), and the **level shift with its own stars**. The lone **black curve =
single-condition pooled fit** (only in-vivo *or* ex-vivo available), labelled
per panel and in the header.

**Decision — which p drives the stars:**
- ★ **Difference-smooth p top-right + level-shift p in the block** — shows both the shape test and the discriminating offset test.
- Difference-smooth p only — rejected: saturates at *** everywhere.
- Level-shift p only — misses shape-only differences (crossing curves).

## 8. Write the Materials & Methods section for the statistical study.   _(2026-07-15)_

**A.** `article/material_and_method.tex` — dataset/selection, the GAMM (Eq. 1),
AR(1), effect size on the overlap (Eq. 2), generalisation to planches, and
software (R 4.3.3, mgcv 1.9.1). Uses `natbib`; `references.bib` stub carries the
8 citation keys.

## 9. Compile / run the LaTeX.   _(2026-07-15)_

**A.** Added `article/main.tex` preview wrapper; builds a clean 3-page PDF with
**0 undefined citations** (`pdflatex` → `bibtex` → `pdflatex` ×2). Equations and
citations render correctly.

---

## 10. Where is the spline coded, what are its coefficients, how are they fitted?   _(2026-07-15)_

**A.** One formula: `Y ~ condO + s(TIME) + s(TIME, by=condO) + s(ID, bs="re")`
([analysis_shoulder.R:143](analysis_shoulder.R#L143)). `s(TIME)` is a thin-plate
spline `f(x)=ΣβₖBₖ(x)` (default `k=10` basis functions); the **coefficients βₖ
are estimated, not set** — we only choose structure/basis/`k`. `bam()` solves a
**penalised least-squares** problem `‖y−Xβ‖² + ΣλⱼβᵀSⱼβ`; λ (smoothness) is picked
by fREML, then AR(1) is added ([:144-149](analysis_shoulder.R#L144-L149)).

**Decision — basis & k:**
- ★ **defaults (thin-plate, k=10)** — flexible, auto-penalised (k is just a ceiling).
- cubic regression spline (`bs="cr"`) — cheaper, but you manage knots.
- larger k — more wiggle, more compute; unneeded here.

## 11. What does every symbol in Eq. (1) mean?   _(2026-07-15)_

**A.** `yᵢⱼ`=angle of point `j` in shoulder `i`; `β₀`=ex-vivo baseline;
`𝟙(c=in vivo)`=1 for in-vivo rows only; `β₁`=constant in-vivo level shift;
`f(x)`=ex-vivo reference trajectory (spline); `g(x)`=in-vivo shape difference (so
in-vivo = f+g+β₁); `bᵢ~N(0,σ²_b)`=per-shoulder random intercept; `εᵢⱼ`=residual
(AR(1)-correlated). Read: *ex-vivo = β₀+f+bᵢ+noise; in-vivo = that + β₁ + g(x).*

## 12. What is a GAMM and why does it handle irregularly sampled data?   _(2026-07-15)_

**A.** GAM = mean is a sum of **smooths**, each a penalised basis expansion
`f=ΣβₖBₖ(x)` with a wiggliness penalty `λ∫f''²` (λ learned). GAMM = GAM + random
effects (same math — a penalty is a Gaussian prior). It handles irregular
sampling because `f` is a function of **continuous x**, not a parameter per
timepoint: every point evaluates the basis at its own x and pools into one `f` —
**no grid, no aligned columns** (exactly what LGCM requires and this avoids).
Sparse regions borrow strength via the penalty.

## 13. Teach me fREML.   _(2026-07-15)_

**A.** REML estimates the **variance/smoothing parameters** (λ) from *residual
contrasts* with the fixed effects projected out, removing ML's downward variance
bias; since `λ=σ²_ε/σ²_f`, choosing smoothness = estimating a variance component.
It is more **stable** than GCV/AIC. **fREML** = a fast, numerically stable
reformulation (Wood 2011) used by `bam()` to scale to large `n`.

## 14. Teach me the AR(1) residual.   _(2026-07-15)_

**A.** `εₜ = ρ·εₜ₋₁ + noise`: successive residuals along a curve are correlated
(`ρ≈0.996` here). Ignoring it treats correlated points as independent →
effective N hugely overstated → SEs too small → everything "significant."
`bam(rho=ρ, AR.start=…)` discounts that redundancy; `AR.start` resets the chain
at each shoulder's first row so correlation doesn't leak across shoulders.
Approximate (assumes even spacing, constant ρ) but far better than ignoring it.

## 15. Why a Wald t-test for β₁ but an F-test for g?   _(2026-07-15)_

**A.** Different objects. `β₁` is **one scalar** → `estimate/SE ~ t` (Wald).
`g(·)` is a **whole smooth** (several basis coefficients + effective df) → a
**joint/multi-df** question, so an approximate **F-test** on its edf. A single t
cannot test a multidimensional term.

**Decision — how to test the two effects:**
- ★ **t for β₁ + F for g** — each matched to the term's dimension, native to `mgcv`.
- LRT between nested models — cleaner theory, but extra refits × 72 combos.
- one omnibus "any condition effect" test — simplest, can't separate level vs shape.

## 16. Is the effect size classical, and what is `|𝒳|` (a norm)?   _(2026-07-15)_

**A.** `ES = (1/|𝒳|)∫_𝒳 |β₁+g(x)| dx`. `𝒳` = the x-range where **both**
conditions have data; **`|𝒳|` is its length (e.g. 130°), NOT a vector norm** —
just the width used to average. `β₁+g(x)` = in-vivo − ex-vivo difference at `x`;
`|·|` = absolute value. So ES = **average absolute curve difference in degrees**
over the shared range — in code `mean(abs(dif))`. It's a **functional** effect
size in native units, *not* a standardized one (Cohen's d).

**Decision — effect-size metric:**
- ★ **mean |Δ| (degrees) on overlap** — interpretable, sign-robust.
- standardized (÷ residual SD) — unit-free, but hides physical magnitude.
- peak |Δ| (max) — flags the worst point, ignores the rest.

## 17. Why fit a descriptive single-group model when only one condition exists?   _(2026-07-15)_

**A.** With only in-vivo *or* only ex-vivo shoulders there's nothing to contrast
(can't estimate β₁ or g). We fit `y=β₀+f(x)+bᵢ+ε` (AR(1)) — mean trajectory +
random intercepts, no condition terms ([analyze_all.R:160](analyze_all.R#L160))
— so the planche still shows a fitted population curve (the lone **black** curve)
for that DoF, making **no** in-vivo/ex-vivo claim. Keeps the atlas complete.

**Decision — solo-condition combos:**
- ★ **descriptive single-group smooth** — completes the atlas, labelled "no contrast".
- skip them — cleaner, but 30 of 72 cells go blank.
- pool both labels into one curve — never: fakes a comparison that isn't there.
