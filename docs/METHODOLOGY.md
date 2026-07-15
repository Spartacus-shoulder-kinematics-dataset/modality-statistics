# Statistical Methodology — Comparing In-Vivo vs Ex-Vivo Scapulothoracic Rhythm (Spartacus)

> **Audience:** the M2 biostatistics intern.
> **Purpose:** decide *how* to test whether the scapulothoracic kinematic
> trajectory differs between in-vivo and ex-vivo conditions in the Spartacus
> dataset, and *how* to implement it. The **statistical method** (Part 1) is
> deliberately kept separate from the **software** (Parts 2 & 3): the method is
> what you must understand and defend; the software is just a tool that
> implements it.

---

## 0. The scientific question and the data

**Question.** During arm elevation, does the scapulothoracic (ST) angle follow a
different trajectory as a function of the humerothoracic (HT) elevation angle
depending on condition (in-vivo, i.e. with muscle activation, vs ex-vivo, i.e.
cadaveric / passive)?

**Feature-example subset** used throughout this document (a single well-defined
slice; the same recipe generalizes to other joints/motions/DoFs):

| filter | value |
| --- | --- |
| `joint` | `scapulothoracic` |
| `humeral_motion` | `frontal plane elevation` |
| `degree_of_freedom` | `2` |
| `unit` | `rad` |

**What the data actually looks like** (measured on `corrected_confident_data.csv`):

- **22,102 observations** across **44 `article × shoulder` units**.
- **11 to 2,987 points per unit** — extreme imbalance.
- x = `humerothoracic_angle` is **irregularly sampled and differs across
  studies**: ex-vivo cadaver rigs produce dense continuous curves (18,302 rows),
  in-vivo protocols produce sparse points (3,800 rows).
- Each row is one (unit, x, y) triple: `y = value` (ST angle),
  `x = humerothoracic_angle`, unit = `article` + `shoulder_id`.

This shape — **unbalanced, irregularly sampled, few grouping units, repeated
measures within unit** — is the single most important fact for choosing a
method.

---

## 1. The method (software-agnostic)

### 1.1 Chosen family: Nonlinear Mixed-Effects Model (NLME)

Also called a **population model** or a **nonlinear multilevel/hierarchical
model**. It has three pieces:

**(a) A structural model** — a smooth function `f` describing the *shape* of one
unit's curve:

```
y_ij = f(x_ij ; φ_i) + ε_ij
```

for observation *j* of unit *i*, where `φ_i` is that unit's vector of curve
parameters and `ε_ij` is residual noise.

Candidate shapes for the scapulohumeral rhythm (typically S-shaped: a low-slope
"setting phase" at low elevation, then a steeper near-linear rise, sometimes a
plateau near end-range):

1. **Sigmoid / 4-parameter logistic (Emax-type)** — *recommended starting
   point*. Interpretable parameters (baseline, amplitude, midpoint, steepness)
   and it maps directly onto Monolix's built-in model library.
   `f(x) = base + amplitude · x^h / (mid^h + x^h)`
2. **Segmented / bilinear with a breakpoint** — biomechanically explicit:
   a setting-phase intercept, a breakpoint elevation, and two slopes.
3. **Quadratic polynomial** — the shape used in `PolynomialFit.R`. Note this is
   *linear in its parameters*, so it is technically a **linear** mixed model
   (LMM), not nonlinear. Keep it only as a **simple baseline** to benchmark the
   nonlinear forms against.

> Do not pick the final shape dogmatically. Fit 2–3 candidates and let the
> selection criteria in §1.4 decide.

**(b) Random effects** — every unit deviates from the population-typical curve.
Put random effects on the structural parameters, at the `article × shoulder`
level:

```
φ_i = φ_pop ⊙ exp(η_i) ,   η_i ~ N(0, Ω)
```

(log-normal parameters keep them positive; this is the Monolix default). `Ω` is
the between-study/between-shoulder covariance — it explicitly models the fact
that the 44 units come from different studies.

**(c) A residual error model** — start with a **constant** (additive) error and
also try **proportional**; compare.

### 1.2 Why NLME and not LGCM (why we gave up on Latent Growth Curve Models)

LGCM / SEM growth curves (`lavaan::growth`) require the data in **wide format
with aligned timepoints**: every subject measured at the *same* set of x-values,
laid out as columns `y1, y2, …, yk`. Our data violates every precondition:

- x is **irregular and study-specific** — there is no common grid of elevation
  angles. Forcing a wide format via `trial_number` (as the old `LGM_spartacus.R`
  did) silently pretends "trial 1" of one study is the same elevation as "trial
  1" of another. It is not — this **fabricates structure**.
- With only **44 units** and a possibly large number of pseudo-timepoints, the
  SEM covariance matrix is unstable / underidentified.
- LGCM's latent intercept/slope factors assume a (piecewise-)linear growth basis,
  which is a poor description of an S-shaped rhythm.

NLME sidesteps all of this: it uses **every raw point at its true x**, absorbs
between-study variability into random effects, and needs no common grid.

### 1.3 The hypothesis test: `in_vivo` as a covariate

This is exactly the "covariate selection on the model" that Mélanie described.
The condition `in_vivo` (categorical: True/False) enters the **population**
value of one or more structural parameters:

```
φ_pop,i = φ_ref · exp( β_φ · 1[in_vivo_i = True] )
```

- `β_φ = 0` for a parameter ⇒ condition does **not** change that aspect of the
  trajectory.
- `β_φ ≠ 0` ⇒ condition shifts it. Because the parameters are interpretable, a
  significant `β` answers biomechanical questions directly, e.g. *"in-vivo
  increases the steepness of scapular rotation after mid-range"* (β on the
  steepness parameter) or *"in-vivo raises the plateau amplitude"* (β on
  amplitude).

This is the NLME equivalent of the LGCM "conditional/multi-group model", but
valid for unbalanced data.

### 1.4 Model & covariate selection

1. Fit a **base model** (no covariate).
2. Add `in_vivo` on candidate parameters → **full model(s)**.
3. Compare using:
   - **BICc** (Monolix) / **BIC** — lower is better; primary criterion.
   - **Likelihood-ratio test** between nested models — fit with **ML, not
     REML**, when the comparison concerns fixed effects.
   - **Wald test / p-value** on each `β_φ`.
4. Check the fit, don't just trust a number:
   - **VPC** (Visual Predictive Check), per group.
   - **Individual fits** for a sample of units.
   - **Residual diagnostics**: IWRES / NPDE.
   - **Shrinkage** on the random effects (high shrinkage ⇒ that random effect is
     not supported by the data).

### 1.5 Practical cautions specific to this dataset

- **Imbalance (18k ex-vivo vs 4k in-vivo rows, 11–2987 pts/unit).** Densely
  sampled units must not dominate. The per-unit random effect helps, but also
  sanity-check by (optionally) thinning ultra-dense curves or verifying results
  are stable when down-weighting them.
- **Units of x vs y.** In this subset `y` (`value`) is labelled `rad` but its
  magnitudes look like **degrees** (x ~0–180, y in the tens). Flag this for the
  data owner; keep x in its own scale and be explicit in the model and plots.
- **Small number of units (44).** Standard errors on `Ω` will be modest —
  interpret between-study variance components cautiously.
- **Within-curve autocorrelation (important).** Each unit is a dense, smooth
  trajectory, so consecutive residuals are strongly correlated (measured lag-1
  ≈ **0.996**). Treating points as independent makes every p-value look
  absurdly significant. Use an **AR(1)** residual model (or equivalent) and
  report the corrected inference.
- **The trajectory does not saturate (model-shape lesson).** Over the observed
  elevation range the ST angle keeps changing monotonically without a plateau,
  so a **sigmoid/logistic is a poor, non-identifiable choice here** (its
  asymptotes fall outside the data). A flexible spline (or low-order polynomial)
  fits far better — see Part 3.

---

## 2. How to do it in MonolixSuite (primary tool)

MonolixSuite is a GUI for exactly this class of model (nonlinear mixed effects,
SAEM estimation) with built-in covariate selection — no R expertise required.

**Step 1 — Prepare the data.** Export a long-format CSV for the feature subset
with Monolix-friendly columns. Use the provided script
[`prepare_monolix_data.py`](prepare_monolix_data.py):

```
python3 prepare_monolix_data.py
```

It writes `monolix_st_frontal_dof2.csv` with columns:

| column | role in Monolix | source |
| --- | --- | --- |
| `ID` | subject identifier | `article` + `_` + `shoulder_id` |
| `TIME` | the regressor / independent variable (x) | `humerothoracic_angle` |
| `Y` | observation | `value` |
| `in_vivo` | **categorical covariate** | `in_vivo` (True/False) |

> In Monolix "TIME" is simply the independent variable — here it is elevation
> angle, not clock time.

**Step 2 — Datxplore.** Open the CSV in Datxplore to visualize the raw curves
and *see* the sampling imbalance before modelling.

**Step 3 — Monolix project.**
1. New project → load `monolix_st_frontal_dof2.csv`.
2. Assign column types: `ID` → *id*, `TIME` → *time*, `Y` → *observation*,
   `in_vivo` → *categorical covariate*.
3. **Structural model**: pick one from the model library (e.g. a sigmoid Emax /
   `Emax` model) or write a short custom `.txt` (mlxtran) model of `f(x; φ)`.
4. **Statistical model**: set parameter distributions to **lognormal**, enable
   **random effects** on the structural parameters, choose an **error model**
   (start `constant`, also try `proportional`).

**Step 4 — Run tasks (top toolbar), in order:**
1. *Population parameters* (SAEM).
2. *Standard errors* (Fisher information).
3. *Log-likelihood* (importance sampling) → gives −2LL, AIC, **BICc**.
4. *Plots* → VPC, individual fits, residual (NPDE) plots.

**Step 5 — Covariate selection (the hypothesis test).** Two routes:
- **Manual:** add `in_vivo` on a parameter, re-run, compare **BICc** and the
  Wald p-value of `β` against the base model.
- **Automatic:** use the **covariate model building assistant (SAMBA / COSSAC)**
  to let Monolix search which parameters `in_vivo` should act on.

**Step 6 — Report.** Population parameters, the `in_vivo` covariate coefficients
with p-values, a **BICc comparison table** (base vs covariate models), and
**per-group VPCs**.

**Resources.** Lixoft channel <https://www.youtube.com/c/Lixoft> and the Monolix
documentation (<https://monolixsuite.slp-software.com/>).

---

## 3. How to do it in R (fully worked — since Monolix is not yet licensed)

There is a **complete, runnable analysis** in
[`analysis_shoulder.R`](analysis_shoulder.R). It uses only preinstalled packages
(`mgcv`, `nlme`, `dplyr`) and writes every figure to `figures/` (organised as
iteration sets — see [`figures/README.md`](figures/README.md)). Run it with:

```bash
python3 prepare_monolix_data.py   # writes monolix_st_frontal_dof2.csv
Rscript analysis_shoulder.R       # data-exploration + models + all figures
```

### 3.1 What the data forced us to choose

We first tried the **parametric sigmoid NLME** (below). On this data it is
**non-identifiable** (the curve never plateaus, so the logistic asymptotes and
their per-unit random effects blow up) and fits worst. We therefore use a
**penalised spline mixed model** — still a nonlinear mixed model — as the
workhorse:

```r
library(mgcv)
# condO = ordered factor(cond) so the by-smooth is a direct difference curve.
# AR.start marks each unit's first row; rho = estimated lag-1 autocorrelation.
m <- bam(Y ~ condO + s(TIME) + s(TIME, by = condO) + s(ID, bs = "re"),
         data = d, method = "fREML", discrete = TRUE,
         rho = 0.996, AR.start = d$ar_start)      # AR(1) for autocorrelation
summary(m)                                        # s(TIME):condO... = the test
```

- `s(TIME, by = condO)` is the in-vivo-vs-ex-vivo **shape difference** (the
  hypothesis); `condO` (parametric) is the **level shift**.
- **Result on the feature subset:** both are significant *after* the AR(1)
  correction (level shift t ≈ −4.3; shape difference p < 0.001) — in-vivo sits
  ~10–15° below ex-vivo through low-to-mid elevation. See
  `figures/02_spline_mixed_model/`.

### 3.2 The parametric NLME route (kept for the eventual Monolix work)

Still useful to understand, and closest to what Monolix does — just expect the
non-identifiability above with a sigmoid, so prefer a non-saturating structural
form there too.

**Option A — `saemix`** (closest to Monolix: same SAEM algorithm and outputs).

```r
library(saemix)
d <- read.csv("monolix_st_frontal_dof2.csv")

# structural model: sigmoid Emax  f(x) = base + amp * x^h / (mid^h + x^h)
model <- function(psi, id, x) {
  base <- psi[id, 1]; amp <- psi[id, 2]; mid <- psi[id, 3]; h <- psi[id, 4]
  tt <- x[, 1]
  base + amp * tt^h / (mid^h + tt^h)
}

sm  <- saemixData(name.data = d, name.group = "ID",
                  name.predictors = "TIME", name.response = "Y",
                  name.covariates = "in_vivo")

# base model (no covariate)
m0 <- saemixModel(model = model, psi0 = c(base=0, amp=1, mid=60, h=2),
                  transform.par = c(1,1,1,1))            # lognormal params
fit0 <- saemix(m0, sm, list(seed = 1))

# full model: in_vivo acts on the structural parameters
m1 <- saemixModel(model = model, psi0 = c(base=0, amp=1, mid=60, h=2),
                  transform.par = c(1,1,1,1),
                  covariate.model = matrix(1, nrow = 1, ncol = 4))
fit1 <- saemix(m1, sm, list(seed = 1))

# compare
compare.saemix(fit0, fit1)     # BIC / AIC / -2LL
```

Read `β` estimates and Wald p-values from `summary(fit1)`. Note: with a sigmoid
you will hit the non-identifiability described in §3.1; use a non-saturating
structural function.

**Option B — `nlme::nlme`** for a Gaussian NLME with the same structure
(`fixed = ...`, `random = ...`, covariate via the `fixed` formula), comparing
models with `anova()` / `BIC()`. The worked script fits this as the labelled
"attempt 1" comparison model.

---

## 4. What is *not* the inferential method: the descriptive plot

`PolynomialFit.R` fits polynomials to **group-mean curves**. It is useful only
for **visualization**: it aggregates away the multi-study structure and the
within-group variability, so it must **not** be used to test the in-vivo /
ex-vivo hypothesis. Use it to eyeball the shapes; use Parts 1–3 for inference.

The old LGCM scripts (`LGM_spartacus.R`, `first_comparison.R`) are **retired** by
the decision to drop LGCM.

---

## 5. Suggested work plan for the intern

**Now (R, no license needed):**
1. Run `prepare_monolix_data.py`, then `analysis_shoulder.R`; study the figures
   in order (`figures/README.md`) to understand the data and the modelling.
2. Reproduce and interrogate the spline mixed model: check the AR(1) correction,
   the diagnostics, and the difference figure. Confirm the conclusion is robust
   (e.g. try a quadratic-polynomial mixed model; down-weight the densest units).
3. For the whole shoulder at once, run **`analyze_all.R`**: it applies the same
   model to **every** joint x movement x DoF and assembles **planches** (one
   plate per movement, joints x DoF) plus a master effect-size table. See
   `figures/generalized/00_SUMMARY.md`.

**Later (when MonolixSuite is licensed):**
4. Repeat in Monolix with a **non-saturating** structural model (§3.1), using
   SAMBA covariate search, and cross-check against the R result.
5. Report population curves per condition (predicted trajectory + CI) over the
   raw data, plus the covariate-effect table.
