# modality-statistics

**Do in-vivo and ex-vivo shoulder kinematics agree?**

Statistical analysis of the pooled [Spartacus](https://github.com/Spartacus-shoulder-kinematics-dataset)
shoulder-kinematics dataset. The data pool 20 source articles; each observation is a joint
angle (`Y`) measured at a given humerothoracic elevation angle (`TIME`), for one shoulder,
under one of two conditions — **in vivo** or **ex vivo**.

The analysis fits a penalised-spline mixed model (GAMM) per
`joint × humeral_motion × degree_of_freedom` cell, and estimates the in-vivo − ex-vivo
difference curve with a confidence band. Two levels of depth:

- a **worked example** on one cell (scapulothoracic / frontal plane elevation / DoF 2),
  kept as visible iteration sets so the modelling reasoning stays legible;
- a **generalized sweep** over all 72 cells, assembled into one plate per movement.

---

## Read the results correctly

The condition contrast is **exploratory, not causal**. Condition is perfectly confounded
with source study — no study measured both in vivo and ex vivo — so any difference mixes
muscle activation with protocol, measurement and population effects. The effective
replication is the number of **studies**, not rows or shoulders, so p-values are
pseudo-replicated. Read the **effect size in degrees** and the study counts, not the stars.

The full critique and the response to it:
[`docs/STATISTICAL_METHODS_REVIEW.md`](docs/STATISTICAL_METHODS_REVIEW.md) ·
[`docs/STATISTICAL_METHODS_RESPONSE.md`](docs/STATISTICAL_METHODS_RESPONSE.md).

---

## Requirements

There is no dependency manifest in this repo (no `renv.lock`, `DESCRIPTION`, `.Rprofile`,
`pyproject.toml`). The environment this was actually run in:

| | version / source |
| --- | --- |
| R | **4.3.3** (`r-base-core 4.3.3-2build2`) |
| `mgcv` 1.9.1, `nlme` 3.1.164 | apt (`r-recommended` / `r-cran-*`) |
| `dplyr` 1.1.4 | `install.packages()`, user library |
| Python | 3.13 (any Python 3 works) |

`prepare_monolix_data.py` uses **only the standard library** (`csv`, `sys`) — no pip
install, no virtualenv, no `requirements.txt` needed.

The three live R scripts need exactly **`mgcv`, `nlme`, `dplyr`**. Other packages
(`tidyr`, `lme4`, `here`, `lattice`, `MASS`) are used only by the retired/tutorial
scripts; `ggplot2` and `lavaan` are needed by some of those and are **not installed** on
the reference machine — those scripts will not run without installing them first.

## Installation

```bash
# R 4.3.3 plus mgcv, nlme, lattice, MASS
sudo apt install r-base r-recommended

# the rest of the pipeline
Rscript -e 'install.packages("dplyr", repos = "https://cloud.r-project.org")'
```

That is everything needed to reproduce the figures. Optional extras, only if you want to
run the retired / tutorial scripts:

```bash
Rscript -e 'install.packages(c("tidyr","lme4","here","ggplot2","lavaan"), repos = "https://cloud.r-project.org")'
```

Check your setup:

```bash
Rscript -e 'sapply(c("mgcv","nlme","dplyr"), requireNamespace, quietly = TRUE)'
# all three must be TRUE
```

**Monolix is not required.** No Monolix project file lives in this repo — the R path is
fully worked out. Monolix appears only as an export CSV and a GUI recipe in
[`docs/METHODOLOGY.md`](docs/METHODOLOGY.md) §2, for whoever gets a license.

---

## Data

**No data file is versioned here** — the raw dataset is ~300 MB and the derived files are
reproducible from it. Download the raw dataset first, then generate the rest locally.

### Step 1 — get the raw dataset

The easiest way: grab the published release asset from the Spartacus dataset repository
and save it under the name the scripts expect.

```bash
curl -L -o corrected_confident_data.csv \
  https://github.com/Spartacus-shoulder-kinematics-dataset/shoulder-kinematics/releases/download/0.3.1/spartacus.csv
```

Alternatively, regenerate it from source by running the data-building code in
[Spartacus-shoulder-kinematics-dataset/shoulder-kinematics](https://github.com/Spartacus-shoulder-kinematics-dataset/shoulder-kinematics),
which produces `spartacus/dataset/corrected_confident_data.csv`; copy that file into this
repository's root. Use this route if you need a version other than the released one, or if
you want to trace how the pooled dataset is assembled.

### Step 2 — build the derived files

```bash
python3 prepare_monolix_data.py
```

### The three files

| file | role |
| --- | --- |
| `corrected_confident_data.csv` | **raw input**, ~306 MB, 1,381,818 rows × 26 columns — downloaded in step 1, never committed. Only rows with `unit == "rad"` (angular data) are analysed; `unit == "mm"` rows are translations and are dropped. |
| `spartacus_angles_long.csv` | derived, ~97 MB — all angular data in tidy long form. Input to `analyze_all.R`. |
| `monolix_st_frontal_dof2.csv` | derived, ~1.7 MB — the single worked-example subset. Input to `analysis_shoulder.R`. |

All three are gitignored. Column roles of the derived files:

| column | meaning |
| --- | --- |
| `ID` | unit = one shoulder, built as `article_shoulderid` |
| `TIME` | the regressor **x** — humerothoracic elevation angle (not clock time; the name follows the Monolix convention) |
| `Y` | the observation — the joint angle |
| `in_vivo` | condition, `True` / `False` |
| `study` | source article, kept separate from `ID` so study-level effects can be fitted |

`spartacus_angles_long.csv` additionally carries `joint`, `humeral_motion` and
`degree_of_freedom`. Joints: sternoclavicular, acromioclavicular, scapulothoracic,
glenohumeral. Movements: frontal / sagittal / scapular plane elevation, horizontal
flexion, internal-external rotation at 0° and 90° abduction.

---

## Running the analysis

All scripts use relative paths — **run them from the repository root**, after downloading
the raw dataset as described under [Data](#data).

```bash
python3 prepare_monolix_data.py     # raw CSV -> the two derived long CSVs
Rscript analysis_shoulder.R         # worked example -> figures/00_, 01_, 02_
Rscript analyze_all.R               # 72-cell sweep -> figures/generalized/  (long-running)
Rscript review_response_analysis.R  # reviewer-response evidence, console output only
```

`analyze_all.R` fits a GAMM per cell over the 97 MB long file; expect it to run for a
while and to rewrite every PNG under `figures/generalized/`.

The scripts create their own `figures/` subdirectories, but they do **not** regenerate the
hand-written [`figures/README.md`](figures/README.md) — keep that file if you ever
`rm -rf figures/`.

---

## The model

The retained specification, validated in `analysis_shoulder.R` and reused for every cell:

```r
bam(Y ~ condO + s(TIME) + s(TIME, by = condO) + s(ID, bs = "re"),
    method = "fREML", discrete = TRUE, rho = rho, AR.start = ar_start)
```

- `s(TIME)` — the common trajectory; `s(TIME, by = condO)` — the **shape** difference
  between conditions; `condO` — a constant **level** shift; `s(ID, bs = "re")` — a random
  intercept per shoulder.
- Residuals within a curve are massively autocorrelated (lag-1 ρ ≈ 0.996 on the worked
  example), so an AR(1) correction is applied **within shoulder** via `rho` / `AR.start`.
  Without it the inference is wildly over-optimistic — BIC −31,454 with AR(1) versus
  111,037 without.
- A parametric 4-parameter logistic NLME was tried first and abandoned: the data do not
  saturate over the observed range, so the plateau parameters are non-identifiable
  (BIC 125,717). The abandoned fit is kept in `figures/01_sigmoid_nlme/` on purpose.
- The difference curve is evaluated **only on the x-overlap** of the two conditions, and
  the 33 compared cells get Benjamini–Hochberg FDR adjustment.

The full statistical rationale — model family, covariate selection, diagnostics, and why
latent growth curve modelling was rejected — is in
[`docs/METHODOLOGY.md`](docs/METHODOLOGY.md).

---

## Repository map

### Live pipeline

| script | what it does | reads | writes |
| --- | --- | --- | --- |
| [`prepare_monolix_data.py`](prepare_monolix_data.py) | filters to angular data and reshapes to long format, in one streaming pass | `corrected_confident_data.csv` | `spartacus_angles_long.csv`, `monolix_st_frontal_dof2.csv` |
| [`analysis_shoulder.R`](analysis_shoulder.R) | the pedagogical worked example: exploration → sigmoid NLME (rejected) → spline mixed model + AR(1) → difference curve | `monolix_st_frontal_dof2.csv` | `figures/00_data_exploration/`, `figures/01_sigmoid_nlme/`, `figures/02_spline_mixed_model/` |
| [`analyze_all.R`](analyze_all.R) | the same model over all 72 joint × movement × DoF cells, plus one plate per movement and FDR adjustment | `spartacus_angles_long.csv` | `figures/generalized/` (plates, per-cell drill-downs, `00_master_summary.csv`, `00_SUMMARY.md`) |
| [`review_response_analysis.R`](review_response_analysis.R) | evidence for the reviewer response: study random effect, corrected ρ, signed effect size, leave-one-study-out | `spartacus_angles_long.csv` | console only |

### Retired / scratch — not part of the pipeline

| script | status |
| --- | --- |
| [`PolynomialFit.R`](PolynomialFit.R) | descriptive polynomial fits on group-mean curves. **Not valid for inference** — see `docs/METHODOLOGY.md` §4. Needs `ggplot2`. |
| [`first_comparison.R`](first_comparison.R), [`LGM_spartacus.R`](LGM_spartacus.R) | the abandoned latent-growth-curve route. Both carry a bare `filter(...)` syntax bug and are kept for the record only. |
| [`LGM_test.R`](LGM_test.R), [`latent_growth_curve_tuto.R`](latent_growth_curve_tuto.R) | teaching scripts on **simulated** data; read and write nothing. |
| [`load_a_csv.R`](load_a_csv.R) | scratch loader with a stale hard-coded path. |

`article/` holds the LaTeX manuscript sources and is not documented here.

---

## Documentation

| document | contents |
| --- | --- |
| [`figures/README.md`](figures/README.md) | the figure-by-figure narrative — what each iteration set shows and why the model moved on |
| [`figures/generalized/00_SUMMARY.md`](figures/generalized/00_SUMMARY.md) | generated index of the sweep: 33 compared / 30 single-condition / 9 skipped cells, sorted by effect size |

---

## Known repository issues

- The values labelled `rad` in the raw data have degree-like magnitudes; the provenance of
  that unit label is still unresolved (noted in the review and in the results digest).

## License

MIT — see [`LICENSE`](LICENSE).
