# =============================================================================
# GENERALIZED SWEEP v3 — Melanie Prague's inference, applied to every
#   joint x humeral_motion x degree_of_freedom cell
#
# What changed from analyze_all_v2.R (iteration 04's model on every cell):
#   * ns df is CHOSEN PER CELL by BIC over K = 1..20, not fixed at 4
#   * the df is chosen on the REDUCED model (no condition terms), fitted on the
#     REFERENCE condition's own rows — in vivo is the reference now, so the
#     shape of the reference trajectory is settled from the reference's own data
#   * IN VIVO is the reference level, so the departure term is ex vivo and the
#     difference curve is EX VIVO - IN VIVO — mirrored against v2
#   * the global test is a LIKELIHOOD-RATIO TEST, df = K + 1, replacing v2's
#     joint Wald chi2 on the gamma block
#   * multiplicity is BONFERRONI over the whole plate, not BH
#   * NO WALD TESTS ANYWHERE, and no knots drawn on any figure
#   * one folder of figures PER CELL, plus a plate showing how every K was chosen
#
# Melanie's plan, verbatim (docs/meeting-questions.md):
#   1) Best BIC donc le degre de spline change en fonction du marqueur, dans un
#      modele reduit sans in vivo departure.  (a) verification model predict vs
#      data, BLUP.  (b) choix des knots au quantile.
#   2) Construire le modele de spline avec le MEME degre que 1) en ajoutant le
#      in vivo departure.
#   3) Test du rapport de vraisemblance, pour chaque marqueur — ajustement
#      multiple a faire ici, sur la planche au complet. Prendre bonferroni.
#   4) PAS de test de wald, on abandonne !
#
# THE SWEEP IS SPLIT BY STEP, AND FITTING IS SEPARATE FROM DRAWING. Step 1 is
# the long pole — 20 fits per cell against steps 2-3's two — so it is the step
# worth sharding. Each step caches one RDS per cell as that cell finishes, so a
# killed run resumes rather than restarts, and each step reads only the step
# before it: changing the LRT or the multiplicity adjustment costs seconds.
#
#   Rscript analyze_all_v3_step1.R     # choose K   (hours; shard it)
#   Rscript analyze_all_v3_step2.R     # full model (minutes)
#   Rscript analyze_all_v3_step3.R     # LRT + Bonferroni + the master CSVs
#   Rscript analyze_all_v3_figures.R   # seconds — as often as you like
#
# Sharding across processes, since hlme(nproc=) already forks internally:
#   for s in 0 1 2 3; do V3_NSHARD=4 V3_SHARD=$s Rscript analyze_all_v3_step1.R & done; wait
# One cell only:
#   V3_ONLY='scapulothoracic|frontal plane elevation|2' Rscript analyze_all_v3_step1.R
#
# Shared constants, the per-cell recipe and the renderer live in
# analyze_all_v3_common.R; edit them there, not in one of the steps.
#
# Run (this file does all four, from the repository root):
#       python3 prepare_monolix_data.py
#       Rscript analyze_all_v3.R
#
# Outputs (figures/generalized_v3/):
#   cells/<cell>/01_choosing_k .. 06_diagnostics   one folder per cell
#   planche_<motion>, planche_diff_<motion>        curves and differences
#   00_k_map                the df selection for the whole plate
#   00_forest, 00_significance_map, 00_sigmap_elevation, 00_sigmap_rotation
#   00_master_summary.csv   one row per cell, incl. K, LR, p, Bonferroni
#   00_bic_by_k.csv         every cell x every K — the selection evidence
#   00_individual_fit.csv   every cell x every shoulder
#   00_coefficients.csv     estimates and standard errors, NO Wald column
#   00_SUMMARY.md           human-readable index
# =============================================================================

source("analyze_all_v3_step1.R")
source("analyze_all_v3_step2.R")
source("analyze_all_v3_step3.R")
source("analyze_all_v3_figures.R")
