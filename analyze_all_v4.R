# =============================================================================
# GENERALIZED SWEEP v4 — Melanie Prague's inference on every
#   joint x humeral_motion x degree_of_freedom cell, with K chosen on all data
#
# What changed from analyze_all_v3.R — step 1 only, steps 2 to 4 keep their method:
#   * K is chosen on EVERY row of the cell, ex vivo and in vivo together — the
#     reduced model still has no departure, but keeps its random effects. v3
#     chose K on the in-vivo rows alone.
#   * the K grid stops at 10 (v3: 20). Past 10 the extra df made no sense and
#     thin cells asked for more random variances than they had shoulders.
#   * 01_choosing_k draws a bar at the first K whose BIC gain falls below 5 %
#     and below 1 % of the total BIC drop over the grid.
#   * the K used downstream is selected_K in figures/generalized_v4/K_selection.csv,
#     which defaults to the 5 % K and can be edited by hand after looking at
#     the figures. Re-running step 1 never overwrites a selected_K.
#   * step 3 reuses step 1's reduced fit as H0 (same rows, same basis) instead of
#     refitting it.
#   * a cell that is not compared draws ONE grey pooled curve, not the same curve
#     twice in both condition colours (v3 OBS5).
#
# Workflow:
#   ./run_v4_sweep.sh step1     # choose K, write K_selection.csv, draw 01_choosing_k
#   # inspect cells/*/01_choosing_k, edit selected_K
#   ./run_v4_sweep.sh           # step 2 -> 3 -> figures; refits only cells whose K changed
#
# Or step by step, from the repository root:
#   Rscript analyze_all_v4_step1.R     # choose K (shard it: V4_NSHARD / V4_SHARD)
#   Rscript analyze_all_v4_step2.R     # full model at selected_K
#   Rscript analyze_all_v4_step3.R     # LRT + Bonferroni + the master CSVs
#   Rscript analyze_all_v4_figures.R   # seconds — as often as you like
# One cell only:
#   V4_ONLY='scapulothoracic|frontal plane elevation|2' Rscript analyze_all_v4_step1.R
#
# Outputs (figures/generalized_v4/):
#   K_selection.csv         the K decision per cell — edit selected_K
#   cells/<cell>/01_choosing_k .. 06_diagnostics
#   planche_<motion>, planche_diff_<motion>, 00_k_map, 00_forest,
#   00_significance_map, 00_sigmap_elevation, 00_sigmap_rotation
#   00_master_summary.csv, 00_bic_by_k.csv, 00_individual_fit.csv,
#   00_coefficients.csv, 00_SUMMARY.md
# =============================================================================

source("analyze_all_v4_step1.R")
source("analyze_all_v4_step2.R")
source("analyze_all_v4_step3.R")
source("analyze_all_v4_figures.R")
