# =============================================================================
# GENERALIZED SWEEP v2 — iteration 04's model applied to every
#   joint x humeral_motion x degree_of_freedom
#
# What changed from analyze_all.R (v1, the mgcv GAMM):
#   * lcmm::hlme instead of mgcv::bam
#   * ns(TIME, df = 4) natural spline, 10 fixed effects instead of 64 coefficients
#   * random = ~ns(...)  -> a random CURVE per shoulder, not just an intercept
#   * idiag = TRUE       -> D DIAGONAL (no off-diagonal covariances: 21 covariance
#                           parameters from ~44 shoulders was never comfortable)
#   * data TRIMMED to the x-overlap of the two conditions, per cell
#   * WALD test per coefficient + a JOINT Wald on the whole difference block
#
# Per-cell recipe (validated in figures/04_natural_spline_hlme/):
#   thin to <=100 pts/shoulder -> trim to the overlap -> ns basis with knots from
#   the cell's own quantiles and boundary knots from its UNTRIMMED range
#   (tying the boundary knots to the overlap makes the residual tails worse).
#
# Outputs (figures/generalized_v2/):
#   planche_<motion>.png        population curves, rows = joints, cols = DoF
#   planche_diff_<motion>.png   the in-vivo - ex-vivo difference, same layout
#   00_master_summary.csv       one row per cell, incl. Wald counts and joint test
#   00_wald_tests.csv           one row per COEFFICIENT, every cell
#   00_SUMMARY.md               human-readable index
#
# THE SWEEP AND THE FIGURES ARE SEPARATE SCRIPTS. Fitting 63 hlme models takes
# ~12 minutes; drawing takes seconds. The fit writes everything the figures need
# to figures/generalized_v2/cache/v2_fit.rds, so a change to a colour, a label or
# a layout never costs a refit:
#
#   Rscript analyze_all_v2_fit.R       # ~12 min — only when the MODEL changes
#   Rscript analyze_all_v2_figures.R   # seconds — as often as you like
#
# Shared constants and helpers live in analyze_all_v2_common.R; edit them there,
# not in one half.
#
# Run (this file does both, from the repository root):
#       python3 prepare_monolix_data.py
#       Rscript analyze_all_v2.R
# =============================================================================

source("analyze_all_v2_fit.R")
source("analyze_all_v2_figures.R")
