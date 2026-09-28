# =============================================================================
# GENERALIZED SWEEP v3 — STEP 3: the likelihood-ratio test, and the plate
#
#   H0  Y ~ ns1 + ... + nsK                  (no condition terms)
#   H1  Y ~ (ns1 + ... + nsK) * cond         (step 2)
#   LR = 2(LL1 - LL0) ~ chi2 with df = K + 1  (gamma_0 plus K shape terms)
#
# H0 is refitted HERE, on all shoulders, at the same K. Step 1's model chose K
# but was fitted on the reference condition alone; differencing its
# log-likelihood against a full-data H1 would be meaningless. Both models must
# sit on exactly the same rows, so the row and subject counts are asserted equal
# before the log-likelihoods are differenced — a silent row drop would make the
# statistic meaningless rather than merely wrong.
#
# MULTIPLICITY IS APPLIED HERE, over the whole plate at once, as Melanie asked:
# Bonferroni is the headline, BH kept beside it for comparability with v2. It is
# applied after every compared cell has landed, and written BEFORE the cache the
# figures read — the plates draw their stars from p_bonf, so a cache saved
# earlier would be missing exactly the column they need.
#
# NO WALD TESTS ANYWHERE. The reporting object is the difference curve.
#
# Run:  Rscript analyze_all_v3_step3.R      (V3_ONLY / V3_SHARD / V3_REFIT as before)
# Outputs:
#   cache/step3/<cell>.rds   H0 loglik and the cell's own LRT
#   cache/v3_fit.rds         everything the figures read
#   00_master_summary.csv    one row per cell
#   00_coefficients.csv      estimates and standard errors, no Wald column
# =============================================================================

source("analyze_all_v3_common.R")
suppressMessages(library(lcmm))

combos <- select_combos()
cat(sprintf("STEP 3 — H0 and the likelihood-ratio test, %d cells\n", nrow(combos)))
D <- load_long()
t_start <- Sys.time()

for (i in seq_len(nrow(combos))) {
  jt <- combos$joint[i]; mo <- combos$motion[i]; dof <- combos$dof[i]
  tag <- sprintf("[%2d/%d] %-18s %-30s dof%d", i, nrow(combos), jt, substr(mo, 1, 30), dof)

  s2 <- cell_load(2, jt, mo, dof)
  if (is.null(s2) || !isTRUE(s2$converged)) { cat(tag, "  no usable step-2 cache\n", sep = ""); next }
  if (!isTRUE(s2$compare)) { cat(sprintf("%s  single — no test\n", tag)); next }
  if (!REFIT && !is.null(cell_load(3, jt, mo, dof))) {
    cat(tag, " cached, skipped\n", sep = ""); next }

  P <- prepare_cell(D, jt, mo, dof)
  if (is.null(P$st)) { cat(sprintf("%s  SKIP (%s)\n", tag, P$row$note)); next }
  K <- s2$K
  b <- cell_basis(P$st, K, P$ov)
  if (is.null(b)) { cat(sprintf("%s  SKIP (basis failed at K=%d)\n", tag, K)); next }

  t0 <- Sys.time()
  m0 <- try(hlme(fixed = f_red(b$nm), random = f_rand(b$nm), subject = "IDnum",
                 ng = 1, idiag = TRUE, data = b$df, verbose = FALSE, nproc = NPROC),
            silent = TRUE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!fit_ok(m0)) {
    cat(sprintf("%s  H0 DID NOT CONVERGE at K=%d (%.0fs)\n", tag, K, secs))
    cell_save(3, jt, mo, dof, list(K = K, converged = FALSE)); next
  }

  # Identical rows and identical subjects, or the difference is not a statistic.
  stopifnot(nrow(b$df) == s2$n_rows, m0$ns == s2$n_subjects)

  LR <- 2 * (s2$loglik_full - m0$loglik)
  DF <- K + 1L
  Pv <- pchisq(LR, df = DF, lower.tail = FALSE)

  coefs0 <- data.frame(model = "reduced", parameter = names(m0$best),
                       estimate = as.vector(m0$best), se = sqrt(diag(VarCov(m0))),
                       row.names = NULL)

  cell_save(3, jt, mo, dof, list(
    K = K, converged = TRUE, loglik_reduced = m0$loglik, loglik_full = s2$loglik_full,
    LR = LR, df = DF, p_lrt = Pv, BIC_reduced = m0$BIC, BIC_full = s2$BIC_full,
    n_rows = nrow(b$df), n_subjects = m0$ns, coefs0 = coefs0, seconds = secs))

  cat(sprintf("%s  K=%2d  %5.0fs  LR %8.2f  df %2d  p %9.2e  %s\n", tag, K, secs, LR, DF, Pv,
              if (Pv < 0.05) "reject H0" else "cannot reject"))
  flush.console()
}

# =============================================================================
# ASSEMBLE — every cell on disk, the plate-wide adjustment, then the cache
# =============================================================================
cat("\n---- assembling the plate ----\n")
A <- collect_cells()          # shared with the figures — see analyze_all_v3_common.R
CUR <- A$cells; master <- A$master; cfl <- A$coefs

write.csv(master, file.path(OUT, "00_master_summary.csv"), row.names = FALSE)
if (!is.null(cfl)) write.csv(cfl, file.path(OUT, "00_coefficients.csv"), row.names = FALSE)
saveRDS(list(cells = CUR, master = master, meta = fit_meta()), CACHE)

cmp <- master$mode == "compare" & !is.na(master$p_lrt)
nt <- sum(cmp)
cat(sprintf("\nSTEP 3 done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat(sprintf("%d cells | compare %d | single %d | skip %d\n", nrow(master),
            sum(master$mode == "compare"), sum(master$mode == "single"),
            sum(master$mode == "skip")))
cat(sprintf("%d tested | significant: %d raw, %d Bonferroni, %d BH (alpha 0.05)\n", nt,
            sum(master$p_lrt[cmp] < 0.05), sum(master$p_bonf[cmp] < 0.05),
            sum(master$p_bh[cmp] < 0.05)))
cat(sprintf("cached %d fitted cells -> %s (%.1f MB)\n", length(CUR), CACHE,
            file.size(CACHE)/1024^2))
cat("Next:  Rscript analyze_all_v3_figures.R\n")
