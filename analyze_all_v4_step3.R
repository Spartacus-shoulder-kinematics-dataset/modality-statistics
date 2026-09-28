# =============================================================================
# GENERALIZED SWEEP v4 — STEP 3: the likelihood-ratio test, and the plate
#
#   H0  Y ~ ns1 + ... + nsK                  (no condition terms)
#   H1  Y ~ (ns1 + ... + nsK) * cond         (step 2)
#   LR = 2(LL1 - LL0) ~ chi2 with df = K + 1  (gamma_0 plus K shape terms)
#
# K is the one step 2 fitted, i.e. the selected_K of K_selection.csv.
#
# H0 IS ALREADY FITTED. Step 1 fits the reduced model on every row of the cell
# at every K, which is exactly H0: same rows, same basis, same random part. Its
# log-likelihood is reused here after checking that the row and subject counts
# match H1's — a silent row drop would make the statistic meaningless rather
# than merely wrong. If step 1 did not converge at that K, or V4_H0_REFIT=1 is
# set, H0 is refitted instead (the latter exists to check the two agree).
#
# MULTIPLICITY IS APPLIED HERE, over the whole plate at once, as Melanie asked:
# Bonferroni is the headline, BH kept beside it. It is written BEFORE the cache
# the figures read, since the plates draw their stars from p_bonf.
#
# NO WALD TESTS ANYWHERE. The reporting object is the difference curve.
#
# Run:  Rscript analyze_all_v4_step3.R      (V4_ONLY / V4_SHARD / V4_REFIT as before)
# Outputs:
#   cache/step3/<cell>.rds   H0 loglik and the cell's own LRT
#   cache/v4_fit.rds         everything the figures read
#   00_master_summary.csv    one row per cell
#   00_coefficients.csv      estimates and standard errors, no Wald column
# =============================================================================

source("analyze_all_v4_common.R")
suppressMessages(library(lcmm))

H0_REFIT <- nzchar(Sys.getenv("V4_H0_REFIT"))
combos <- select_combos()
cat(sprintf("STEP 3 — H0 and the likelihood-ratio test, %d cells%s\n", nrow(combos),
            if (H0_REFIT) " (H0 refitted, not reused)" else ""))
# The raw data is only needed when H0 has to be refitted, so it is read lazily.
D <- NULL
get_D <- function() { if (is.null(D)) D <<- load_long(); D }
t_start <- Sys.time()

for (i in seq_len(nrow(combos))) {
  jt <- combos$joint[i]; mo <- combos$motion[i]; dof <- combos$dof[i]
  tag <- sprintf("[%2d/%d] %-18s %-30s dof%d", i, nrow(combos), jt, substr(mo, 1, 30), dof)

  s2 <- cell_load(2, jt, mo, dof)
  if (is.null(s2) || !isTRUE(s2$converged)) { cat(tag, "  no usable step-2 cache\n", sep = ""); next }
  if (!isTRUE(s2$compare)) { cat(sprintf("%s  single — no test\n", tag)); next }
  K <- as.integer(s2$K)
  s3c <- cell_load(3, jt, mo, dof)
  if (!REFIT && !is.null(s3c) && identical(as.integer(s3c$K), K)) {
    cat(tag, " cached, skipped\n", sep = ""); next }

  s1 <- cell_load(1, jt, mo, dof)
  bk <- if (is.null(s1$bic)) NULL else s1$bic[s1$bic$K == K & s1$bic$converged, , drop = FALSE]
  pk <- if (is.null(s1)) NULL else s1$perK[[as.character(K)]]
  reuse <- !H0_REFIT && !is.null(bk) && nrow(bk) == 1 && !is.null(pk) &&
           bk$n_rows == s2$n_rows && identical(as.integer(bk$n_subjects), as.integer(s2$n_subjects))

  t0 <- Sys.time()
  if (reuse) {
    ll0 <- bk$loglik; bic0 <- bk$BIC; coefs0 <- pk$coefs
    n_rows0 <- bk$n_rows; n_sub0 <- bk$n_subjects; src <- "step 1"
  } else {
    P <- prepare_cell(get_D(), jt, mo, dof)
    if (is.null(P$st)) { cat(sprintf("%s  SKIP (%s)\n", tag, P$row$note)); next }
    b <- cell_basis(P$st, K, P$ov)
    if (is.null(b)) { cat(sprintf("%s  SKIP (basis failed at K=%d)\n", tag, K)); next }
    m0 <- try(hlme(fixed = f_red(b$nm), random = f_rand(b$nm), subject = "IDnum",
                   ng = 1, idiag = TRUE, data = b$df, verbose = FALSE, nproc = NPROC),
              silent = TRUE)
    if (!fit_ok(m0)) {
      cat(sprintf("%s  H0 DID NOT CONVERGE at K=%d\n", tag, K))
      cell_save(3, jt, mo, dof, list(K = K, converged = FALSE)); next
    }
    ll0 <- m0$loglik; bic0 <- m0$BIC; n_rows0 <- nrow(b$df); n_sub0 <- m0$ns; src <- "refit"
    coefs0 <- data.frame(model = "reduced", parameter = names(m0$best),
                         estimate = as.vector(m0$best), se = sqrt(diag(VarCov(m0))),
                         row.names = NULL)
  }
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Identical rows and identical subjects, or the difference is not a statistic.
  stopifnot(n_rows0 == s2$n_rows, n_sub0 == s2$n_subjects)

  LR <- 2 * (s2$loglik_full - ll0)
  DF <- K + 1L
  Pv <- pchisq(LR, df = DF, lower.tail = FALSE)

  cell_save(3, jt, mo, dof, list(
    K = K, converged = TRUE, loglik_reduced = ll0, loglik_full = s2$loglik_full,
    LR = LR, df = DF, p_lrt = Pv, BIC_reduced = bic0, BIC_full = s2$BIC_full,
    n_rows = n_rows0, n_subjects = n_sub0, coefs0 = coefs0, h0_source = src,
    seconds = secs))

  cat(sprintf("%s  K=%2d  H0 from %-6s  LR %8.2f  df %2d  p %9.2e  %s\n", tag, K, src, LR, DF, Pv,
              if (Pv < 0.05) "reject H0" else "cannot reject"))
  flush.console()
}

# =============================================================================
# ASSEMBLE — every cell on disk, the plate-wide adjustment, then the cache
# =============================================================================
cat("\n---- assembling the plate ----\n")
write_step1_tables()          # individual fit follows the selected K
A <- collect_cells()          # shared with the figures — see analyze_all_v4_common.R
CUR <- A$cells; master <- A$master; cfl <- A$coefs

write.csv(master, file.path(OUT, "00_master_summary.csv"), row.names = FALSE)
if (!is.null(cfl)) write.csv(cfl, file.path(OUT, "00_coefficients.csv"), row.names = FALSE)
saveRDS(list(cells = CUR, master = master, meta = fit_meta()), CACHE)

cmp <- master$mode == "compare" & !is.na(master$p_lrt)
cat(sprintf("\nSTEP 3 done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat(sprintf("%d cells | compare %d | single %d | skip %d\n", nrow(master),
            sum(master$mode == "compare"), sum(master$mode == "single"),
            sum(master$mode == "skip")))
cat(sprintf("%d tested | significant: %d raw, %d Bonferroni, %d BH (alpha 0.05)\n", sum(cmp),
            sum(master$p_lrt[cmp] < 0.05), sum(master$p_bonf[cmp] < 0.05),
            sum(master$p_bh[cmp] < 0.05)))
pend <- master[isTRUE_vec(master$refit_pending), ]
if (nrow(pend))
  cat(sprintf("!! %d cell(s) have a selected_K different from their fitted K — re-run step 2:\n%s\n",
              nrow(pend), paste(sprintf("   %s | %s | dof%d", pend$joint, pend$motion, pend$dof),
                                collapse = "\n")))
cat(sprintf("cached %d fitted cells -> %s (%.1f MB)\n", length(CUR), CACHE,
            file.size(CACHE)/1024^2))
cat("Next:  Rscript analyze_all_v4_figures.R\n")
