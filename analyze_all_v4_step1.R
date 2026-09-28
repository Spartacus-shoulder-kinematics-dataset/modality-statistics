# =============================================================================
# GENERALIZED SWEEP v4 — STEP 1: choose K on ALL the data, then check the fit
#
#   1  fit the REDUCED model (no condition terms, random part mirroring K) at
#      K = 1..10 on EVERY row of the cell — ex vivo and in vivo together. The
#      departure is not in this model, but both conditions' shoulders inform
#      how many degrees of freedom the trajectory needs. v3 chose K on the
#      in-vivo rows alone and let K run to 20, which left thin cells asking for
#      more random variances than they had shoulders (00_human_summary.md of
#      v3, OBS1 and OBS2). Past 10 the extra df made no sense anyway.
#   1a per-shoulder RMSE, robust flag, marginal vs subject-specific R2 — kept at
#      EVERY converged K, so a K changed by hand needs no refit to be redrawn.
#
# Three candidate K per cell, all written to K_selection.csv:
#   K_bic    the BIC minimum
#   K_5pct   first K whose BIC gain BIC(K-1) - BIC(K) falls below 5 % of the
#            total BIC drop over the grid, BIC(first K) - min BIC
#   K_1pct   same, below 1 %
# selected_K defaults to K_5pct (K_bic when the 5 % threshold is never met).
# Look at cells/*/01_choosing_k, then edit selected_K by hand: re-running
# step 1 never overwrites a selected_K already in the file.
#
# Run (from the repository root):
#   Rscript analyze_all_v4_step1.R                          # every cell
#   V4_NSHARD=4 V4_SHARD=0 Rscript analyze_all_v4_step1.R   # one shard of four
#   V4_ONLY='scapulothoracic|frontal plane elevation|2' Rscript analyze_all_v4_step1.R
#   V4_REFIT=1 Rscript analyze_all_v4_step1.R               # ignore existing caches
#
# Outputs:
#   cache/step1/<cell>.rds   BIC table, candidate K, per-K individual fit
#   K_selection.csv          one row per cell — THE FILE TO EDIT
#   00_bic_by_k.csv          every cell x every K
#   00_individual_fit.csv    every cell x every shoulder, at the selected K
# =============================================================================

source("analyze_all_v4_common.R")
suppressMessages(library(lcmm))

combos <- select_combos()
cat(sprintf("STEP 1 — BIC over K = %d..%d on the reduced model, ALL rows, %d cells\n",
            min(KS), max(KS), nrow(combos)))
D <- load_long()
t_start <- Sys.time()

for (i in seq_len(nrow(combos))) {
  jt <- combos$joint[i]; mo <- combos$motion[i]; dof <- combos$dof[i]
  tag <- sprintf("[%2d/%d] %-18s %-30s dof%d", i, nrow(combos), jt, substr(mo, 1, 30), dof)

  if (!REFIT && !is.null(cell_load(1, jt, mo, dof))) {
    cat(tag, " cached, skipped\n", sep = ""); next
  }

  P <- prepare_cell(D, jt, mo, dof)
  if (is.null(P$st)) {
    cat(sprintf("%s  SKIP (%s)\n", tag, P$row$note))
    cell_save(1, jt, mo, dof, list(row = P$row, bic = NULL)); next
  }

  # Every row chooses K, compared cell or not: the model has no condition term,
  # and it is the same model step 3 uses as H0 — which is why step 3 can reuse
  # this log-likelihood instead of refitting it.
  st <- P$st
  rows <- list(); perK <- list()
  for (K in KS) {
    b <- cell_basis(st, K, P$ov)
    if (is.null(b)) {
      rows[[length(rows)+1]] <- data.frame(K = K, n_fixed = K+1L, n_random = K+1L,
        converged = FALSE, loglik = NA, BIC = NA, npm = NA,
        n_rows = nrow(st), n_subjects = NA, rmse_max = NA, rmse_median = NA,
        rmse_ratio = NA, r2_ss = NA, r2_marg = NA, n_var_at_zero = NA, seconds = 0)
      next
    }
    dat <- b$df
    t0 <- Sys.time()
    m <- try(hlme(fixed = f_red(b$nm), random = f_rand(b$nm), subject = "IDnum",
                  ng = 1, idiag = TRUE, data = dat, verbose = FALSE, nproc = NPROC),
             silent = TRUE)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ok <- fit_ok(m)
    ps <- if (ok) per_shoulder(m, dat) else NULL
    # variance components pinned at the boundary: the model's own signal that it
    # has been given more random dimensions than the data can support
    nvz <- if (ok) sum(abs(m$best[grep("varcov|std", names(m$best), ignore.case = TRUE)]) < 1e-4) else NA
    rows[[length(rows)+1]] <- data.frame(
      K = K, n_fixed = b$Kc + 1L, n_random = b$Kc + 1L, converged = ok,
      loglik = if (ok) m$loglik else NA, BIC = if (ok) m$BIC else NA,
      npm = if (ok) length(m$best) else NA,
      n_rows = nrow(dat), n_subjects = if (ok) m$ns else NA,
      rmse_max    = if (ok) max(ps$rmse) else NA,
      rmse_median = if (ok) median(ps$rmse) else NA,
      rmse_ratio  = if (ok) max(ps$rmse)/median(ps$rmse) else NA,
      r2_ss   = if (ok) r2_conc(dat$Y, m$pred$pred_ss) else NA,
      r2_marg = if (ok) r2_conc(dat$Y, m$pred$pred_m) else NA,
      n_var_at_zero = nvz, seconds = secs)
    if (ok) perK[[as.character(K)]] <- list(
      ps = ps, pred = pred_rows(m, dat), knots = b$knots, Kc = b$Kc,
      r2_marg = r2_conc(dat$Y, m$pred$pred_m), r2_ss = r2_conc(dat$Y, m$pred$pred_ss),
      # the reduced model's coefficients, so 00_coefficients.csv still carries
      # H0 when step 3 reuses this fit rather than refitting it
      coefs = data.frame(model = "reduced", parameter = names(m$best),
                         estimate = as.vector(m$best), se = sqrt(diag(VarCov(m))),
                         row.names = NULL))
    # Per-K progress, flushed immediately, so a hard cell and a hung one can be
    # told apart from the log.
    if (VERBOSE_K)
      cat(sprintf("      %s dof%d  K=%2d  %6.1fs  %s\n", substr(mo, 1, 18), dof, K, secs,
                  if (ok) sprintf("BIC %9.1f", m$BIC) else "NOT CONVERGED"))
    flush.console()
  }
  bic <- do.call(rbind, rows); rownames(bic) <- NULL
  th  <- k_thresholds(bic)
  n_conv <- sum(bic$converged)

  if (is.na(th$K_bic)) {
    cat(sprintf("%s  NO K CONVERGED\n", tag))
    P$row$note <- "no K converged"
    cell_save(1, jt, mo, dof, list(row = P$row, bic = bic, n_converged = 0L)); next
  }

  cell_save(1, jt, mo, dof, list(
    row = P$row, bic = bic, k_basis = "all rows (ex + in vivo)",
    K_bic = th$K_bic, K_5pct = th$K_5pct, K_1pct = th$K_1pct,
    default_K = th$default_K,
    # a BIC minimum on an end of the grid reports the edge rather than a choice
    bic_edge = th$K_bic %in% range(KS), dBIC_next = th$dBIC_next,
    n_converged = n_conv, perK = perK,
    seconds = sum(bic$seconds, na.rm = TRUE)))

  cat(sprintf("%s  %-7s K_bic=%s  K_5%%=%s  K_1%%=%s  -> default %s  %2d/%d conv  %5.0fs\n",
              tag, P$row$mode, fmtK(th$K_bic), fmtK(th$K_5pct), fmtK(th$K_1pct),
              fmtK(th$default_K), n_conv, length(KS),
              sum(bic$seconds, na.rm = TRUE)))
  flush.console()
}

# Written from the cache rather than from this process's own loop, so a sharded
# run produces the complete files as soon as the last shard lands. Merges into
# K_selection.csv without touching a selected_K already there.
write_step1_tables()

cat(sprintf("\nSTEP 1 done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("Next:  inspect cells/*/01_choosing_k, edit selected_K in ", SEL_CSV,
    "\n       then  Rscript analyze_all_v4_step2.R\n", sep = "")
