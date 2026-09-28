# =============================================================================
# GENERALIZED SWEEP v3 — STEP 1: choose K, then check it captures individuals
#
#   1  fit the REDUCED model (no condition terms) at K = 1..20 and take the BIC
#      minimum, on the REFERENCE CONDITION'S OWN ROWS. in vivo is the reference,
#      so the shape of the reference trajectory is settled from the reference's
#      own data — the ex-vivo shoulders get no vote on how many degrees of
#      freedom the common curve needs, since they are what will be tested
#      against it.
#   1a at the retained K, check no shoulder is being left behind: per-shoulder
#      RMSE, a robust flag, and marginal vs subject-specific concordance R2.
#
# The random part MIRRORS K, as in iterations 04 and 06: K = 10 means 11 random
# effects per shoulder. Selection runs over the CONVERGED rows only.
#
# This is the long pole of the sweep — 20 fits per cell against step 2/3's two —
# so it is the step worth sharding across processes.
#
# Run (from the repository root):
#   Rscript analyze_all_v3_step1.R                       # every cell
#   V3_NSHARD=4 V3_SHARD=0 Rscript analyze_all_v3_step1.R   # one shard of four
#   V3_ONLY='scapulothoracic|frontal plane elevation|2' Rscript analyze_all_v3_step1.R
#   V3_REFIT=1 Rscript analyze_all_v3_step1.R            # ignore existing caches
#
# Outputs:
#   cache/step1/<cell>.rds   K table, K_star, individual fit, reference predictions
#   00_bic_by_k.csv          every cell x every K — the selection evidence
#   00_individual_fit.csv    every cell x every shoulder
# =============================================================================

source("analyze_all_v3_common.R")
suppressMessages(library(lcmm))

combos <- select_combos()
cat(sprintf("STEP 1 — BIC over K = %d..%d on the reduced model, %d cells\n",
            min(KS), max(KS), nrow(combos)))
cat(sprintf("reference condition: %s (departure term = %s)\n", LEVELS[1], LEVELS[2]))
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
    cell_save(1, jt, mo, dof, list(row = P$row, bic = NULL, K_star = NA_integer_)); next
  }

  # Which rows choose K, recorded in k_basis so the figure can say so.
  #
  # Compared cell: the reference condition alone. in vivo is the reference, so
  # the shape of the reference trajectory is settled from the reference's own
  # data — the ex-vivo shoulders get no vote on how many degrees of freedom the
  # common curve needs, since they are what will be tested against it.
  #
  # NOT compared: every row. The reduced model has no condition term and step 2
  # fits it on the whole cell, so K has to be chosen on the whole cell too or
  # the df would be selected from data the model is not fitted to. Restricting
  # to one condition here is not a smaller version of the same idea — for
  # sternoclavicular / scapular plane elevation / DoF1 (1 ex vivo, 25 in vivo,
  # "single" only because 1 < MIN_PER_COND) it would pick the lone ex-vivo
  # shoulder and choose the df from one curve.
  st <- P$st
  if (P$compare) {
    k_basis <- LEVELS[1]
    dref_i  <- st$cond == k_basis
  } else {
    k_basis <- "all rows (not compared)"
    dref_i  <- rep(TRUE, nrow(st))
  }
  if (!any(dref_i)) { cat(sprintf("%s  SKIP (no %s rows)\n", tag, k_basis))
    P$row$note <- paste("no", k_basis, "rows")
    cell_save(1, jt, mo, dof, list(row = P$row, bic = NULL, K_star = NA_integer_)); next }

  rows <- list(); keep <- list()
  for (K in KS) {
    b <- cell_basis(st, K, P$ov)
    if (is.null(b)) {
      rows[[length(rows)+1]] <- data.frame(K = K, n_fixed = K+1L, n_random = K+1L,
        converged = FALSE, loglik = NA, BIC = NA, npm = NA, rmse_max = NA,
        rmse_median = NA, rmse_ratio = NA, r2_ss = NA, r2_marg = NA,
        n_var_at_zero = NA, seconds = 0)
      next
    }
    dref <- b$df[dref_i, ]
    dref$IDnum <- as.integer(droplevels(dref$ID))
    t0 <- Sys.time()
    m <- try(hlme(fixed = f_red(b$nm), random = f_rand(b$nm), subject = "IDnum",
                  ng = 1, idiag = TRUE, data = dref, verbose = FALSE, nproc = NPROC),
             silent = TRUE)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ok <- fit_ok(m)
    ps <- if (ok) per_shoulder(m, dref) else NULL
    # variance components pinned at the boundary: the model's own signal that it
    # has been given more random dimensions than the data can support
    nvz <- if (ok) sum(abs(m$best[grep("varcov|std", names(m$best), ignore.case = TRUE)]) < 1e-4) else NA
    rows[[length(rows)+1]] <- data.frame(
      K = K, n_fixed = b$Kc + 1L, n_random = b$Kc + 1L, converged = ok,
      loglik = if (ok) m$loglik else NA, BIC = if (ok) m$BIC else NA,
      npm = if (ok) length(m$best) else NA,
      rmse_max    = if (ok) max(ps$rmse) else NA,
      rmse_median = if (ok) median(ps$rmse) else NA,
      rmse_ratio  = if (ok) max(ps$rmse)/median(ps$rmse) else NA,
      r2_ss   = if (ok) r2_conc(dref$Y, m$pred$pred_ss) else NA,
      r2_marg = if (ok) r2_conc(dref$Y, m$pred$pred_m) else NA,
      n_var_at_zero = nvz, seconds = secs)
    if (ok) keep[[as.character(K)]] <- list(ps = ps, pred = pred_rows(m, dref),
                                            Y = dref$Y, knots = b$knots, Kc = b$Kc)
    # Per-K progress, flushed immediately. A cell is 20 fits and the slow ones
    # take tens of minutes, so without this a shard is silent for an hour and
    # there is no way to tell a hard cell from a hung one.
    if (VERBOSE_K)
      cat(sprintf("      %s dof%d  K=%2d  %6.1fs  %s\n", substr(mo, 1, 18), dof, K, secs,
                  if (ok) sprintf("BIC %9.1f", m$BIC) else "NOT CONVERGED"))
    flush.console()
  }
  bic <- do.call(rbind, rows); rownames(bic) <- NULL

  okb <- bic[bic$converged, ]
  if (!nrow(okb)) {
    cat(sprintf("%s  NO K CONVERGED\n", tag))
    P$row$note <- "no K converged"
    cell_save(1, jt, mo, dof, list(row = P$row, bic = bic, K_star = NA_integer_)); next
  }
  K_star <- okb$K[which.min(okb$BIC)]
  # A minimum sitting on an edge of the grid is BIC reporting the edge, not
  # making a choice. Recorded so the figure can say so rather than pass it off
  # as a selection.
  bic_edge  <- K_star == min(KS) || K_star == max(KS)
  dBIC_next <- if (nrow(okb) > 1) sort(okb$BIC)[2] - min(okb$BIC) else NA_real_

  kk  <- keep[[as.character(K_star)]]
  ind <- flag_shoulders(kk$ps)
  R2_MARG <- r2_conc(kk$Y, kk$pred$pred_m)
  R2_SS   <- r2_conc(kk$Y, kk$pred$pred_ss)

  cell_save(1, jt, mo, dof, list(
    row = P$row, bic = bic, K_star = K_star, k_basis = k_basis,
    bic_edge = bic_edge, dBIC_next = dBIC_next, n_converged = nrow(okb),
    knots = kk$knots, Kc = kk$Kc, ind = ind, pred_ref = kk$pred,
    r2_marg = R2_MARG, r2_ss = R2_SS,
    seconds = sum(bic$seconds, na.rm = TRUE)))

  cat(sprintf("%s  %-7s K*=%2d%s  BIC %8.1f (+%.0f)  %2d/%d conv  rmse max %.2f med %.2f  flagged %d/%d  %5.0fs\n",
              tag, P$row$mode, K_star, if (bic_edge) " EDGE" else "     ",
              min(okb$BIC), dBIC_next, nrow(okb), length(KS),
              max(ind$rmse), median(ind$rmse), sum(ind$flagged), nrow(ind),
              sum(bic$seconds, na.rm = TRUE)))
  flush.console()
}

# ---- gather every cell on disk into the two long tables ---------------------
# Written from the cache rather than from this process's own loop, so a sharded
# run produces the complete tables as soon as the last shard lands.
cb <- all_combos(); bl <- list(); il <- list()
for (i in seq_len(nrow(cb))) {
  o <- cell_load(1, cb$joint[i], cb$motion[i], cb$dof[i]); if (is.null(o)) next
  idc <- data.frame(joint = cb$joint[i], motion = cb$motion[i], dof = cb$dof[i])
  if (!is.null(o$bic)) bl[[length(bl)+1]] <- cbind(idc, K_star = o$K_star, o$bic)
  if (!is.null(o$ind)) il[[length(il)+1]] <- cbind(idc, K = o$K_star, o$ind)
}
if (length(bl)) write.csv(do.call(rbind, bl), file.path(OUT, "00_bic_by_k.csv"), row.names = FALSE)
if (length(il)) write.csv(do.call(rbind, il), file.path(OUT, "00_individual_fit.csv"), row.names = FALSE)

cat(sprintf("\nSTEP 1 done in %.1f min | %d cells cached | %d with a K\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins")),
            length(bl), sum(sapply(bl, function(z) !is.na(z$K_star[1])))))
cat("Next:  Rscript analyze_all_v3_step2.R\n")
