# =============================================================================
# GENERALIZED SWEEP v3 — THE FIGURE HALF
#
# Draws everything from figures/generalized_v3/cache/v3_fit.rds. FITS NOTHING,
# reads no raw data, loads no lcmm, takes seconds. Edit a colour or a label and
# re-run this; the sweep in analyze_all_v3_step{1,2,3}.R only has to happen when
# the MODEL changes.
#
# Run (from the repository root):  Rscript analyze_all_v3_figures.R
#
# Outputs (figures/generalized_v3/):
#   cells/<cell>/01_choosing_k       .. 06_diagnostics   one folder per cell
#   planche_<motion>, planche_diff_<motion>              population curves / difference
#   00_k_map                          the df selection for the whole plate
#   00_forest                         every compared cell on one difference axis
#   00_significance_map, 00_sigmap_elevation, 00_sigmap_rotation
#   00_SUMMARY.md
#
# NO WALD ANYWHERE, and NO KNOTS DRAWN ANYWHERE.
# =============================================================================

source("analyze_all_v3_common.R")

# ---- drawing in parallel ----------------------------------------------------
# Every per-cell folder and every planche is independent: a task opens its own
# devices and writes its own files, so forking is safe and there is nothing to
# merge afterwards. ~63 cells x 6 figures x 2 devices is the bulk of the work.
#
# mclapply swallows an error into a try-error object instead of throwing, so a
# broken panel would otherwise show up as a silently missing file. par_do()
# collects them and stops with the cell named.
FIGCORES <- max(1L, min(8L, parallel::detectCores() - 1L))
if (nzchar(Sys.getenv("V3_FIGCORES")))
  FIGCORES <- max(1L, as.integer(Sys.getenv("V3_FIGCORES")))

par_do <- function(xs, f, what) {
  res <- if (FIGCORES > 1L && .Platform$OS.type == "unix")
    parallel::mclapply(xs, f, mc.cores = min(FIGCORES, length(xs)), mc.preschedule = FALSE)
  else lapply(xs, f)
  bad <- vapply(res, function(z) inherits(z, "try-error"), logical(1))
  if (any(bad)) {
    for (i in which(bad))
      message("  FAILED ", what, " [", xs[[i]], "]: ", conditionMessage(attr(res[[i]], "condition")))
    stop(sum(bad), " of ", length(xs), " ", what, " failed — see the messages above.")
  }
  invisible(res)
}

# Normally the plate comes from the cache step 3 assembled. When that does not
# exist yet — or V3_PARTIAL=1 is set — read the per-cell RDS files directly, so a
# sweep still in progress can be drawn from whatever has landed. A cell with only
# step 1 has its BIC-vs-K table and its individual fit, which is enough for
# 01_choosing_k and 05_individual_fit; the four figures that need fitted curves
# are skipped until step 2 reaches it.
# The cache can also be STALE rather than absent: step 3 assembles it from the
# cells that existed when it ran, so a sweep that carried on afterwards leaves
# per-cell files newer than the snapshot. Reading the snapshot then silently
# drops those cells. Compare mtimes and prefer the per-cell files when they are
# ahead — the same silent-staleness trap the per-cell figure wipe closes.
newer_cells <- function() {
  if (!file.exists(CACHE)) return(TRUE)
  f <- list.files(CACHE_DIR, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  f <- f[basename(f) != basename(CACHE)]
  if (!length(f)) return(FALSE)
  max(file.mtime(f)) > file.mtime(CACHE)
}
STALE   <- newer_cells()
PARTIAL <- nzchar(Sys.getenv("V3_PARTIAL")) || !file.exists(CACHE) || STALE
if (STALE && file.exists(CACHE))
  cat("cache is behind the per-cell files — reading those instead\n")
if (PARTIAL) {
  A <- collect_cells()
  CUR <- A$cells; master <- A$master
  if (!length(CUR))
    stop("nothing to draw: no per-cell caches under ", CACHE_DIR,
         "\n  Run:  Rscript analyze_all_v3_step1.R")
  n2 <- sum(master$has_step2, na.rm = TRUE); n3 <- sum(master$has_step3, na.rm = TRUE)
  cat(sprintf("PARTIAL VIEW — read %d cells from %s/step*/\n", length(CUR), CACHE_DIR))
  cat(sprintf("  %d have step 1 (K chosen) · %d have step 2 (curves) · %d have step 3 (LRT)\n",
              length(CUR), n2, n3))
  if (n3 < n2 || n2 < length(CUR))
    cat("  figures needing a stage a cell has not reached are skipped, not faked.\n")
} else {
  cache <- readRDS(CACHE)
  check_meta(cache$meta)
  CUR <- cache$cells; master <- cache$master
  cat(sprintf("loaded %d fitted cells from %s (built %s)\n", length(CUR), CACHE,
              format(cache$meta$fitted_at, "%Y-%m-%d %H:%M")))
}

# range() over a plate row that may be entirely empty. 9 of the 72 cells are
# never modelled and whole motions have no compared cell at all, so a planche
# row with nothing in it is normal — range(NULL) would warn on every one.
rng_or <- function(v, default) {
  v <- unlist(v); v <- v[is.finite(v)]
  if (!length(v)) default else range(v)
}

mrow <- function(jt, mo, dof) {
  m <- master[master$joint == jt & master$motion == mo & master$dof == dof, ]
  if (nrow(m) == 1) m else NULL
}
cell_title <- function(jt, mo, dof)
  sprintf("%s · DoF%d — %s\n%s", capj(jt), dof, DOF_LEGEND[[jt]][dof], mo)
# p_bonf is the headline everywhere: Melanie asked for Bonferroni over the whole
# plate, and every star on every plate is drawn from it.
pstar <- function(ri) if (is.null(ri)) "" else stars(ri$p_bonf)

# =============================================================================
# PER-CELL FIGURES — one folder per cell, one file per figure type
# =============================================================================
draw_choosing_k <- function(e, ri, ttl) {
  b <- e$bic; ok <- b[b$converged, ]
  op <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 3.0, 1.0), mgp = c(2.5, 0.7, 0),
            oma = c(0, 0, 3.4, 0)); on.exit(par(op), add = TRUE)

  # Log y: at K = 1 the reduced model is a straight line and its BIC dwarfs the
  # rest of the grid, which on a linear axis squashes the minimum — the one part
  # of the curve the figure exists to show — into the bottom few pixels. The
  # ticks stay real BIC values. Linear fallback if any BIC is non-positive.
  logy <- if (all(ok$BIC > 0)) "y" else ""
  plot(b$K, b$BIC, type = "n", log = logy, xlab = "K (ns columns)", ylab = "BIC",
       main = "BIC over the grid", cex.main = 0.95)
  lines(ok$K, ok$BIC, col = BLUE, lwd = 2)
  points(ok$K, ok$BIC, pch = 16, col = BLUE, cex = 0.8)
  bad <- b[!b$converged, ]
  if (nrow(bad)) {
    # A K that did not converge is not a K with a bad BIC — it has none. Marked
    # on the axis so the gap in the curve reads as a failure, not as a dip.
    u <- par("usr")
    points(bad$K, rep(u[3] + 0.02*diff(u[3:4]), nrow(bad)), pch = 4, cex = 0.9, col = "grey45")
    text(u[1] + 0.03*diff(u[1:2]), u[3] + 0.075*diff(u[3:4]),
         sprintf("× = did not converge (%d of %d)", nrow(bad), nrow(b)),
         adj = c(0, 0), cex = 0.62, col = "grey40")
  }
  points(e$K_star, ok$BIC[ok$K == e$K_star], pch = 21, bg = BLUE, col = "white",
         cex = 2.0, lwd = 1.6)
  # A minimum on an edge of the grid is BIC reporting the edge rather than
  # choosing. Said on the face of the figure, since the dot alone looks identical.
  lab <- if (isTRUE(e$bic_edge))
    sprintf("K = %d — AT THE GRID EDGE, not an interior minimum", e$K_star)
  else sprintf("K = %d — BIC minimum, next best %.0f higher", e$K_star, e$dBIC_next)
  mtext(lab, side = 3, line = 0.2, cex = 0.68,
        col = if (isTRUE(e$bic_edge)) "#b2182b" else "grey25",
        font = if (isTRUE(e$bic_edge)) 2 else 1)
  # The random part mirrors K, so this fit estimates K+1 variances from
  # nrow(e$ind) shoulders. At a ratio of 1 or more there are as many variance
  # parameters as curves, and the retained K describes how much flexibility the
  # optimiser would bear rather than what the trajectory needs. Said on the
  # figure, since BIC's subject-count penalty will not have resisted it.
  if (isTRUE(ri$overparam))
    mtext(sprintf("%d random variances from %d shoulders — K is not identifiable here",
                  e$K_star + 1, ri$n_shoulders_ref),
          side = 3, line = -0.75, cex = 0.62, col = "#b2182b", font = 2)

  yl <- range(c(ok$rmse_max, ok$rmse_median), na.rm = TRUE)
  if (!all(is.finite(yl)) || yl[1] <= 0) yl <- c(1e-3, 1)
  plot(ok$K, ok$rmse_max, type = "b", log = "y", ylim = yl, pch = 16, cex = 0.8,
       col = COL[["ex vivo"]], lwd = 2, xlab = "K (ns columns)",
       ylab = "per-shoulder RMSE (°)", main = "What the df buys each shoulder",
       cex.main = 0.95)
  lines(ok$K, ok$rmse_median, type = "b", pch = 16, cex = 0.8,
        col = COL[["in vivo"]], lwd = 2)
  abline(v = e$K_star, lty = 2, col = "grey55")
  legend("topright", bty = "n", cex = 0.68, lwd = 2, pch = 16,
         col = c(COL[["ex vivo"]], COL[["in vivo"]]),
         legend = c("worst shoulder", "median shoulder"))
  basis_txt <- if (identical(e$k_basis, LEVELS[1]))
    sprintf("the %s data alone (the reference condition)", e$k_basis)
  else "every row — this cell is not compared, so there is no reference condition"
  # Two lines, not one: joined they overrun the device and get clipped at BOTH
  # ends, losing the joint name on the left and the convergence count on the right.
  mtext(ttl, outer = TRUE, line = 1.0, cex = 0.74, font = 2)
  mtext(sprintf("K chosen by BIC on %s — %d of %d converged",
                basis_txt, e$n_converged, nrow(e$bic)),
        outer = TRUE, line = 0.1, cex = 0.62, col = "grey30")
}

draw_population <- function(e, ri, ttl) {
  op <- par(mar = c(4.2, 4.4, 4.2, 1.0), mgp = c(2.5, 0.7, 0)); on.exit(par(op), add = TRUE)
  raw <- e$raw
  plot(raw$TIME, raw$Y, col = adjustcolor(COL[as.character(raw$cond)], 0.16),
       pch = 16, cex = 0.35, xlab = "thoracohumeral elevation (°)", ylab = "joint angle (°)",
       main = "")
  for (cc in unique(e$curves$cond)) { g <- e$curves[e$curves$cond == cc, ]
    polygon(c(g$x, rev(g$x)), c(g$lo, rev(g$hi)), col = adjustcolor(COL[cc], 0.28), border = NA)
    lines(g$x, g$fit, col = COL[cc], lwd = 2.6) }
  legend("topleft", bty = "n", cex = 0.72, lwd = 2.6, col = COL[names(COL)],
         legend = names(COL))
  title(main = ttl, cex.main = 0.86, line = 2.4)
  mtext(sprintf("population curves with 95%% bands, K = %d   ·   %d ex / %d in shoulders",
                e$K, ri$n_ex, ri$n_in), side = 3, line = 0.4, cex = 0.68, col = "grey30")
}

draw_difference <- function(e, ri, ttl) {
  op <- par(mar = c(4.2, 4.4, 4.6, 1.0), mgp = c(2.5, 0.7, 0)); on.exit(par(op), add = TRUE)
  z <- e$dif
  yl <- range(0, z$lo, z$hi); yl <- yl + c(-1, 1)*0.06*diff(yl)
  plot(z$x, z$fit, type = "n", ylim = yl, xlab = "thoracohumeral elevation (°)",
       ylab = "ex vivo − in vivo (°)", main = "")
  sg <- sig_segments(e)
  if (!is.null(sg) && nrow(sg))
    for (q in seq_len(nrow(sg)))
      rect(sg$lo[q], yl[1], sg$hi[q], yl[2], col = adjustcolor(BLUE, 0.07), border = NA)
  polygon(c(z$x, rev(z$x)), c(z$lo, rev(z$hi)), col = adjustcolor(BLUE, 0.25), border = NA)
  abline(h = 0, lty = 2, col = "grey40"); lines(z$x, z$fit, lwd = 2.8, col = BLUE)
  title(main = ttl, cex.main = 0.86, line = 2.8)
  mtext(sprintf("LR = %.2f, df = %d, p = %.2e   ·   Bonferroni p = %.2e %s",
                ri$LR, ri$lrt_df, ri$p_lrt, ri$p_bonf, pstar(ri)),
        side = 3, line = 0.9, cex = 0.70, col = "grey20", font = 2)
  mtext(sprintf("mean |Δ| = %.2f°, max %.2f°  ·  band excludes 0 over %.0f%% of the fitted range (shaded)",
                ri$diff_mean_abs, ri$diff_max_abs, 100*ri$sig_frac_x),
        side = 3, line = 0.1, cex = 0.64, col = "grey35")
}

draw_obs_pred <- function(e, ri, ttl) {
  p <- e$pred_full
  op <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 3.0, 1.0), mgp = c(2.5, 0.7, 0),
            oma = c(0, 0, 2.6, 0)); on.exit(par(op), add = TRUE)
  rg <- range(c(p$obs, p$pred_m, p$pred_ss))
  # 1:1 line, not a regression line: the question is agreement, not correlation.
  # And concordance R2 rather than a squared correlation, so a shoulder fitted
  # with a constant offset cannot score 0.99 while sitting off the diagonal.
  for (w in c("pred_m", "pred_ss")) {
    plot(p[[w]], p$obs, pch = 16, cex = 0.3, xlim = rg, ylim = rg,
         col = adjustcolor(COL[p$cond], 0.20), xlab = "predicted (°)", ylab = "observed (°)",
         main = if (w == "pred_m") "marginal — population curve only"
                else "subject-specific — + each shoulder's BLUP", cex.main = 0.88)
    abline(0, 1, col = "grey30", lwd = 1.4, lty = 2)
    legend("topleft", bty = "n", cex = 0.74,
           legend = sprintf("concordance R² = %.3f", r2_conc(p$obs, p[[w]])))
  }
  mtext(sprintf("%s   ·   the gap between the panels is what the random curves buy", ttl),
        outer = TRUE, line = 0.2, cex = 0.72, font = 2)
}

draw_individual <- function(e, ri, ttl) {
  ind <- e$ind[order(e$ind$rmse), ]
  nw  <- min(4L, nrow(ind))
  worst <- rev(tail(ind$shoulder, nw))
  op <- par(mfrow = c(1, 1 + nw), mar = c(4.2, 4.4, 3.0, 0.8), mgp = c(2.5, 0.7, 0),
            oma = c(0, 0, 2.8, 0)); on.exit(par(op), add = TRUE)

  cols <- ifelse(ind$flagged, "#b2182b", BLUE)
  bp <- barplot(ind$rmse, horiz = TRUE, col = cols, border = NA, xlab = "RMSE (°)",
                main = "per-shoulder fit", cex.main = 0.9, cex.axis = 0.75)
  # The threshold is a robust z on log(rmse), so it maps back to a curved cut on
  # the raw scale; drawn where it actually falls rather than at a round number.
  lg <- log(ind$rmse); m <- mad(lg, constant = 1)
  if (m > 0) { thr <- exp(median(lg) + 3*1.4826*m)
    if (thr < max(ind$rmse)*3) abline(v = thr, lty = 2, col = "#b2182b", lwd = 1.3) }
  mtext(sprintf("worst/median = %.2f×  ·  %d of %d flagged (robust z > 3)",
                max(ind$ratio_to_median), sum(ind$flagged), nrow(ind)),
        side = 3, line = 0.1, cex = 0.6, col = "grey35")

  p <- e$pred_ref
  for (sh in worst) {
    q <- p[p$ID == sh, ]
    rg <- range(c(q$obs, q$pred_ss))
    plot(q$TIME, q$obs, pch = 16, cex = 0.5, col = adjustcolor("grey35", 0.5), ylim = rg,
         xlab = "elevation (°)", ylab = "angle (°)",
         main = sprintf("%s\nRMSE %.2f°", substr(sh, 1, 22), ind$rmse[ind$shoulder == sh]),
         cex.main = 0.72, cex.axis = 0.75)
    q <- q[order(q$TIME), ]; lines(q$TIME, q$pred_ss, col = BLUE, lwd = 2)
  }
  # These are step 1's shoulders, i.e. the rows K was chosen on — the reference
  # condition for a compared cell, not every shoulder in it. Said explicitly,
  # because "per-shoulder fit" otherwise reads as all of them.
  mtext(sprintf("%s   ·   worst-fitted of the %d shoulders K was chosen on (%s), subject-specific prediction, K = %d",
                ttl, nrow(e$ind), e$k_basis, e$K),
        outer = TRUE, line = 0.4, cex = 0.68, font = 2)
}

draw_diagnostics <- function(e, ri, ttl) {
  p <- e$pred_full; r <- p$resid_ss
  op <- par(mfrow = c(2, 2), mar = c(4.0, 4.2, 2.8, 1.0), mgp = c(2.4, 0.7, 0),
            oma = c(0, 0, 2.8, 0)); on.exit(par(op), add = TRUE)
  plot(p$pred_ss, r, pch = 16, cex = 0.3, col = adjustcolor("grey30", 0.25),
       xlab = "fitted (°)", ylab = "residual (°)", main = "residuals vs fitted", cex.main = 0.9)
  abline(h = 0, col = BLUE, lty = 2)
  qqnorm(r, pch = 16, cex = 0.3, col = adjustcolor("grey30", 0.25), main = "normal Q-Q",
         cex.main = 0.9); qqline(r, col = BLUE, lwd = 1.6)
  hist(r, breaks = 60, col = "grey85", border = "white", xlab = "residual (°)",
       main = "residual distribution", cex.main = 0.9)
  plot(p$TIME, r, pch = 16, cex = 0.3, col = adjustcolor(COL[p$cond], 0.25),
       xlab = "elevation (°)", ylab = "residual (°)", main = "residuals vs elevation",
       cex.main = 0.9)
  abline(h = 0, col = BLUE, lty = 2)
  mtext(sprintf("%s   ·   sd %.3f°  skew %.2f  excess kurtosis %.2f",
                ttl, ri$resid_sd, ri$skew, ri$kurtosis),
        outer = TRUE, line = 0.4, cex = 0.72, font = 2)
}

# Same wipe as the per-cell folders, for the plate. 00_forest and
# 00_significance_map are drawn only when some cell has a difference curve, so a
# step-1-only run skips them entirely and leaves whatever was there before —
# which is how a forest from a superseded basis survived three re-runs looking
# current. Clearing first means a plate figure on disk was drawn by THIS run or
# it is not there at all. CSVs and 00_SUMMARY.md are left alone.
unlink(list.files(OUT, pattern = "^(00_|planche_).*\\.(png|pdf)$", full.names = TRUE))

cat(sprintf("\n---- per-cell figures (%d cells, %d cores) ----\n", length(CUR), FIGCORES))
one_cell <- function(kk) {
  e <- CUR[[kk]]
  pr <- strsplit(kk, "|", fixed = TRUE)[[1]]
  jt <- pr[1]; mo <- pr[2]; dof <- as.integer(pr[3])
  ri <- mrow(jt, mo, dof); if (is.null(ri)) return("")
  dir <- file.path(CELL_DIR, cell_slug(jt, mo, dof))
  ttl <- gsub("\n", " — ", cell_title(jt, mo, dof))
  # Clear the folder before redrawing. Each figure below is conditional on the
  # stage its cell has reached, so a figure drawn under earlier settings — or at
  # a stage since invalidated — would otherwise sit here looking current. That
  # is not hypothetical: a 03_difference left over from a superseded basis is
  # what made two runs of the same cell disagree. What is on disk should always
  # be what the cache currently supports.
  unlink(list.files(dir, pattern = "\\.(png|pdf)$", full.names = TRUE))

  # step 1 is enough for these two — the df selection and who it fits badly
  if (!is.null(e$bic))
    figout("01_choosing_k", function() draw_choosing_k(e, ri, ttl), w = 1250, h = 620, dir = dir)
  if (!is.null(e$ind) && nrow(e$ind))
    figout("05_individual_fit", function() draw_individual(e, ri, ttl), w = 1500, h = 560, dir = dir)
  # the rest need step 2's fitted curves and predictions
  if (!is.null(e$curves))
    figout("02_population_curves", function() draw_population(e, ri, ttl), w = 1000, h = 760, dir = dir)
  # and 03 additionally needs step 3's test to annotate
  if (isTRUE(e$compare) && !is.null(e$dif) && isTRUE(ri$has_step3))
    figout("03_difference", function() draw_difference(e, ri, ttl), w = 1000, h = 760, dir = dir)
  if (!is.null(e$pred_full)) {
    figout("04_observed_vs_predicted", function() draw_obs_pred(e, ri, ttl), w = 1250, h = 660, dir = dir)
    figout("06_diagnostics", function() draw_diagnostics(e, ri, ttl), w = 1100, h = 900, dir = dir)
  }
  # returned, not printed: forked workers interleave their output, so the lines
  # are collected and written in cache order once they are all back
  sprintf("  %-52s K=%2d %s", cell_slug(jt, mo, dof), e$K,
          if (isTRUE(e$compare)) sprintf("p_bonf %.1e %s", ri$p_bonf, pstar(ri)) else "single")
}
lines_out <- par_do(names(CUR), one_cell, "cell figures")
for (z in unlist(lines_out)) if (nzchar(z)) cat(z, "\n", sep = "")
cat(sprintf("  %d cell folders under %s/\n", sum(nzchar(unlist(lines_out))), CELL_DIR))

# =============================================================================
# PLATE — planches, one per motion
# =============================================================================
annot <- function(ri, e) {
  u <- par("usr"); dx <- u[2]-u[1]; dy <- u[4]-u[3]
  if (is.null(ri) || nrow(ri) != 1) return(invisible())
  if (identical(ri$mode, "compare")) {
    text(u[2]-0.03*dx, u[4]-0.04*dy, pstar(ri), adj = c(1,1), font = 2, cex = 1.2, col = "grey10")
    ln <- c(sprintf("%d/%d sh · %d/%d st", ri$n_ex, ri$n_in, ri$st_ex, ri$st_in),
            sprintf("K = %d%s%s", ri$K_star,
                    if (isTRUE(ri$bic_edge)) " (edge)" else "",
                    if (isTRUE(ri$overparam)) " (not identifiable)" else ""),
            # NA until step 2 has fitted this cell — say "not fitted yet" rather
            # than printing "|Δ| NA°"
            if (is.finite(ri$diff_mean_abs))
              sprintf("|Δ| %.0f° (max %.0f°)", ri$diff_mean_abs, ri$diff_max_abs)
            else "not fitted yet")
    for (k in seq_along(ln))
      text(u[1]+0.03*dx, u[3]+0.05*dy + (length(ln)-k)*0.08*dy, ln[k],
           adj = c(0,0), cex = 0.60, col = "grey30")
  } else if (identical(ri$mode, "single")) {
    text(u[1]+0.03*dx, u[3]+0.05*dy,
         sprintf("single group · K = %d", ri$K_star), adj = c(0,0), cex = 0.62, col = "grey30")
  } else {
    text(mean(u[1:2]), mean(u[3:4]), "not modelled", cex = 0.8, col = "grey55")
  }
}

planche_key <- function(diff = FALSE) {
  l1 <- if (diff)
    "blue line = EX VIVO − IN VIVO, band = 95% CI  ·  squares along the bottom = band excludes 0.  Note the sign: in vivo is the reference, so this is mirrored against v2."
  else
    "orange = ex vivo, green = in vivo  ·  line = population curve, band = 95% CI.  No knots are drawn."
  l2 <- "stars = LIKELIHOOD-RATIO TEST of H0 'no ex-vivo departure' against H1, Bonferroni-adjusted across every compared cell on the plate: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant."
  # The support rule is a constant, so say what it actually is rather than
  # hardcoding a number the footer would go on claiming after MIN_UNITS_X moved.
  supp <- if (MIN_UNITS_X <= 1)
    "Curves span the full abscissa — every x any shoulder reaches is fitted, so the thinnest stretches rest on a single shoulder and the 95% band is what shows it."
  else
    sprintf("Curves are drawn only where ≥%d shoulders have data.", MIN_UNITS_X)
  l3 <- paste("bottom-left: shoulders ex/in · studies ex/in  |  K, this cell's own spline df chosen by BIC  |  mean and max |Δ| in degrees. ", supp)
  mtext(l1, side = 1, outer = TRUE, line = 0.0, cex = 0.55, col = "grey25")
  mtext(l2, side = 1, outer = TRUE, line = 0.9, cex = 0.55, col = "grey25")
  mtext(l3, side = 1, outer = TRUE, line = 1.8, cex = 0.55, col = "grey25")
}

cat(sprintf("\n---- planches (%d motions x 2, %d cores) ----\n", length(MOTION_ORDER), FIGCORES))
one_planche <- function(mo) {
  key_of <- function(jt, k) cell_key(jt, mo, k)
  xall <- rng_or(lapply(JOINT_ORDER, function(jt) lapply(1:3, function(k) {
    e <- CUR[[key_of(jt,k)]]; if (is.null(e)) NULL else range(e$curves$x) })), c(0, 150))
  xall <- xall + c(-1, 1)*0.03*diff(xall)

  figout(paste0("planche_", slug(mo)), function() {
    op <- par(mfrow = c(4, 3), mar = c(3.4, 3.4, 2.4, 0.8), mgp = c(2.1, 0.7, 0),
              oma = c(3.6, 0, 2.2, 0)); on.exit(par(op), add = TRUE)
    for (jt in JOINT_ORDER) {
      yl <- rng_or(lapply(1:3, function(k) { e <- CUR[[key_of(jt,k)]]
        if (is.null(e)) NULL else range(e$raw$Y) }), c(-1, 1))
      for (k in 1:3) {
        e <- CUR[[key_of(jt,k)]]; ri <- mrow(jt, mo, k)
        ttl <- sprintf("%s — %s", capj(jt), DOF_LEGEND[[jt]][k])
        # e$curves is NULL until step 2 reaches this cell — draw the empty panel
        # rather than dying, so a partial sweep still produces a readable plate
        if (is.null(e) || is.null(e$curves)) {
          plot.new(); title(main = ttl, cex.main = 0.78); annot(ri, e); next }
        plot(e$raw$TIME, e$raw$Y, col = adjustcolor(COL[as.character(e$raw$cond)], 0.16),
             pch = 16, cex = 0.3, xlim = xall, ylim = yl, xlab = "elevation (°)",
             ylab = "angle (°)", main = ttl, cex.main = 0.78, cex.axis = 0.8)
        for (cc in unique(e$curves$cond)) { g <- e$curves[e$curves$cond == cc, ]
          polygon(c(g$x, rev(g$x)), c(g$lo, rev(g$hi)), col = adjustcolor(COL[cc], 0.28), border = NA)
          lines(g$x, g$fit, col = COL[cc], lwd = 2.2) }
        annot(ri, e)
      }
    }
    mtext(sprintf("%s — random spline, per-cell df by BIC, D diagonal", mo),
          outer = TRUE, cex = 0.85, font = 2)
    planche_key()
  }, w = 1300, h = 1800, res = 135)

  figout(paste0("planche_diff_", slug(mo)), function() {
    op <- par(mfrow = c(4, 3), mar = c(3.4, 3.4, 2.4, 0.8), mgp = c(2.1, 0.7, 0),
              oma = c(3.6, 0, 2.2, 0)); on.exit(par(op), add = TRUE)
    for (jt in JOINT_ORDER) {
      yl <- rng_or(c(0, unlist(lapply(1:3, function(k) { e <- CUR[[key_of(jt,k)]]
        if (is.null(e) || is.null(e$dif)) NULL else c(e$dif$lo, e$dif$hi) }))), c(-1, 1))
      for (k in 1:3) {
        e <- CUR[[key_of(jt,k)]]; ri <- mrow(jt, mo, k)
        ttl <- sprintf("%s — %s", capj(jt), DOF_LEGEND[[jt]][k])
        if (is.null(e) || is.null(e$dif)) { plot.new(); title(main = ttl, cex.main = 0.78)
          annot(ri, e); next }
        z <- e$dif; sig <- (z$lo > 0) | (z$hi < 0)
        plot(z$x, z$fit, type = "n", xlim = xall, ylim = yl, xlab = "elevation (°)",
             ylab = "ex vivo − in vivo (°)", main = ttl, cex.main = 0.78, cex.axis = 0.8)
        polygon(c(z$x, rev(z$x)), c(z$lo, rev(z$hi)), col = adjustcolor(BLUE, 0.25), border = NA)
        abline(h = 0, lty = 2, col = "grey40"); lines(z$x, z$fit, lwd = 2.4, col = BLUE)
        if (any(sig)) points(z$x[sig], rep(yl[1], sum(sig)), pch = 15, cex = 0.3, col = BLUE)
        annot(ri, e)
      }
    }
    mtext(sprintf("%s — ex vivo − in vivo difference, 95%% CI (overlap only)", mo),
          outer = TRUE, cex = 0.85, font = 2)
    planche_key(diff = TRUE)
  }, w = 1300, h = 1800, res = 135)
  sprintf("  saved planche_%s and planche_diff_%s", slug(mo), slug(mo))
}
for (z in unlist(par_do(MOTION_ORDER, one_planche, "planches"))) cat(z, "\n", sep = "")

# =============================================================================
# 00_k_map — how K was chosen, for the whole plate at once
#
# One tile per cell: K as a number on a sequential fill, and inside it the
# cell's own BIC-vs-K curve as a sparkline with its minimum marked. A tile whose
# minimum sits on an edge of the grid, or that lost K to non-convergence, is
# hatched — those are not selections and must not read as ones.
# =============================================================================
cat("\n---- 00_k_map ----\n")
kmap_rows <- expand.grid(dof = 1:3, joint = JOINT_ORDER, stringsAsFactors = FALSE)
kmap_rows <- kmap_rows[, c("joint", "dof")]
KS_RANGE  <- rng_or(master$K_star, range(KS))
kfill <- function(k) {
  if (is.na(k)) return("grey93")
  ramp <- colorRampPalette(c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))(100)
  if (diff(KS_RANGE) == 0) return(ramp[50])
  ramp[1 + round(99 * (k - KS_RANGE[1]) / diff(KS_RANGE))]
}

figout("00_k_map", function() {
  op <- par(mar = c(1.0, 15.5, 6.0, 1.0)); on.exit(par(op), add = TRUE)
  nr <- nrow(kmap_rows); nc <- length(MOTION_ORDER)
  plot(NA, xlim = c(0, nc), ylim = c(nr, 0), axes = FALSE, xlab = "", ylab = "")
  for (ci in seq_len(nc)) {
    mo <- MOTION_ORDER[ci]
    # MOTION_SHORT, not the full name wrapped: wrapping "internal-external
    # rotation 0 degree-abducted" and its 90-degree twin renders both as
    # "internal-external" with the distinguishing half clipped, so two different
    # columns end up labelled identically.
    text(ci - 0.5, -0.15, MOTION_SHORT[[mo]], xpd = NA,
         adj = c(0.5, 1), cex = 0.72, font = 3, col = "grey20")
    for (rj in seq_len(nr)) {
      jt <- kmap_rows$joint[rj]; dof <- kmap_rows$dof[rj]
      ri <- mrow(jt, mo, dof); e <- CUR[[cell_key(jt, mo, dof)]]
      x0 <- ci - 1; y0 <- rj - 1
      kst <- if (is.null(ri)) NA else ri$K_star
      rect(x0 + 0.03, y0 + 0.05, x0 + 0.97, y0 + 0.95, col = kfill(kst),
           border = "white", lwd = 1.4)
      if (is.na(kst)) {
        text(x0 + 0.5, y0 + 0.5, if (is.null(ri)) "—" else substr(ri$note, 1, 14),
             cex = 0.5, col = "grey55"); next
      }
      # the sparkline: this cell's own BIC curve, scaled inside the tile
      if (!is.null(e$bic)) {
        b <- e$bic[e$bic$converged, ]
        if (nrow(b) > 1) {
          sx <- x0 + 0.10 + 0.80 * (b$K - min(KS)) / max(1, diff(range(KS)))
          sy <- y0 + 0.78 - 0.34 * (b$BIC - min(b$BIC)) / max(1e-9, diff(range(b$BIC)))
          lines(sx, sy, col = adjustcolor("grey15", 0.55), lwd = 1.1)
          im <- which.min(b$BIC)
          points(sx[im], sy[im], pch = 21, bg = "#b2182b", col = "white", cex = 0.7, lwd = 0.8)
        }
      }
      # High K fills dark, so both the number and the flag have to flip to a
      # light ink or they vanish into the tile — and the flag is precisely what
      # a dark tile most needs to show.
      dark <- !is.na(kst) && kst > mean(KS_RANGE)
      flag_col <- if (dark) "#fddbc7" else "#b2182b"
      text(x0 + 0.5, y0 + 0.36, kst, cex = 1.25, font = 2,
           col = if (dark) "white" else "grey10")
      # not a selection: say so on the tile
      if (isTRUE(ri$bic_edge) || isTRUE(ri$overparam) ||
          (!is.na(ri$n_converged) && ri$n_converged < length(KS))) {
        mk <- c(if (isTRUE(ri$bic_edge)) "edge",
                if (isTRUE(ri$overparam)) sprintf("%d var/%d sh", ri$K_star + 1, ri$n_shoulders_ref),
                if (!is.na(ri$n_converged) && ri$n_converged < length(KS))
                  sprintf("%d/%d conv", ri$n_converged, length(KS)))
        text(x0 + 0.5, y0 + 0.88, paste(mk, collapse = " · "), cex = 0.46,
             col = flag_col, font = 2)
      }
      if (identical(ri$mode, "single"))
        text(x0 + 0.5, y0 + 0.62, "single", cex = 0.46,
             col = if (dark) "grey85" else "grey35", font = 3)
    }
  }
  for (rj in seq_len(nr)) {
    jt <- kmap_rows$joint[rj]; dof <- kmap_rows$dof[rj]
    mtext(sprintf("%s · DoF%d — %s", capj(jt), dof, DOF_LEGEND[[jt]][dof]),
          side = 2, at = rj - 0.5, las = 1, line = 0.4, cex = 0.60, col = "grey20")
    if (dof == 1 && rj > 1) segments(0, rj-1, nc, rj-1, col = "grey60", lwd = 1.2, xpd = NA)
  }
  mtext("How K was chosen — the spline df each cell selected by BIC on its own reduced model",
        side = 3, line = 4.6, cex = 0.92, font = 2, adj = 0)
  # One long mtext runs off the right edge and is silently clipped, so the key is
  # broken into lines that fit the device.
  key <- c(sprintf("fill and number = the retained K over the grid K = %d..%d   ·   sparkline = that cell's own BIC-vs-K curve, red dot = its minimum",
                   min(KS), max(KS)),
           "'edge' = the minimum sits on a grid endpoint — BIC reported the edge rather than turning   ·   'n/20 conv' = the other K did not converge",
           "'N var/M sh' = the random part asks for N variances from M shoulders, so K is NOT identifiable: it describes the flexibility the optimiser bore, not the trajectory")
  for (i in seq_along(key))
    mtext(key[i], side = 3, line = 3.4 - 0.9*(i-1), cex = 0.56, col = "grey30", adj = 0)
}, w = 1450, h = 1180, res = 135)
message("  saved 00_k_map")

# =============================================================================
# 00_forest and the significance maps
# =============================================================================
forest_frame <- function(master) {
  fo <- master[master$mode == "compare" & !is.na(master$diff_mean_signed), ]
  if (!nrow(fo)) return(list(fo = fo, lab = character(0), hgt = 600))
  fo$mkey <- match(fo$motion, MOTION_ORDER); fo$jkey <- match(fo$joint, JOINT_ORDER)
  fo <- fo[order(fo$mkey, fo$jkey, fo$dof), ]
  fo <- fo[rev(seq_len(nrow(fo))), ]                     # first motion at the top
  lab <- sprintf("%s · DoF%d — %s", capj(fo$joint), fo$dof,
                 mapply(function(j,k) DOF_LEGEND[[j]][k], fo$joint, fo$dof))
  list(fo = fo, lab = lab, hgt = max(600, 150 + 30*nrow(fo)))
}
FR <- forest_frame(master); fo <- FR$fo; lab <- FR$lab; hgt <- FR$hgt

if (nrow(fo)) {
  cat("\n---- 00_forest ----\n")
  figout("00_forest", function() {
    # THREE REAL PANELS, not one plot with fat margins. v2's TODO recorded this
    # plate as unreadable — right-hand columns overprinting, motion labels
    # clipped off the left edge — and widening `mar` only makes it worse: the
    # margins eat the plot region, so the columns (placed at fractions of the
    # DATA range) crowd together even harder. sigmap_plate already learned this:
    # "a panel cannot clip". Same fix here — labels, forest and numbers each get
    # their own panel with their own coordinate system.
    op <- par(oma = c(4.6, 0, 3.6, 0)); on.exit(par(op), add = TRUE)
    layout(matrix(1:3, nrow = 1), widths = c(3.5, 5.0, 1.6))
    ylim <- c(0.5, nrow(fo) + 0.5)
    band <- function(xl) for (mi in unique(fo$mkey)) { rr <- which(fo$mkey == mi)
      if (mi %% 2 == 0) rect(xl[1], min(rr)-0.5, xl[2], max(rr)+0.5,
                             col = adjustcolor("grey85", 0.35), border = NA, xpd = NA) }

    # ---- panel 1: motion group + cell label ----
    par(mar = c(3.4, 0.4, 0.6, 0.4))
    plot(NA, xlim = c(0, 1), ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    band(c(0, 1))
    for (mi in unique(fo$mkey)) { rr <- which(fo$mkey == mi)
      text(0.02, mean(rr), paste(strwrap(MOTION_ORDER[mi], 14), collapse = "\n"),
           adj = c(0, 0.5), cex = 0.62, col = "grey25", font = 3) }
    text(1, seq_len(nrow(fo)), lab, adj = c(1, 0.5), cex = 0.66)

    # ---- panel 2: the forest ----
    par(mar = c(3.4, 0.3, 0.6, 0.3))
    xr <- range(0, fo$diff_lo_curve, fo$diff_hi_curve, na.rm = TRUE)
    xr <- xr + c(-1, 1)*0.04*diff(xr)
    plot(NA, xlim = xr, ylim = ylim, yaxt = "n", xlab = "", ylab = "", cex.axis = 0.8)
    band(xr)
    abline(v = 0, lty = 2, col = "grey45")
    for (i in seq_len(nrow(fo))) { r <- fo[i, ]
      sig <- !is.na(r$p_bonf) && r$p_bonf < 0.05
      col <- if (sig) BLUE else "grey55"
      segments(r$diff_lo_curve, i, r$diff_hi_curve, i, lwd = 7,
               col = adjustcolor(col, 0.28), lend = 1)          # span across elevation
      segments(r$diff_mean_lo, i, r$diff_mean_hi, i, lwd = 2.6, col = col, lend = 1)
      points(r$diff_mean_signed, i, pch = 21, bg = col, col = "white", cex = 1.05, lwd = 1.1)
    }
    mtext("ex vivo − in vivo (°)", side = 1, line = 2.1, cex = 0.72)

    # ---- panel 3: the numeric columns ----
    par(mar = c(3.4, 0.4, 0.6, 0.4))
    plot(NA, xlim = c(0, 1), ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    band(c(0, 1))
    cx <- c(0.16, 0.50, 0.78, 0.94)
    hd <- c("mean |Δ|", "elevation", "K", "LRT")
    for (q in seq_along(cx))
      text(cx[q], nrow(fo) + 1.0, hd[q], xpd = NA, cex = 0.64, font = 2, adj = c(0.5, 0.5))
    for (i in seq_len(nrow(fo))) { r <- fo[i, ]
      text(cx[1], i, sprintf("%.1f°", r$diff_mean_abs), cex = 0.64, adj = c(0.5, 0.5))
      text(cx[2], i, sprintf("%.0f–%.0f°", r$x_lo, r$x_hi), cex = 0.62, adj = c(0.5, 0.5), col = "grey30")
      text(cx[3], i, sprintf("%d%s", r$K_star, if (isTRUE(r$bic_edge)) "*" else ""),
           cex = 0.64, adj = c(0.5, 0.5), col = if (isTRUE(r$bic_edge)) "#b2182b" else "grey20")
      text(cx[4], i, stars(r$p_bonf), cex = 0.70, adj = c(0.5, 0.5), font = 2,
           col = if (!is.na(r$p_bonf) && r$p_bonf < 0.05) "grey10" else "grey55") }

    mtext(sprintf("All %d compared cells — per-cell spline df by BIC, LRT with Bonferroni", nrow(fo)),
          outer = TRUE, side = 3, line = 1.0, cex = 0.95, font = 2)
    # Legend and footer belong to the plate, not a panel: placed by device
    # position, as in sigmap_plate.
    omb <- par("omd")[3]; xc <- grconvertX(0.5, from = "ndc", to = "user")
    legend(xc, grconvertY(omb*0.62, from = "ndc", to = "user"),
           xjust = 0.5, yjust = 0.5, xpd = NA, horiz = TRUE, bty = "n", cex = 0.64,
           lwd = c(7, 2.6), col = c(adjustcolor(BLUE, 0.28), BLUE),
           legend = c("range across elevation", "95% CI of the mean difference"))
    text(xc, grconvertY(omb*0.22, from = "ndc", to = "user"),
         "K* = the BIC minimum sat on a grid endpoint   ·   stars = likelihood-ratio test, Bonferroni-adjusted across all compared cells",
         cex = 0.58, col = "grey30", adj = c(0.5, 0.5), xpd = NA)
  }, w = 1500, h = hgt, res = 135)
  message("  saved 00_forest")

  cat("\n---- 00_significance_map ----\n")
  segs <- lapply(seq_len(nrow(fo)), function(i)
    sig_segments(CUR[[cell_key(fo$joint[i], fo$motion[i], fo$dof[i])]]))
  figout("00_significance_map", function() {
    op <- par(mar = c(4.4, 22, 3.2, 12), mgp = c(2.4, 0.7, 0)); on.exit(par(op), add = TRUE)
    xr <- range(fo$x_lo, fo$x_hi, na.rm = TRUE); xr <- xr + c(-1,1)*0.03*diff(xr)
    plot(NA, xlim = xr, ylim = c(0.5, nrow(fo)+0.5), yaxt = "n",
         xlab = "thoracohumeral elevation (°)", ylab = "",
         main = "Where in the movement do the conditions differ?", cex.main = 0.95)
    for (mi in unique(fo$mkey)) { rr <- which(fo$mkey == mi)
      if (mi %% 2 == 0) rect(xr[1], min(rr)-0.5, xr[2], max(rr)+0.5,
                             col = adjustcolor("grey85", 0.35), border = NA)
      mtext(paste(strwrap(MOTION_ORDER[mi], 14), collapse = "\n"), side = 2,
            at = mean(rr), line = 20.5, las = 1, cex = 0.6, col = "grey25", font = 3) }
    abline(v = 0, lty = 3, col = "grey60")
    for (i in seq_len(nrow(fo))) { r <- fo[i, ]
      segments(r$x_lo, i, r$x_hi, i, lwd = 6, col = adjustcolor("grey60", 0.35), lend = 1)
      sg <- segs[[i]]
      if (!is.null(sg) && nrow(sg))
        segments(sg$lo, i, sg$hi, i, lwd = 6, col = adjustcolor(BLUE, 0.9), lend = 1) }
    axis(2, at = seq_len(nrow(fo)), labels = lab, las = 1, cex.axis = 0.66, tick = FALSE, line = -0.6)
    u <- par("usr"); c1 <- u[2] + 0.04*diff(u[1:2]); c2 <- u[2] + 0.20*diff(u[1:2])
    text(c1, nrow(fo)+1.4, "mean |Δ|",   xpd = NA, cex = 0.66, font = 2, adj = 0)
    text(c2, nrow(fo)+1.4, "% of range", xpd = NA, cex = 0.66, font = 2, adj = 0)
    for (i in seq_len(nrow(fo))) { r <- fo[i, ]
      frac <- 100 * r$sig_frac_x
      text(c1, i, sprintf("%.1f°", r$diff_mean_abs), xpd = NA, cex = 0.64, adj = 0)
      text(c2, i, sprintf("%.0f%%", frac), xpd = NA, cex = 0.64, adj = 0,
           col = if (frac >= 50) "grey10" else "grey50") }
    legend("bottomright", inset = c(0, -0.085), xpd = NA, bty = "n", cex = 0.62, horiz = TRUE,
           lwd = 6, col = c(adjustcolor("grey60", 0.35), adjustcolor(BLUE, 0.9)),
           legend = c("fitted range (x-overlap)", "95% band excludes zero"))
  }, w = 1500, h = hgt, res = 135)
  message("  saved 00_significance_map")
}

# ---- DoF-grouped significance plates ----------------------------------------
# 00_significance_map orders its rows motion -> joint -> DoF, which answers
# "where in THIS movement does THIS cell differ". It cannot answer the other
# question: does the difference on one DoF depend on which plane the arm moves
# in? — because the three planes of a given joint x DoF end up scattered across
# the plate. These plates invert the nesting: rows are joint x DoF, and the
# motions are the INNER axis, so the planes sit adjacent on one shared abscissa.
HDR_H <- 2.10   # header band, tall enough for joint / DoF / description stacked
GAP_H <- 0.55
TOP_H <- 1.45

sigmap_frame <- function(groups, motions) {
  rownames(groups) <- NULL
  groups$y_top <- NA_real_; groups$y_bot <- NA_real_; groups$y_hdr <- NA_real_
  y <- 0; rows <- list(); hdr <- numeric(nrow(groups))
  for (g in seq_len(nrow(groups))) {
    hdr[g] <- y + HDR_H/2; y <- y + HDR_H
    for (mo in motions) {
      rows[[length(rows)+1]] <- data.frame(g = g, joint = groups$joint[g],
        dof = groups$dof[g], motion = mo, y = y + 0.5, stringsAsFactors = FALSE)
      y <- y + 1
    }
    groups$y_top[g] <- hdr[g] - HDR_H/2; groups$y_bot[g] <- y
    y <- y + GAP_H
  }
  groups$y_hdr <- hdr
  list(rows = do.call(rbind, rows), groups = groups, ytot = y - GAP_H + 0.25)
}

sigmap_plate <- function(groups, motions, base, title, subtitle) {
  FR <- sigmap_frame(groups, motions); RW <- FR$rows; GR <- FR$groups
  ylim <- c(FR$ytot, -TOP_H)
  MR <- lapply(seq_len(nrow(RW)), function(i) mrow(RW$joint[i], RW$motion[i], RW$dof[i]))
  cells <- lapply(seq_len(nrow(RW)), function(i)
    CUR[[cell_key(RW$joint[i], RW$motion[i], RW$dof[i])]])
  gv <- function(i, f) { m <- MR[[i]]; if (is.null(m)) NA else m[[f]] }

  xr <- rng_or(lapply(seq_len(nrow(RW)), function(i)
    c(gv(i,"x_lo"), gv(i,"x_hi"), gv(i,"x_ex_lo"), gv(i,"x_ex_hi"),
      gv(i,"x_in_lo"), gv(i,"x_in_hi"))), c(0, 150))
  xr <- xr + c(-1, 1)*0.03*diff(xr)
  dr <- rng_or(c(0, unlist(lapply(seq_len(nrow(RW)), function(i)
    c(gv(i,"diff_mean_lo"), gv(i,"diff_mean_hi"), gv(i,"diff_mean_signed"))))), c(-1, 1))
  if (diff(dr) == 0) dr <- dr + c(-1, 1)
  dr <- dr + c(-1, 1)*0.10*diff(dr)

  shade <- function(xl) for (g in seq_len(nrow(GR))) if (g %% 2 == 0)
    rect(xl[1], GR$y_top[g], xl[2], GR$y_bot[g],
         col = adjustcolor("grey85", 0.35), border = NA, xpd = NA)

  figout(base, function() {
    op <- par(oma = c(5.0, 0, 3.4, 0)); on.exit(par(op), add = TRUE)
    layout(matrix(1:4, nrow = 1), widths = c(3.0, 5.2, 1.7, 1.4))
    MAR <- c(3.6, 0.3, 0.6, 0.3)

    # ---- panel 1: labels ----
    par(mar = c(MAR[1], 0.4, MAR[3], 0.6))
    plot(NA, xlim = c(0,1), ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    shade(c(0,1))
    for (g in seq_len(nrow(GR))) {
      if (g > 1) segments(0.02, GR$y_top[g] - GAP_H/2, 1.0, GR$y_top[g] - GAP_H/2,
                          col = "grey80", lwd = 0.8, xpd = NA)
      yc <- GR$y_hdr[g]
      text(0.55, yc - 0.60, capj(GR$joint[g]), adj = c(0.5, 0.5), cex = 0.80, font = 2)
      text(0.55, yc, sprintf("DoF%d", GR$dof[g]), adj = c(0.5, 0.5), cex = 0.72, col = "grey25")
      text(0.55, yc + 0.60, DOF_LEGEND[[GR$joint[g]]][GR$dof[g]], adj = c(0.5, 0.5),
           cex = 0.68, font = 3, col = "grey20")
    }
    for (i in seq_len(nrow(RW)))
      text(0.98, RW$y[i], MOTION_SHORT[[RW$motion[i]]], adj = c(1, 0.5), cex = 0.70, col = "grey30")

    # ---- panel 2: the significance map ----
    par(mar = MAR)
    plot(NA, xlim = xr, ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    shade(xr)
    axis(1, cex.axis = 0.72, mgp = c(2.1, 0.55, 0))
    mtext("thoracohumeral elevation (°)", side = 1, line = 1.9, cex = 0.66)
    abline(v = 0, lty = 3, col = "grey65")
    for (i in seq_len(nrow(RW))) {
      y <- RW$y[i]
      # coverage of each condition ON ITS OWN, above and below the fitted band.
      # Two stripes rather than one: the point is precisely where they DON'T
      # coincide — the overlap below is their intersection.
      if (is.finite(gv(i,"x_ex_lo")))
        rect(gv(i,"x_ex_lo"), y-0.40, gv(i,"x_ex_hi"), y-0.17,
             col = adjustcolor(COL[["ex vivo"]], 0.55), border = NA)
      if (is.finite(gv(i,"x_in_lo")))
        rect(gv(i,"x_in_lo"), y+0.17, gv(i,"x_in_hi"), y+0.40,
             col = adjustcolor(COL[["in vivo"]], 0.55), border = NA)
      if (identical(gv(i,"mode"), "compare")) {
        rect(gv(i,"x_lo"), y-0.13, gv(i,"x_hi"), y+0.13,
             col = adjustcolor("grey55", 0.40), border = NA)
        sg <- sig_segments(cells[[i]])
        if (!is.null(sg) && nrow(sg))
          rect(sg$lo, y-0.13, sg$hi, y+0.13, col = adjustcolor(BLUE, 0.95), border = NA)
      } else {
        # Why this row is empty, on the face of the figure. Scapular plane has no
        # compared cell anywhere in the sweep and IER 90° has no ex-vivo shoulder
        # in any joint; both are findings, not gaps in the plotting.
        ne <- gv(i,"n_ex"); ni <- gv(i,"n_in")
        msg <- if (is.na(ne)) "not modelled"
               else if (ne == 0) "no ex-vivo data"
               else if (ni == 0) "no in-vivo data"
               else sprintf("%d ex / %d in shoulders (need ≥%d each) — not compared",
                            ne, ni, MIN_PER_COND)
        text(mean(xr), y, msg, cex = 0.60, font = 3, col = "grey45")
      }
    }
    # ---- panel 3: mean difference ----
    par(mar = c(MAR[1], 0.5, MAR[3], 0.5))
    plot(NA, xlim = dr, ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    shade(dr)
    tk <- axTicks(1)
    # a shared difference axis is what makes the groups comparable, but it also
    # crushes a 10° effect next to a 60° one — the gridlines are what keep the
    # small bars readable without giving each group its own scale
    abline(v = tk, col = adjustcolor("grey70", 0.45), lwd = 0.6)
    axis(1, at = tk, cex.axis = 0.72, mgp = c(2.1, 0.55, 0))
    mtext("ex vivo − in vivo (°)", side = 1, line = 1.9, cex = 0.62)
    abline(v = 0, lty = 2, col = "grey45")
    for (i in seq_len(nrow(RW))) {
      if (!identical(gv(i,"mode"), "compare")) next
      y <- RW$y[i]; d <- gv(i,"diff_mean_signed")
      if (!is.finite(d)) next          # step 2 has not fitted this cell yet
      p <- gv(i,"p_bonf"); sig <- !is.na(p) && p < 0.05
      cl <- if (sig) BLUE else "grey55"
      rect(0, y-0.30, d, y+0.30, col = adjustcolor(cl, 0.55), border = NA)
      segments(gv(i,"diff_mean_lo"), y, gv(i,"diff_mean_hi"), y, lwd = 1.6, col = cl)
      points(d, y, pch = 21, bg = cl, col = "white", cex = 0.62, lwd = 0.8)
    }

    # ---- panel 4: inference — K and the Bonferroni LRT ----
    par(mar = c(MAR[1], 0.5, MAR[3], 0.4))
    plot(NA, xlim = c(0,1), ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    shade(c(0,1))
    text(0.30, -0.95, "K",   cex = 0.66, font = 2, adj = c(0.5, 0.5))
    text(0.72, -0.95, "LRT", cex = 0.66, font = 2, adj = c(0.5, 0.5))
    for (i in seq_len(nrow(RW))) {
      if (!identical(gv(i,"mode"), "compare")) next
      y <- RW$y[i]; p <- gv(i,"p_bonf"); ke <- gv(i,"K_star")
      text(0.30, y, sprintf("%d%s", ke, if (isTRUE(gv(i,"bic_edge"))) "*" else ""),
           cex = 0.68, adj = c(0.5, 0.5),
           col = if (isTRUE(gv(i,"bic_edge"))) "#b2182b" else "grey20")
      text(0.72, y, stars(p), cex = 0.72, font = 2, adj = c(0.5, 0.5),
           col = if (!is.na(p) && p < 0.05) "grey10" else "grey55")
    }

    mtext(title, outer = TRUE, side = 3, line = 1.0, cex = 0.92, font = 2)
    mtext(subtitle, outer = TRUE, side = 3, line = 0.0, cex = 0.64, col = "grey30")
    # The legend and the footer belong to the plate, not to any one panel, so both
    # are placed by DEVICE position and drawn with xpd = NA from wherever we are.
    # Anchored to fractions of omd[3] — the NDC band the outer bottom margin
    # occupies — rather than of the whole device, because the plate's height
    # varies while that margin stays five lines.
    omb <- par("omd")[3]
    xc  <- grconvertX(0.5, from = "ndc", to = "user")
    legend(xc, grconvertY(omb*0.60, from = "ndc", to = "user"),
           xjust = 0.5, yjust = 0.5, xpd = NA, horiz = TRUE, bty = "n",
           cex = 0.62, border = NA,
           fill = c(adjustcolor(COL[["ex vivo"]], 0.55), adjustcolor(COL[["in vivo"]], 0.55),
                    adjustcolor("grey55", 0.40), adjustcolor(BLUE, 0.95)),
           legend = c("ex-vivo coverage", "in-vivo coverage", "fitted overlap (the model)",
                      "95% band on the difference excludes 0"))
    text(xc, grconvertY(omb*0.20, from = "ndc", to = "user"),
         paste("bars = 95% CI of the mean difference  ·  K = this cell's own spline df,",
               "chosen by BIC (* = the minimum sat on a grid endpoint)  ·  stars =",
               "likelihood-ratio test, Bonferroni-adjusted across the whole plate"),
         cex = 0.58, col = "grey30", adj = c(0.5, 0.5), xpd = NA)
  }, w = 1500, h = max(430, 150 + 33*FR$ytot), res = 135)
  message("  saved ", base)
}

cat("\n---- DoF-grouped significance plates ----\n")
ELEV_SUB <- paste("rows = joint × DoF, three elevation planes each, one shared abscissa.",
                  "Scapular plane has NO compared cell in the whole sweep (≤ 2 ex-vivo",
                  "shoulders everywhere).")
grid_elev <- expand.grid(dof = 1:3, joint = JOINT_ORDER, stringsAsFactors = FALSE)
grid_elev <- grid_elev[, c("joint", "dof")]
sigmap_plate(grid_elev, ELEV_MOTIONS, "00_sigmap_elevation",
             "Where the conditions differ, by DoF across the three elevation planes", ELEV_SUB)
# IER 90° is kept deliberately: it has zero ex-vivo shoulders in every joint, and
# a row saying so is more honest than the motion's absence.
sigmap_plate(grid_elev, ROT_MOTIONS, "00_sigmap_rotation",
             "Where the conditions differ, internal-external rotation at 0° and 90°",
             paste("IER 90° has NO ex-vivo shoulder in any joint, so no cell there can be",
                   "compared — the empty rows are the finding."))

# =============================================================================
# 00_SUMMARY.md
# =============================================================================
cm <- master[master$mode == "compare" & !is.na(master$p_lrt), ]
cm <- cm[order(-cm$diff_mean_abs), ]
nsig_b <- sum(cm$p_bonf < 0.05); nsig_h <- sum(cm$p_bh < 0.05)
edge   <- master[!is.na(master$bic_edge) & master$bic_edge, ]

md <- c(
  "# Generalized sweep v3 — per-cell spline df, likelihood-ratio test, Bonferroni",
  "",
  sprintf("[Iteration 06](../06_natural_spline_hlme/)'s inference applied to all %d", nrow(master)),
  "joint × humeral motion × degree-of-freedom cells.",
  "",
  "```",
  "1  K by BIC on the REDUCED model (no condition terms), K = 1..20,",
  "   fitted on the REFERENCE condition's own rows.  Random part mirrors K.",
  "1a per-shoulder RMSE, robust flag, marginal vs subject-specific R2.",
  "2  full model at the SAME K, all shoulders, in vivo as the reference.",
  "3  LR = 2(LL1 - LL0) ~ chi2(K+1), Bonferroni across every compared cell.",
  "4  no Wald tests anywhere — the difference curve is the reporting object.",
  "```",
  "",
  "## What came out",
  "",
  sprintf("| cells | %d |", nrow(master)),
  "| --- | ---: |",
  sprintf("| compared | %d |", sum(master$mode == "compare")),
  sprintf("| single condition | %d |", sum(master$mode == "single")),
  sprintf("| skipped | %d |", sum(master$mode == "skip")),
  sprintf("| tested | %d |", nrow(cm)),
  sprintf("| significant, Bonferroni | **%d** |", nsig_b),
  sprintf("| significant, BH | %d |", nsig_h),
  sprintf("| K retained, range | %d – %d |",
          min(master$K_star, na.rm = TRUE), max(master$K_star, na.rm = TRUE)),
  sprintf("| cells whose BIC minimum sat on a grid edge | %d |", nrow(edge)),
  "",
  "**Bonferroni is the headline** — Melanie asked for the conservative adjustment,",
  "applied to the LRT p-values of every compared cell at once. BH is reported",
  "beside it only so v3 can be read against v2.",
  "",
  "## The difference curve is the reporting object",
  "",
  "In vivo is the reference, so every difference is **ex vivo − in vivo** — mirrored",
  "against v2 and iterations 03–05. `03_difference` in each cell folder is the figure",
  "that carries the finding.",
  "",
  "## Figures",
  "",
  "| file | what it shows |",
  "| --- | --- |",
  "| `00_k_map` | the spline df every cell chose, with its own BIC-vs-K curve inside each tile |",
  "| `00_forest` | every compared cell on one difference axis |",
  "| `00_significance_map` | where in the movement each cell differs |",
  "| `00_sigmap_elevation`, `00_sigmap_rotation` | the same, grouped by DoF so the planes sit adjacent |",
  "| `planche_<motion>` | population curves, rows = joints, cols = DoF |",
  "| `planche_diff_<motion>` | the ex vivo − in vivo difference, same layout |",
  "| `cells/<cell>/01_choosing_k` | that cell's df selection, shown rather than asserted |",
  "| `cells/<cell>/02_population_curves` … `06_diagnostics` | the per-cell drill-down |",
  "",
  "## Compared cells, largest mean |Δ| first",
  "",
  "| joint | motion | DoF | K | shoulders ex/in | mean \\|Δ\\| | LR | df | p | Bonferroni | |",
  "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |")
for (i in seq_len(nrow(cm))) { r <- cm[i, ]
  md <- c(md, sprintf("| %s | %s | %d | %d%s | %d/%d | %.1f° | %.1f | %d | %.1e | %.1e | %s |",
    capj(r$joint), r$motion, r$dof, r$K_star, if (isTRUE(r$bic_edge)) "\\*" else "",
    r$n_ex, r$n_in, r$diff_mean_abs, r$LR, r$lrt_df, r$p_lrt, r$p_bonf, stars(r$p_bonf))) }
md <- c(md, "",
  "\\* the BIC minimum sat on an endpoint of the K grid — the edge was reported, not a choice.",
  "",
  "## Standing caveat",
  "",
  "**Condition is perfectly confounded with source study** — no study measured both.",
  "Whatever the LRT says, the result is descriptive and not causal. See",
  "[`../../docs/STATISTICAL_METHODS_REVIEW.md`](../../docs/STATISTICAL_METHODS_REVIEW.md).")

writeLines(md, file.path(OUT, "00_SUMMARY.md"))
cat(sprintf("\nDONE. See %s/00_SUMMARY.md\n", OUT))
