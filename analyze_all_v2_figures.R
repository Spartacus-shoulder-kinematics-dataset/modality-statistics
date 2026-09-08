# =============================================================================
# GENERALIZED SWEEP v2 — THE FIGURE HALF
#
# Draws every plate from figures/generalized_v2/cache/v2_fit.rds. Fits nothing,
# reads no raw data, takes seconds. Edit a colour or a label and re-run this;
# the 12-minute sweep in analyze_all_v2_fit.R only has to happen when the MODEL
# changes.
#
# Run (from the repository root):  Rscript analyze_all_v2_figures.R
#
# Outputs (figures/generalized_v2/):
#   planche_<motion>.png/.pdf       population curves, rows = joints, cols = DoF
#   planche_diff_<motion>.png/.pdf  the in-vivo - ex-vivo difference, same layout
#   00_forest.png/.pdf              every compared cell on one difference axis
#   00_significance_map.png/.pdf    WHERE in the movement each cell differs
#   00_SUMMARY.md                   human-readable index
# =============================================================================

source("analyze_all_v2_common.R")

if (!file.exists(CACHE))
  stop("No fit cache at ", CACHE, "\n  Run:  Rscript analyze_all_v2_fit.R  (~12 min) first.")
cache <- readRDS(CACHE)
check_meta(cache$meta)
CUR <- cache$cells; master <- cache$master
cat(sprintf("loaded %d fitted cells from %s (fitted %s)\n", length(CUR), CACHE,
            format(cache$meta$fitted_at, "%Y-%m-%d %H:%M")))

# ---- planche annotation -----------------------------------------------------
annot <- function(ri) {
  u <- par("usr"); dx <- u[2]-u[1]; dy <- u[4]-u[3]
  if (nrow(ri) != 1) return(invisible())
  if (identical(ri$mode, "compare")) {
    text(u[2]-0.03*dx, u[4]-0.04*dy, stars(ri$joint_p_fdr), adj = c(1,1), font = 2,
         cex = 1.2, col = "grey10")
    # glyph strip with the coefficient index printed under each dot, so there is no
  # guessing which gamma a dot refers to. Left to right: gamma_0, gamma_1, ...
  # drawn as POINTS, not text glyphs: the filled/hollow circle characters have
  # no glyph in the legacy PDF font set and would silently disappear there.
  gl <- strsplit(ri$sym_diff, "")[[1]]; ng <- length(gl); stp <- 0.052*dx
  x0 <- u[2] - 0.03*dx - (ng-1)*stp
  for (q in seq_len(ng)) {
    points(x0 + (q-1)*stp, u[4]-0.175*dy, pch = 21, cex = 0.95, lwd = 1.2,
           bg = if (gl[q] == "●") BLUE else "white", col = BLUE)
    text(x0 + (q-1)*stp, u[4]-0.235*dy, q-1, adj = c(0.5,1), cex = 0.5, col = "grey45")
  }
  text(x0 - 0.035*dx, u[4]-0.19*dy, expression(gamma), adj = c(1,0.5), cex = 0.6, col = "grey45")
    ln <- c(sprintf("%d/%d sh · %d/%d st", ri$n_ex, ri$n_in, ri$st_ex, ri$st_in),
            sprintf("Wald ref %d/%d · diff %d/%d", ri$wald_pass_ref, ri$wald_total_ref,
                    ri$wald_pass_diff, ri$wald_total_diff),
            sprintf("|Δ| %.0f° (max %.0f°)", ri$diff_mean_abs, ri$diff_max_abs))
    for (k in seq_along(ln))
      text(u[1]+0.03*dx, u[3]+0.05*dy + (length(ln)-k)*0.08*dy, ln[k],
           adj = c(0,0), cex = 0.60, col = "grey30")
  } else if (identical(ri$mode, "single")) {
    text(u[1]+0.03*dx, u[3]+0.05*dy, "single group", adj = c(0,0), cex = 0.62, col = "grey30")
  } else {
    text(mean(u[1:2]), mean(u[3:4]), "not modelled", cex = 0.8, col = "grey55")
  }
}
knot_marks <- function(kn, curves = NULL, ycol = "fit") {
  u <- par("usr"); ik <- kn[kn > u[1] & kn < u[2]]
  if (length(ik)) abline(v = ik, lty = 2, lwd = 0.8, col = adjustcolor("grey35", 0.5))
  if (!is.null(curves) && length(ik)) {
    for (g in unique(curves$cond)) { cc <- curves[curves$cond == g, ]; cc <- cc[order(cc$x), ]
      ins <- ik >= min(cc$x) & ik <= max(cc$x); if (!any(ins)) next
      points(ik[ins], approx(cc$x, cc[[ycol]], xout = ik[ins])$y, pch = 21,
             bg = COL[[as.character(g)]], col = "white", cex = 0.7, lwd = 0.9) }
  }
}

# The x-intervals where the 95% band on the difference excludes zero — the solid
# bars of every significance map. A run-length encoding of the pointwise flag, so
# a cell whose significance is broken into several stretches draws as several
# bars rather than one bar spanning the gaps. NULL for a cell that was not compared.
sig_segments <- function(e) {
  if (is.null(e) || is.null(e$dif)) return(NULL)
  z <- e$dif; sig <- (z$lo > 0) | (z$hi < 0)
  r <- rle(sig); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  k <- which(r$values)
  if (!length(k)) return(data.frame(lo = numeric(0), hi = numeric(0)))
  data.frame(lo = z$x[starts[k]], hi = z$x[ends[k]])
}

# Footer key — every glyph on a planche, spelled out. Drawn in the outer bottom
# margin so it applies to all twelve panels.
planche_key <- function(diff = FALSE) {
  l1 <- if (diff)
    "blue line = in vivo − ex vivo, band = 95% CI  ·  squares along the bottom = band excludes 0  ·  dashed verticals = interior knots, dots = value at each knot"
  else
    "orange = ex vivo, green = in vivo  ·  line = population curve, band = 95% CI  ·  dashed verticals = interior knots, dots = fitted value at each knot"
  l2 <- "top-right: one dot per in-vivo difference coefficient, read LEFT TO RIGHT as gamma_0, gamma_1, ... — the index is printed under each dot.  gamma_0 = constant level shift (condin vivo);  gamma_k = ns_k:condin vivo, a shape term."
  l2b <- "filled = that coefficient passed |coef/se| >= 1.96, hollow = did not.    stars = JOINT test on ALL difference coefficients at once, BH-FDR: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant"
  l3 <- "bottom-left: shoulders ex/in · studies ex/in  |  Wald ref/diff counts  |  mean and max |Δ| in degrees.  Curves are drawn only where ≥2 shoulders have data."
  mtext(l1,  side = 1, outer = TRUE, line = 0.0, cex = 0.55, col = "grey25")
  mtext(l2,  side = 1, outer = TRUE, line = 0.9, cex = 0.55, col = "grey25")
  mtext(l2b, side = 1, outer = TRUE, line = 1.8, cex = 0.55, col = "grey25")
  mtext(l3,  side = 1, outer = TRUE, line = 2.7, cex = 0.55, col = "grey25")
}

# ---- planches: population curves, and differences ---------------------------
for (mo in MOTION_ORDER) {
  key_of <- function(jt, k) cell_key(jt, mo, k)

  # one abscissa for the whole plate, from the SUPPORTED ranges of its cells
  xall <- range(unlist(lapply(JOINT_ORDER, function(jt) lapply(1:3, function(k) {
    e <- CUR[[key_of(jt,k)]]; if (is.null(e)) NULL else range(e$curves$x) }))), na.rm = TRUE)
  if (!all(is.finite(xall))) xall <- c(0, 1)
  xall <- xall + c(-1, 1)*0.03*diff(xall)

  figout(paste0("planche_", slug(mo)), function() {
    op <- par(mfrow = c(4, 3), mar = c(3.4, 3.4, 2.4, 0.8), mgp = c(2.1, 0.7, 0),
              oma = c(4.4, 0, 2.2, 0)); on.exit(par(op), add = TRUE)
    for (jt in JOINT_ORDER) {
      yl <- range(unlist(lapply(1:3, function(k) { e <- CUR[[key_of(jt,k)]]
        if (is.null(e)) NULL else range(e$raw$Y) })), na.rm = TRUE)
      if (!all(is.finite(yl))) yl <- c(-1, 1)
      for (k in 1:3) {
        e <- CUR[[key_of(jt,k)]]; ri <- master[master$joint==jt & master$motion==mo & master$dof==k, ]
        ttl <- sprintf("%s — %s", capj(jt), DOF_LEGEND[[jt]][k])
        if (is.null(e)) { plot.new(); title(main = ttl, cex.main = 0.78); annot(ri); next }
        plot(e$raw$TIME, e$raw$Y, col = adjustcolor(COL[as.character(e$raw$cond)], 0.16),
             pch = 16, cex = 0.3, xlim = xall, ylim = yl, xlab = "elevation (°)", ylab = "angle (°)",
             main = ttl, cex.main = 0.78, cex.axis = 0.8)
        for (cc in unique(e$curves$cond)) { g <- e$curves[e$curves$cond == cc, ]
          polygon(c(g$x, rev(g$x)), c(g$lo, rev(g$hi)), col = adjustcolor(COL[cc], 0.28), border = NA)
          lines(g$x, g$fit, col = COL[cc], lwd = 2.2) }
        knot_marks(e$knots, e$curves); annot(ri)
      }
    }
    mtext(sprintf("%s — random spline (ns df=%d, D diagonal), trimmed to overlap", mo, KDF),
          outer = TRUE, cex = 0.85, font = 2)
    planche_key()
  }, w = 1300, h = 1800, res = 135)
  message("  saved planche_", slug(mo), ".png")

  figout(paste0("planche_diff_", slug(mo)), function() {
    op <- par(mfrow = c(4, 3), mar = c(3.4, 3.4, 2.4, 0.8), mgp = c(2.1, 0.7, 0),
              oma = c(4.4, 0, 2.2, 0)); on.exit(par(op), add = TRUE)
    for (jt in JOINT_ORDER) {
      yl <- range(0, unlist(lapply(1:3, function(k) { e <- CUR[[key_of(jt,k)]]
        if (is.null(e) || is.null(e$dif)) NULL else c(e$dif$lo, e$dif$hi) })), na.rm = TRUE)
      if (!all(is.finite(yl))) yl <- c(-1, 1)
      for (k in 1:3) {
        e <- CUR[[key_of(jt,k)]]; ri <- master[master$joint==jt & master$motion==mo & master$dof==k, ]
        ttl <- sprintf("%s — %s", capj(jt), DOF_LEGEND[[jt]][k])
        if (is.null(e) || is.null(e$dif)) { plot.new(); title(main = ttl, cex.main = 0.78)
          annot(ri); next }
        z <- e$dif; sig <- (z$lo > 0) | (z$hi < 0)
        plot(z$x, z$fit, type = "n", xlim = xall, ylim = yl, xlab = "elevation (°)",
             ylab = "in vivo − ex vivo (°)", main = ttl, cex.main = 0.78, cex.axis = 0.8)
        polygon(c(z$x, rev(z$x)), c(z$lo, rev(z$hi)), col = adjustcolor(BLUE, 0.25), border = NA)
        abline(h = 0, lty = 2, col = "grey40"); lines(z$x, z$fit, lwd = 2.4, col = BLUE)
        if (any(sig)) points(z$x[sig], rep(yl[1], sum(sig)), pch = 15, cex = 0.3, col = BLUE)
        ik <- e$knots[e$knots > min(z$x) & e$knots < max(z$x)]
        if (length(ik)) { abline(v = ik, lty = 2, lwd = 0.8, col = adjustcolor("grey35", 0.5))
          points(ik, approx(z$x, z$fit, xout = ik)$y, pch = 21, bg = BLUE, col = "white",
                 cex = 0.7, lwd = 0.9) }
        annot(ri)
      }
    }
    mtext(sprintf("%s — in vivo − ex vivo difference, 95%% CI (overlap only)", mo),
          outer = TRUE, cex = 0.85, font = 2)
    planche_key(diff = TRUE)
  }, w = 1300, h = 1800, res = 135)
  message("  saved planche_diff_", slug(mo), ".png")
}

# ---- the row ordering shared by the forest and the significance map ----------
# Both plates show the same cells, in the same order, with the same labels and
# the same height. This used to be three loose globals assigned halfway through
# the forest block, which is why the significance map could not be drawn on its
# own. It is a pure function of master, so it is now one.
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

# ---- 00_forest.png — every compared cell on one axis -------------------------
# One row per cell. Two horizontal spans, deliberately different weights:
#   wide pale band = the RANGE the difference curve covers across elevation, i.e.
#                    how much the difference itself varies with arm position
#   narrow dark bar = the 95% CI of the MEAN difference — the inferential span
#   dot            = the mean signed difference
# A cell whose pale band is wide but whose dark bar straddles zero has a
# difference that changes sign or size with elevation and averages to nothing.
if (nrow(fo)) {
  figout("00_forest", function() {
    op <- par(mar = c(4.4, 22, 3.2, 22), mgp = c(2.4, 0.7, 0)); on.exit(par(op), add = TRUE)
    xr <- range(0, fo$diff_lo_curve, fo$diff_hi_curve, na.rm = TRUE)
    xr <- xr + c(-1, 1)*0.04*diff(xr)
    plot(NA, xlim = xr, ylim = c(0.5, nrow(fo)+0.5), yaxt = "n", xlab = "in vivo − ex vivo (°)",
         ylab = "", main = sprintf("All %d compared cells — random spline (ns df=%d, D diagonal), trimmed",
                                   nrow(fo), KDF), cex.main = 0.95)
    abline(v = 0, lty = 2, col = "grey45")
    # alternating motion bands, so the six movements are separable at a glance
    for (mi in unique(fo$mkey)) { rr <- which(fo$mkey == mi)
      if (mi %% 2 == 0) rect(xr[1], min(rr)-0.5, xr[2], max(rr)+0.5,
                             col = adjustcolor("grey85", 0.35), border = NA)
      mtext(sub(" ", "\n", MOTION_ORDER[mi]), side = 2, at = mean(rr), line = 19,
            las = 1, cex = 0.62, col = "grey25", font = 3) }
    for (i in seq_len(nrow(fo))) { r <- fo[i, ]
      sig <- !is.na(r$joint_p_fdr) && r$joint_p_fdr < 0.05
      col <- if (sig) BLUE else "grey55"
      segments(r$diff_lo_curve, i, r$diff_hi_curve, i, lwd = 7,
               col = adjustcolor(col, 0.28), lend = 1)          # span across elevation
      segments(r$diff_mean_lo, i, r$diff_mean_hi, i, lwd = 2.6, col = col, lend = 1)
      points(r$diff_mean_signed, i, pch = 21, bg = col, col = "white", cex = 1.05, lwd = 1.1)
    }
    axis(2, at = seq_len(nrow(fo)), labels = lab, las = 1, cex.axis = 0.66, tick = FALSE, line = -0.6)
    # right-hand columns: effect size, Wald count, significance
    u <- par("usr"); x1 <- u[2] + 0.03*diff(u[1:2])
    x2 <- u[2] + 0.14*diff(u[1:2]); x3 <- u[2] + 0.29*diff(u[1:2])
    x4 <- u[2] + 0.36*diff(u[1:2])
    text(x1, nrow(fo)+1.4, "mean |Δ|", xpd = NA, cex = 0.66, font = 2, adj = 0)
    text(x2, nrow(fo)+1.4, "elevation range", xpd = NA, cex = 0.66, font = 2, adj = 0)
    text(x3, nrow(fo)+1.4, "Wald", xpd = NA, cex = 0.66, font = 2, adj = 0)
    text(x4, nrow(fo)+1.4, "joint", xpd = NA, cex = 0.66, font = 2, adj = 0)
    for (i in seq_len(nrow(fo))) { r <- fo[i, ]
      text(x1, i, sprintf("%.1f°", r$diff_mean_abs), xpd = NA, cex = 0.64, adj = 0)
      text(x2, i, sprintf("%.0f–%.0f°", r$x_lo, r$x_hi), xpd = NA, cex = 0.64, adj = 0, col = "grey30")
      text(x3, i, sprintf("%d/%d", r$wald_pass_diff, r$wald_total_diff), xpd = NA, cex = 0.64, adj = 0)
      text(x4, i, stars(r$joint_p_fdr), xpd = NA, cex = 0.72, adj = 0, font = 2,
           col = if (!is.na(r$joint_p_fdr) && r$joint_p_fdr < 0.05) "grey10" else "grey55") }
    legend("bottomright", inset = c(0, -0.085), xpd = NA, bty = "n", cex = 0.62, horiz = TRUE,
           lwd = c(7, 2.6), col = c(adjustcolor(BLUE, 0.28), BLUE),
           legend = c("range across elevation", "95% CI of the mean"))
  }, w = 1500, h = hgt, res = 135)
  message("  saved 00_forest.png")
}

# ---- 00_significance_map.png — WHERE each cell differs -----------------------
# The forest answers "by how much"; this answers "over which elevations".
# One row per compared cell, x = elevation:
#   pale bar  = the range actually fitted (the x-overlap of the two conditions)
#   solid bar = the sub-ranges where the 95% band on the difference excludes zero
# A cell can have a large mean difference and almost no solid bar, which means the
# difference is real on average but nowhere separable from zero point by point.
if (nrow(fo)) {
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
    axis(2, at = seq_len(nrow(fo)), labels = lab, las = 1, cex.axis = 0.66,
         tick = FALSE, line = -0.6)
    u <- par("usr"); c1 <- u[2] + 0.04*diff(u[1:2]); c2 <- u[2] + 0.20*diff(u[1:2])
    text(c1, nrow(fo)+1.4, "mean |Δ|", xpd = NA, cex = 0.66, font = 2, adj = 0)
    text(c2, nrow(fo)+1.4, "% of range", xpd = NA, cex = 0.66, font = 2, adj = 0)
    for (i in seq_len(nrow(fo))) { r <- fo[i, ]; sg <- segs[[i]]
      frac <- 100 * r$sig_frac_x
      text(c1, i, sprintf("%.1f°", r$diff_mean_abs), xpd = NA, cex = 0.64, adj = 0)
      text(c2, i, sprintf("%.0f%%", frac), xpd = NA, cex = 0.64, adj = 0,
           col = if (frac >= 50) "grey10" else "grey50") }
    legend("bottomright", inset = c(0, -0.085), xpd = NA, bty = "n", cex = 0.62, horiz = TRUE,
           lwd = 6, col = c(adjustcolor("grey60", 0.35), adjustcolor(BLUE, 0.9)),
           legend = c("fitted range (x-overlap)", "95% band excludes zero"))
  }, w = 1500, h = hgt, res = 135)
  message("  saved 00_significance_map.png/.pdf")
}

# =============================================================================
# DoF-GROUPED SIGNIFICANCE MAPS
#
# 00_significance_map above orders its rows motion -> joint -> DoF, which answers
# "where in THIS movement does THIS cell differ". It cannot answer the other
# question: does the difference on one DoF depend on which plane the arm moves
# in? — because the three planes of a given joint x DoF end up scattered across
# the plate. These plates invert the nesting: rows are joint x DoF, and the
# motions are the INNER axis, three to a group, so the planes sit adjacent and
# on one shared abscissa.
#
# Four panels side by side, sharing one row geometry:
#   labels | significance map | mean difference | inference
# =============================================================================

# Row geometry. y runs DOWNWARD (ylim is reversed when plotting), one unit per
# motion row, a taller band per group for its three-line header, a gap between
# groups, and negative space above the first group for the column headings.
HDR_H <- 2.10   # header band, tall enough for joint / DoF / description stacked
GAP_H <- 0.55   # between one group's last row and the next group's header
TOP_H <- 1.45   # column headings above the first group

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

# Every panel is a real plot in the layout, never a margin: the motion labels on
# 00_significance_map are drawn with mtext at line = 20.5 against mar = 22 and are
# clipped off the device. A panel cannot clip.
sigmap_plate <- function(groups, motions, base, title, subtitle) {
  FR <- sigmap_frame(groups, motions); RW <- FR$rows; GR <- FR$groups
  ylim <- c(FR$ytot, -TOP_H)

  # every cell's master row, in row order; NA-padded so indexing is total
  mrow <- lapply(seq_len(nrow(RW)), function(i) {
    m <- master[master$joint == RW$joint[i] & master$motion == RW$motion[i] &
                master$dof == RW$dof[i], ]
    if (nrow(m) == 1) m else NULL })
  cells <- lapply(seq_len(nrow(RW)), function(i)
    CUR[[cell_key(RW$joint[i], RW$motion[i], RW$dof[i])]])
  gv <- function(i, f) { m <- mrow[[i]]; if (is.null(m)) NA else m[[f]] }

  # shared abscissa: every span any row draws, so groups are comparable end to end
  xr <- range(unlist(lapply(seq_len(nrow(RW)), function(i)
    c(gv(i,"x_lo"), gv(i,"x_hi"), gv(i,"x_ex_lo"), gv(i,"x_ex_hi"),
      gv(i,"x_in_lo"), gv(i,"x_in_hi")))), na.rm = TRUE)
  if (!all(is.finite(xr))) xr <- c(0, 150)
  xr <- xr + c(-1, 1)*0.03*diff(xr)
  # shared difference axis, always including zero
  dr <- range(0, unlist(lapply(seq_len(nrow(RW)), function(i)
    c(gv(i,"diff_mean_lo"), gv(i,"diff_mean_hi"), gv(i,"diff_mean_signed")))), na.rm = TRUE)
  if (!all(is.finite(dr)) || diff(dr) == 0) dr <- c(-1, 1)
  dr <- dr + c(-1, 1)*0.10*diff(dr)

  shade <- function(xl) for (g in seq_len(nrow(GR))) if (g %% 2 == 0)
    rect(xl[1], GR$y_top[g], xl[2], GR$y_bot[g],
         col = adjustcolor("grey85", 0.35), border = NA, xpd = NA)

  figout(base, function() {
    op <- par(oma = c(5.0, 0, 3.4, 0)); on.exit(par(op), add = TRUE)
    layout(matrix(1:4, nrow = 1), widths = c(3.0, 5.2, 1.7, 1.8))
    MAR <- c(3.6, 0.3, 0.6, 0.3)

    # ---- panel 1: labels ----------------------------------------------------
    # The three components of the old one-line label get one centred line each,
    # so joint names align with joint names and DoF with DoF down the whole
    # plate. The DoF is written ONCE per group, never repeated per motion row.
    par(mar = c(MAR[1], 0.4, MAR[3], 0.6))
    plot(NA, xlim = c(0,1), ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    shade(c(0,1))
    for (g in seq_len(nrow(GR))) {
      if (g > 1) segments(0.02, GR$y_top[g] - GAP_H/2, 1.0, GR$y_top[g] - GAP_H/2,
                          col = "grey80", lwd = 0.8, xpd = NA)
      yc <- GR$y_hdr[g]
      text(0.55, yc - 0.60, capj(GR$joint[g]),  adj = c(0.5, 0.5), cex = 0.80, font = 2)
      text(0.55, yc,        sprintf("DoF%d", GR$dof[g]), adj = c(0.5, 0.5),
           cex = 0.72, col = "grey25")
      text(0.55, yc + 0.60, DOF_LEGEND[[GR$joint[g]]][GR$dof[g]], adj = c(0.5, 0.5),
           cex = 0.68, font = 3, col = "grey20")
    }
    for (i in seq_len(nrow(RW)))
      text(0.98, RW$y[i], MOTION_SHORT[[RW$motion[i]]], adj = c(1, 0.5), cex = 0.70,
           col = "grey30")

    # ---- panel 2: the significance map --------------------------------------
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
    # ---- panel 3: mean difference, closely stacked ---------------------------
    par(mar = c(MAR[1], 0.5, MAR[3], 0.5))
    plot(NA, xlim = dr, ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    shade(dr)
    tk <- axTicks(1)
    # a shared difference axis is what makes the groups comparable, but it also
    # crushes a 10° effect next to a 60° one — the gridlines are what keep the
    # small bars readable without giving each group its own scale
    abline(v = tk, col = adjustcolor("grey70", 0.45), lwd = 0.6)
    axis(1, at = tk, cex.axis = 0.72, mgp = c(2.1, 0.55, 0))
    mtext("in vivo − ex vivo (°)", side = 1, line = 1.9, cex = 0.62)
    abline(v = 0, lty = 2, col = "grey45")
    for (i in seq_len(nrow(RW))) {
      if (!identical(gv(i,"mode"), "compare")) next
      y <- RW$y[i]; d <- gv(i,"diff_mean_signed")
      sig <- !is.na(gv(i,"joint_p_fdr")) && gv(i,"joint_p_fdr") < 0.05
      cl <- if (sig) BLUE else "grey55"
      rect(0, y-0.30, d, y+0.30, col = adjustcolor(cl, 0.55), border = NA)
      # the 95% CI of the MEAN — computed as xbar' V xbar in the fit, not by
      # averaging pointwise SEs, which would understate it badly
      segments(gv(i,"diff_mean_lo"), y, gv(i,"diff_mean_hi"), y, lwd = 1.6, col = cl)
      points(d, y, pch = 21, bg = cl, col = "white", cex = 0.62, lwd = 0.8)
    }

    # ---- panel 4: inference --------------------------------------------------
    par(mar = c(MAR[1], 0.5, MAR[3], 0.4))
    plot(NA, xlim = c(0,1), ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    shade(c(0,1))
    ng <- KDF + 1; dx <- seq(0.40, 0.94, length.out = ng)
    text(0.13, -0.95, "FDR", cex = 0.64, font = 2, adj = c(0.5, 0.5))
    text(mean(dx), -0.95, expression(gamma), cex = 0.70, font = 2, adj = c(0.5, 0.5))
    for (q in seq_len(ng))
      text(dx[q], -0.42, q-1, cex = 0.52, col = "grey45", adj = c(0.5, 0.5))
    for (i in seq_len(nrow(RW))) {
      if (!identical(gv(i,"mode"), "compare")) next
      y <- RW$y[i]; p <- gv(i,"joint_p_fdr")
      text(0.13, y, stars(p), cex = 0.72, font = 2, adj = c(0.5, 0.5),
           col = if (!is.na(p) && p < 0.05) "grey10" else "grey55")
      # drawn as POINTS, not the filled/hollow circle characters: those have no
      # glyph in the legacy PDF font set and vanish silently there
      s <- gv(i,"sym_diff"); if (is.na(s)) next
      gl <- strsplit(s, "")[[1]]
      for (q in seq_along(gl))
        points(dx[q], y, pch = 21, cex = 0.80, lwd = 1.1,
               bg = if (gl[q] == "●") BLUE else "white", col = BLUE)
    }

    mtext(title, outer = TRUE, side = 3, line = 1.0, cex = 0.92, font = 2)
    mtext(subtitle, outer = TRUE, side = 3, line = 0.0, cex = 0.64, col = "grey30")
    # The legend and the footer belong to the plate, not to any one panel, so both
    # are placed by DEVICE position and drawn with xpd = NA from wherever we are.
    #
    # Not via a par(fig=) overlay: setting `oma` in the same par() call silently
    # resets the figure region, so the overlay lands back inside the layout and
    # the legend prints over the axis. grconvert is immune to that. The anchors
    # are fractions of omd[3] — the NDC band the outer bottom margin occupies —
    # rather than of the whole device, because the plate's height varies sixfold
    # between the one-group POC and the twelve-group sweep while that margin
    # stays five lines, and a fixed device fraction would drift into the rows.
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
         paste("bars = 95% CI of the mean difference  ·  ●/○ = per-coefficient Wald",
               "|coef/se| ≥ 1.96, index printed above  ·  stars = JOINT test on all",
               "γ at once, BH-FDR"),
         cex = 0.58, col = "grey30", adj = c(0.5, 0.5), xpd = NA)
  }, w = 1500, h = max(430, 150 + 33*FR$ytot), res = 135)
  message("  saved ", base, ".png/.pdf")
}

ELEV_SUB <- paste("rows = joint × DoF, three elevation planes each, one shared abscissa.",
                  "Scapular plane has NO compared cell in the whole sweep (≤ 2 ex-vivo",
                  "shoulders everywhere).")

# POC — the single group the layout was designed against
sigmap_plate(data.frame(joint = "scapulothoracic", dof = 2, stringsAsFactors = FALSE),
             ELEV_MOTIONS, "00_sigmap_poc_scapulothoracic_dof2",
             "Scapulothoracic DoF2 across the three elevation planes", ELEV_SUB)

# the full elevation plate, 12 joint x DoF groups
grid_elev <- expand.grid(dof = 1:3, joint = JOINT_ORDER, stringsAsFactors = FALSE)
grid_elev <- grid_elev[, c("joint", "dof")]
sigmap_plate(grid_elev, ELEV_MOTIONS, "00_sigmap_elevation",
             "Where the conditions differ, by DoF across the three elevation planes",
             ELEV_SUB)

# the rotation plate. IER 90° is kept deliberately: it has zero ex-vivo shoulders
# in every joint, and a row saying so is more honest than the motion's absence.
sigmap_plate(grid_elev, ROT_MOTIONS, "00_sigmap_rotation",
             "Where the conditions differ, internal-external rotation at 0° and 90°",
             paste("IER 90° has NO ex-vivo shoulder in any joint, so no cell there can be",
                   "compared — the empty rows are the finding."))

# ---- summary ----------------------------------------------------------------
cm <- master[master$mode == "compare" & !is.na(master$joint_p_fdr), ]
cm <- cm[order(-cm$diff_mean_abs), ]
md <- c("# Generalized sweep v2 — random spline + hlme",
  "",
  sprintf("[Iteration 04](../04_natural_spline_hlme/)'s model applied to all %d", nrow(master)),
  "joint x motion x DoF combinations. Produced by",
  "[`../../analyze_all_v2.R`](../../analyze_all_v2.R).",
  "",
  "## The model",
  "",
  "```r",
  sprintf("hlme(fixed   = Y ~ ns(TIME, knots) * cond,      # knots: %d interior, per cell", KDF - 1),
  "     random  = ~ ns(TIME, knots),               # a random CURVE per shoulder",
  "     subject = \"IDnum\", ng = 1,",
  "     idiag   = TRUE)                            # D diagonal",
  "```",
  "",
  "For a shoulder $i$, observation $j$, in condition $c_i \\in \\{\\text{ex vivo},\\ \\text{in vivo}\\}$,",
  "at elevation $x_{ij}$:",
  "",
  "$$",
  "\\begin{aligned}",
  sprintf("y_{ij} =\\;& \\underbrace{\\beta_0 + \\sum_{k=1}^{%d}\\beta_k B_k(x_{ij})}_{\\text{ex-vivo reference}}", KDF),
  sprintf(" + \\underbrace{\\Big[\\gamma_0 + \\sum_{k=1}^{%d}\\gamma_k B_k(x_{ij})\\Big]\\mathbb{1}(c_i=\\text{in vivo})}_{\\text{in-vivo departure}} \\\\[4pt]", KDF),
  sprintf("&+ \\underbrace{b_{i0} + \\sum_{k=1}^{%d} b_{ik} B_k(x_{ij})}_{\\text{this shoulder's own curve}} + \\varepsilon_{ij},", KDF),
  sprintf("\\qquad \\mathbf{b}_i \\sim \\mathcal{N}_{%d}(\\mathbf{0}, D),", KDF + 1),
  "\\quad \\varepsilon_{ij}\\sim\\mathcal{N}(0,\\sigma_\\varepsilon^2).",
  "\\end{aligned}",
  "$$",
  "",
  sprintf("$D$ is **diagonal** (`idiag = TRUE`): %d variances, no covariances. Per cell that is", KDF + 1),
  sprintf("**%d fixed effects** (intercept, %d basis terms, the level shift, %d interactions),", 2*KDF + 2, KDF, KDF),
  sprintf("%d random-effect variances and 1 residual sd.", KDF + 1),
  "",
  "$\\mathbb{1}(c_i = \\text{in vivo})$ is 1 on in-vivo rows and 0 otherwise, so",
  sprintf("$\\gamma_0$ shifts the curve and $\\gamma_1 \\dots \\gamma_{%d}$ bend it; all $\\gamma = 0$ would mean the", KDF),
  "conditions are indistinguishable. That is exactly the **joint test** tabulated below.",
  "",
  "Cells with only one condition are fitted without the $\\gamma$ block and reported",
  "descriptively.",
  "",
  "**Per cell:** thin to <=100 points/shoulder, trim to the x-overlap, build the basis from",
  "that cell's own quantiles with boundary knots from its untrimmed range. Curves are drawn",
  sprintf("only where at least %d shoulders have data.", MIN_UNITS_X),
  "",
  "## Results",
  "",
  sprintf("- **compare** (both conditions): %d", sum(master$mode=="compare")),
  sprintf("- **single** (one condition):    %d", sum(master$mode=="single")),
  sprintf("- **skipped** (too few units):   %d", sum(master$mode=="skip")),
  sprintf("- **did not converge**:          %d", sum(!is.na(master$converged) & !master$converged)),
  sprintf("- of the compared, **%d** reach BH-FDR p < 0.05 on the JOINT difference test",
          sum(cm$joint_p_fdr < 0.05, na.rm = TRUE)),
  "",
  "> **Exploratory, not causal.** Condition is perfectly confounded with source study",
  "> (no study measured both), so a difference mixes muscle activation with protocol,",
  "> measurement and population. The real replication is the number of **studies**",
  "> (`st ex/in`), not rows or shoulders. Read the effect size in degrees.",
  "",
  "> The **joint** test (all difference coefficients at once) is the one to read, not the",
  "> per-coefficient Wald counts: the basis coefficients are correlated, so which one",
  "> carries the signal shifts with the knots while the joint statistic does not.",
  "",
  "## Reading the planches",
  "",
  "Each motion has two plates — `planche_<motion>.png/.pdf` (population curves) and",
  "`planche_diff_<motion>.png/.pdf` (the in-vivo minus ex-vivo difference). Rows are the",
  "four joints proximal to distal, columns the three rotational DoF, and **all twelve",
  "panels share one abscissa** so cells are directly comparable.",
  "",
  "| glyph | meaning |",
  "| --- | --- |",
  sprintf("| **●** (filled dot, top right) | that in-vivo difference coefficient passed \\|coef/se\\| >= %.2f |", WALD_CRIT),
  "| **○** (hollow dot, top right) | that coefficient did **not** pass |",
  "| **\\*\\*\\*  \\*\\*  \\*  ns** (top right) | the **joint** chi2 test on *all* difference coefficients at once, BH-FDR adjusted: `***` p<0.001, `**` p<0.01, `*` p<0.05, `ns` not significant |",
  "| small dot **on a curve** | the fitted value at an interior knot |",
  "| **dashed vertical** | an interior knot — where the piecewise cubics join |",
  "| orange / green | ex vivo / in vivo; the band is the pointwise 95% CI |",
  "| blue band, difference plate | pointwise 95% CI of the difference; **squares along the bottom** mark where it excludes zero |",
  "| bottom-left text | shoulders ex/in, studies ex/in, Wald counts ref/diff, mean and max \\|delta\\| in degrees |",
  "",
  "### Which dot is which coefficient",
  "",
  sprintf("The strip has **%d dots, read LEFT TO RIGHT**; on the figure the index is printed", KDF + 1),
  "directly under each dot. The order is fixed:",
  "",
  "| position | symbol | model term | meaning |",
  "| ---: | :---: | :--- | :--- |",
  "| **1st** dot | $\\gamma_0$ | `condin vivo` | the **constant level shift** between conditions |",
  sapply(1:KDF, function(k) sprintf(
    "| **%s** dot | $\\gamma_%d$ | `ns%d:condin vivo` | shape term on basis function $B_%d$ |",
    c("2nd","3rd","4th","5th","6th")[k], k, k, k)),
  "",
  "So `○●○○●` reads: level shift **failed**, $\\gamma_1$ **passed**, $\\gamma_2$ failed,",
  "$\\gamma_3$ failed, $\\gamma_4$ **passed**.",
  "",
  "Every strip can be checked against **`00_wald_tests.csv`**, which has one row per",
  "coefficient with its `term`, `wald` and `passed`.",
  "",
  "The ● / ○ strip and the stars answer **different questions**, and they routinely",
  "disagree: the strip says which individual coefficients are separable from zero, the",
  "stars say whether the difference is nonzero *as a whole*. Because the basis",
  "coefficients are correlated, a cell can show `○○○○○` and still be strongly",
  "significant jointly. **Trust the stars, use the strip only as texture.**",
  "",
  "| joint | motion | DoF | sh ex/in | st ex/in | mean abs diff | joint chi2 | joint p (FDR) | shape p (FDR) | Wald ref | Wald diff | symbols |",
  "| --- | --- | ---: | :---: | :---: | ---: | ---: | ---: | ---: | :---: | :---: | :--- |")
for (i in seq_len(nrow(cm))) { r <- cm[i, ]
  md <- c(md, sprintf("| %s | %s | %d | %d/%d | %d/%d | %.1f° | %.1f | %.2g %s | %.2g | %d/%d | %d/%d | `%s` |",
    r$joint, r$motion, r$dof, r$n_ex, r$n_in, r$st_ex, r$st_in, r$diff_mean_abs,
    r$joint_chi2, r$joint_p_fdr, stars(r$joint_p_fdr), r$shape_p_fdr,
    r$wald_pass_ref, r$wald_total_ref, r$wald_pass_diff, r$wald_total_diff, r$sym_diff)) }
writeLines(md, file.path(OUT, "00_SUMMARY.md"))

cat(sprintf("\nDONE. See %s/00_SUMMARY.md\n", OUT))
