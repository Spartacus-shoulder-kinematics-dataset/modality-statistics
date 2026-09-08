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
# Run:  python3 prepare_monolix_data.py
#       Rscript analyze_all_v2.R
# =============================================================================

suppressMessages({ library(lcmm); library(splines); library(dplyr) })

LONG  <- "spartacus_angles_long.csv"
OUT   <- file.path("figures", "generalized_v2")
KDF          <- 4      # ns() columns; 3 interior knots
CAP          <- 100    # max points per shoulder
MIN_PER_COND <- 3      # >=3 shoulders in BOTH conditions -> compare
MIN_TOTAL    <- 4      # >=4 shoulders total              -> at least single
WALD_CRIT    <- 1.96
MIN_UNITS_X  <- 2      # a curve is only drawn where >=2 shoulders have data
NPROC <- max(1L, min(8L, parallel::detectCores() - 1L))
# V2_SUMMARY_ONLY=1 rebuilds 00_SUMMARY.md from the existing master CSV,
# skipping the ~12 min of fitting. Everything else is unchanged.
SUMMARY_ONLY <- nzchar(Sys.getenv("V2_SUMMARY_ONLY"))
COL   <- c("ex vivo" = "#d95f02", "in vivo" = "#1b9e77")
BLUE  <- "#377eb8"

JOINT_ORDER <- c("sternoclavicular", "acromioclavicular", "scapulothoracic", "glenohumeral")
DOF_LEGEND <- list(
  glenohumeral      = c("Plane of elevation", "Elevation(-)/Depression(+)", "Int(+)/ext(-) rotation"),
  scapulothoracic   = c("Protraction(+)/retraction(-)", "Medial(+)/lateral(-) rot", "Posterior(+)/anterior(-) tilt"),
  acromioclavicular = c("Protraction(+)/retraction(-)", "Medial(+)/lateral(-) rot", "Posterior(+)/anterior(-) tilt"),
  sternoclavicular  = c("Protraction(+)/retraction(-)", "Depression(+)/elevation(-)", "Backwards(+)/forward(-) rot"))
MOTION_ORDER <- c("frontal plane elevation", "scapular plane elevation", "sagittal plane elevation",
                  "horizontal flexion", "internal-external rotation 0 degree-abducted",
                  "internal-external rotation 90 degree-abducted")

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
slug  <- function(s) gsub("(^_|_$)", "", gsub("[^a-z0-9]+", "_", tolower(s)))
capj  <- function(j) paste0(toupper(substring(j, 1, 1)), substring(j, 2))
stars <- function(p) if (is.na(p)) "" else if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else "ns"
# Render every figure TWICE: a high-resolution PNG and a vector PDF.
# The draw argument is a FUNCTION, not a braced expression, so it can be called
# once per device — a promise would only evaluate on the first device and the
# second file would come out blank.
FIG_SCALE <- 2.2          # pixel multiplier; res scales with it so text stays proportional
figout <- function(base, draw, w = 1100, h = 750, res = 130) {
  png(file.path(OUT, paste0(base, ".png")),
      width = round(w*FIG_SCALE), height = round(h*FIG_SCALE), res = round(res*FIG_SCALE))
  op <- par(no.readonly = TRUE); tryCatch(draw(), finally = { par(op); dev.off() })
  # cairo_pdf, NOT pdf(): the legacy device uses Type 1 fonts with no glyph for
  # the circles/arrows/Greek used here, so those characters silently vanish.
  if (capabilities("cairo")) cairo_pdf(file.path(OUT, paste0(base, ".pdf")),
                                       width = w/res, height = h/res)
  else pdf(file.path(OUT, paste0(base, ".pdf")), width = w/res, height = h/res)
  op <- par(no.readonly = TRUE); tryCatch(draw(), finally = { par(op); dev.off() })
}

# ---- per-cell helpers -------------------------------------------------------
thin_by_unit <- function(df, cap) do.call(rbind, lapply(split(df, df$ID), function(u)
  if (nrow(u) <= cap) u else u[unique(round(seq(1, nrow(u), length.out = cap))), ]))

# The x-range genuinely SUPPORTED by the data: the widest interval spanned by at
# least MIN_UNITS_X shoulders. Absolute min/max is not usable — e.g. sternoclavicular
# horizontal flexion has one unit (17 of 13,064 rows) reaching 100 deg while the
# other nine stop near 0, which dragged the fitted curve across empty space.
support_range <- function(x, id, min_units = MIN_UNITS_X) {
  rr <- do.call(rbind, lapply(split(x, droplevels(id)), range))
  if (is.null(rr) || nrow(rr) < min_units) return(range(x))
  g <- seq(min(rr[,1]), max(rr[,2]), length.out = 800)
  n <- sapply(g, function(v) sum(rr[,1] <= v & rr[,2] >= v))
  ok <- which(n >= min_units)
  if (!length(ok)) return(range(x))
  c(g[min(ok)], g[max(ok)])
}
# x-range where BOTH conditions are supported
overlap_range <- function(df) {
  cs <- levels(droplevels(df$cond))
  rr <- lapply(cs, function(cc) { k <- df$cond == cc
                                  support_range(df$TIME[k], df$ID[k]) })
  if (length(rr) < 2) return(rr[[1]])
  c(max(sapply(rr, `[`, 1)), min(sapply(rr, `[`, 2)))
}
# interior knots from the cell's own quantiles; never duplicated, never on the edge
knots_for <- function(x, kdf) {
  q <- as.numeric(quantile(x, seq(0, 1, length.out = kdf + 1)[2:kdf]))
  q <- unique(round(q, 1)); q <- q[q > min(x) & q < max(x)]
  if (length(q) < 1) NULL else q
}
wald_from <- function(tb) {
  tm <- rownames(tb); w <- tb[, "coef"] / tb[, "Se"]
  data.frame(curve = ifelse(grepl("cond", tm), "difference", "reference"),
             term = tm, estimate = tb[, "coef"], se = tb[, "Se"], wald = w,
             p_value = tb[, "p-value"], passed = abs(w) >= WALD_CRIT, row.names = NULL)
}
# joint Wald on a block of coefficients: g' V^-1 g ~ chi2_df
joint_wald <- function(g, V) {
  Vi <- try(solve(V), silent = TRUE)
  if (inherits(Vi, "try-error")) return(c(chi2 = NA, df = length(g), p = NA))
  w <- as.numeric(t(g) %*% Vi %*% g)
  c(chi2 = w, df = length(g), p = pchisq(w, length(g), lower.tail = FALSE))
}
# symbol strip: one glyph per coefficient, filled = passed the Wald test
sym_strip <- function(passed) paste(ifelse(passed, "●", "○"), collapse = "")

if (SUMMARY_ONLY) {
  master <- read.csv(file.path(OUT, "00_master_summary.csv"), stringsAsFactors = FALSE)
  cat("summary-only: reloaded", nrow(master), "rows from 00_master_summary.csv\n")
} else {
  if (!file.exists(LONG)) stop("Missing ", LONG, " — run prepare_monolix_data.py first.")
  cat("Reading", LONG, "...\n")
  D <- read.csv(LONG, stringsAsFactors = FALSE)
  D$cond <- factor(ifelse(D$in_vivo %in% c("True","TRUE","1","yes"), "in vivo", "ex vivo"),
                   levels = c("ex vivo","in vivo"))
  D$study <- sub("_[^_]*$", "", D$ID)
  D <- D[is.finite(D$TIME) & is.finite(D$Y), ]

  combos <- expand.grid(joint = JOINT_ORDER, motion = MOTION_ORDER, dof = 1:3,
                        stringsAsFactors = FALSE)
  cat(sprintf("%d combinations to sweep (K = %d, D diagonal, trimmed)\n", nrow(combos), KDF))

  res <- list(); wald_rows <- list(); CUR <- list()
  t_start <- Sys.time()

  for (i in seq_len(nrow(combos))) {
    jt <- combos$joint[i]; mo <- combos$motion[i]; df_ <- combos$dof[i]
    key <- paste(jt, mo, df_, sep = "|")
    sub <- D[D$joint == jt & D$humeral_motion == mo & D$degree_of_freedom == df_, ]
    row <- data.frame(joint = jt, motion = mo, dof = df_, n_rows = nrow(sub),
                      n_units = 0, n_ex = 0, n_in = 0, st_ex = 0, st_in = 0,
                      mode = "skip", x_lo = NA, x_hi = NA, knots = NA, n_fitted = 0,
                      converged = NA, BIC = NA,
                      level_shift = NA, level_wald = NA,
                      joint_chi2 = NA, joint_df = NA, joint_p = NA,
                      shape_chi2 = NA, shape_df = NA, shape_p = NA,
                      wald_pass_ref = NA, wald_total_ref = NA,
                      wald_pass_diff = NA, wald_total_diff = NA, sym_diff = NA,
                      diff_mean_signed = NA, diff_mean_abs = NA, diff_max_abs = NA,
                      diff_lo_curve = NA, diff_hi_curve = NA, sig_frac_x = NA,
                      diff_mean_lo = NA, diff_mean_hi = NA,
                      resid_sd = NA, kurtosis = NA, seconds = NA, note = "",
                      stringsAsFactors = FALSE)

    if (nrow(sub) < 20) { row$note <- "too few rows"; res[[length(res)+1]] <- row; next }
    sub$ID <- factor(sub$ID); sub <- sub[order(sub$ID, sub$TIME), ]
    units_by <- tapply(as.character(sub$cond), sub$ID, function(z) z[1])
    row$n_units <- length(units_by)
    row$n_ex <- sum(units_by == "ex vivo"); row$n_in <- sum(units_by == "in vivo")
    row$st_ex <- length(unique(sub$study[sub$cond == "ex vivo"]))
    row$st_in <- length(unique(sub$study[sub$cond == "in vivo"]))
    compare <- row$n_ex >= MIN_PER_COND && row$n_in >= MIN_PER_COND
    if (!compare && row$n_units < MIN_TOTAL) {
      row$note <- sprintf("too few units (%d ex / %d in)", row$n_ex, row$n_in)
      res[[length(res)+1]] <- row; next
    }

    st <- thin_by_unit(sub, CAP); st <- st[order(st$ID, st$TIME), ]
    full_rng <- range(st$TIME)                      # boundary knots from the UNTRIMMED range
    if (compare) { ov <- overlap_range(st); st <- st[st$TIME >= ov[1] & st$TIME <= ov[2], ] }
    else ov <- support_range(st$TIME, st$ID)
    if (!compare) st <- st[st$TIME >= ov[1] & st$TIME <= ov[2], ]
    st$ID <- droplevels(st$ID); st$IDnum <- as.integer(st$ID)
    if (nrow(st) < 20 || nlevels(st$ID) < MIN_TOTAL) {
      row$note <- "empty after trim"; res[[length(res)+1]] <- row; next }

    kn <- knots_for(st$TIME, KDF)
    Bk <- try(ns(st$TIME, knots = kn, Boundary.knots = full_rng), silent = TRUE)
    if (inherits(Bk, "try-error")) { row$note <- "basis failed"; res[[length(res)+1]] <- row; next }
    Kc <- ncol(Bk); nsv <- paste0("ns", 1:Kc)
    for (j in 1:Kc) st[[nsv[j]]] <- Bk[, j]
    row$x_lo <- ov[1]; row$x_hi <- ov[2]; row$knots <- paste(kn, collapse = "/")
    row$n_fitted <- nrow(st); row$mode <- if (compare) "compare" else "single"

    ff <- if (compare)
      as.formula(paste("Y ~", paste(nsv, collapse=" + "), "+ cond +",
                       paste(paste0("cond:", nsv), collapse=" + ")))
    else as.formula(paste("Y ~", paste(nsv, collapse=" + ")))
    rf <- as.formula(paste("~", paste(nsv, collapse = " + ")))

    t0 <- Sys.time()
    m <- try(hlme(fixed = ff, random = rf, subject = "IDnum", ng = 1, idiag = TRUE,
                  data = st, verbose = FALSE, nproc = NPROC), silent = TRUE)
    row$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    okfit <- !inherits(m, "try-error") && !is.null(m$conv) && m$conv == 1 && all(is.finite(m$best))
    row$converged <- isTRUE(okfit)
    if (!okfit) { row$note <- "did not converge"; res[[length(res)+1]] <- row; next }

    invisible(capture.output(tb <- summary(m)))
    nf <- if (compare) 2*Kc + 2 else Kc + 1
    V  <- VarCov(m); bf <- m$best[1:nf]; Vf <- V[1:nf, 1:nf]
    rr <- m$pred$resid_ss; zz <- (rr - mean(rr))/sd(rr)
    row$BIC <- m$BIC; row$resid_sd <- sd(rr); row$kurtosis <- mean(zz^4) - 3

    wr <- wald_from(tb); wr$curve[!grepl("cond", wr$term)] <- "reference"
    isref <- wr$curve == "reference"
    row$wald_pass_ref  <- sum(wr$passed[isref]);  row$wald_total_ref  <- sum(isref)
    wald_rows[[length(wald_rows)+1]] <- cbind(joint = jt, motion = mo, dof = df_, wr)

    Xrow <- function(x, inv) { b <- predict(Bk, x)
      if (compare) cbind(1, b, inv, b*inv) else cbind(1, b) }
    band <- function(X) { f <- as.vector(X %*% bf)
      s <- sqrt(pmax(0, rowSums((X %*% Vf) * X)))
      data.frame(fit = f, lo = f - 1.96*s, hi = f + 1.96*s) }
    curves <- do.call(rbind, lapply(levels(droplevels(st$cond)), function(cc) {
      k <- st$cond == cc
      r2 <- support_range(st$TIME[k], st$ID[k])          # never draw past this condition's support
      r2 <- c(max(r2[1], ov[1]), min(r2[2], ov[2]))
      if (!(r2[2] > r2[1])) return(NULL)
      xc <- seq(r2[1], r2[2], length.out = 250)
      cbind(cond = cc, x = xc, band(Xrow(xc, if (cc == "in vivo") 1 else 0))) }))
    dif <- NULL
    if (compare) {
      gi <- (Kc+2):nf                                  # gamma_0 .. gamma_K
      jw <- joint_wald(bf[gi], Vf[gi, gi, drop = FALSE])
      sw <- joint_wald(bf[gi[-1]], Vf[gi[-1], gi[-1], drop = FALSE])
      row$joint_chi2 <- jw["chi2"]; row$joint_df <- jw["df"]; row$joint_p <- jw["p"]
      row$shape_chi2 <- sw["chi2"]; row$shape_df <- sw["df"]; row$shape_p <- sw["p"]
      row$level_shift <- bf[Kc+2]; row$level_wald <- bf[Kc+2]/sqrt(Vf[Kc+2, Kc+2])
      row$wald_pass_diff <- sum(wr$passed[!isref]); row$wald_total_diff <- sum(!isref)
      row$sym_diff <- sym_strip(wr$passed[!isref])
      xd <- seq(ov[1], ov[2], length.out = 250)
      bd <- predict(Bk, xd)
      Xd  <- cbind(0, matrix(0, length(xd), Kc), 1, bd)
      dif <- cbind(x = xd, band(Xd))
      row$diff_mean_signed <- mean(dif$fit); row$diff_mean_abs <- mean(abs(dif$fit))
      row$diff_max_abs <- max(abs(dif$fit))
      row$diff_lo_curve <- min(dif$fit); row$diff_hi_curve <- max(dif$fit)
      # fraction of the fitted range where the 95% band excludes zero — the
      # "% of range" column of 00_significance_map, exported so it is data
      # rather than a number read off a picture
      row$sig_frac_x <- mean((dif$lo > 0) | (dif$hi < 0))
      # 95% CI of the MEAN difference: the average design row is itself a contrast,
      # so its variance is xbar' V xbar — exact, not an average of pointwise SEs
      xbar <- colMeans(Xd)
      mse  <- sqrt(max(0, as.numeric(t(xbar) %*% Vf %*% xbar)))
      row$diff_mean_lo <- row$diff_mean_signed - 1.96*mse
      row$diff_mean_hi <- row$diff_mean_signed + 1.96*mse
    }
    CUR[[key]] <- list(curves = curves, dif = dif, knots = kn, raw = st)
    res[[length(res)+1]] <- row
    cat(sprintf("[%2d/%d] %-18s %-44s dof%d  %-7s %5.1fs  %s\n", i, nrow(combos),
                jt, substr(mo, 1, 44), df_, row$mode, row$seconds,
                if (compare) sprintf("joint p=%.2g  %s", row$joint_p, row$sym_diff) else ""))
    flush.console()
  }

  master <- do.call(rbind, res); rownames(master) <- NULL
  cmp <- master$mode == "compare" & !is.na(master$joint_p)
  master$joint_p_fdr <- NA; master$joint_p_fdr[cmp] <- p.adjust(master$joint_p[cmp], "BH")
  master$shape_p_fdr <- NA; master$shape_p_fdr[cmp] <- p.adjust(master$shape_p[cmp], "BH")
  write.csv(master, file.path(OUT, "00_master_summary.csv"), row.names = FALSE)
  if (length(wald_rows))
    write.csv(do.call(rbind, wald_rows), file.path(OUT, "00_wald_tests.csv"), row.names = FALSE)

  cat(sprintf("\nswept %d combos in %.1f min | compare %d | single %d | skip %d | failed %d\n",
              nrow(master), as.numeric(difftime(Sys.time(), t_start, units="mins")),
              sum(master$mode=="compare"), sum(master$mode=="single"),
              sum(master$mode=="skip"), sum(!is.na(master$converged) & !master$converged)))

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
    key_of <- function(jt, k) paste(jt, mo, k, sep = "|")

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

  # ---- 00_forest.png — every compared cell on one axis -------------------------
  # One row per cell. Two horizontal spans, deliberately different weights:
  #   wide pale band = the RANGE the difference curve covers across elevation, i.e.
  #                    how much the difference itself varies with arm position
  #   narrow dark bar = the 95% CI of the MEAN difference — the inferential span
  #   dot            = the mean signed difference
  # A cell whose pale band is wide but whose dark bar straddles zero has a
  # difference that changes sign or size with elevation and averages to nothing.
  fo <- master[master$mode == "compare" & !is.na(master$diff_mean_signed), ]
  if (nrow(fo)) {
    fo$mkey <- match(fo$motion, MOTION_ORDER); fo$jkey <- match(fo$joint, JOINT_ORDER)
    fo <- fo[order(fo$mkey, fo$jkey, fo$dof), ]
    fo <- fo[rev(seq_len(nrow(fo))), ]                     # first motion at the top
    lab <- sprintf("%s · DoF%d — %s", capj(fo$joint), fo$dof,
                   mapply(function(j,k) DOF_LEGEND[[j]][k], fo$joint, fo$dof))
    hgt <- max(600, 150 + 30*nrow(fo))
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

}

# ---- 00_significance_map.png — WHERE each cell differs -----------------------
# The forest answers "by how much"; this answers "over which elevations".
# One row per compared cell, x = elevation:
#   pale bar  = the range actually fitted (the x-overlap of the two conditions)
#   solid bar = the sub-ranges where the 95% band on the difference excludes zero
# A cell can have a large mean difference and almost no solid bar, which means the
# difference is real on average but nowhere separable from zero point by point.
if (nrow(fo)) {
  segs <- lapply(seq_len(nrow(fo)), function(i) {
    e <- CUR[[paste(fo$joint[i], fo$motion[i], fo$dof[i], sep = "|")]]
    if (is.null(e) || is.null(e$dif)) return(NULL)
    z <- e$dif; sig <- (z$lo > 0) | (z$hi < 0)
    r <- rle(sig); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
    k <- which(r$values)
    if (!length(k)) return(data.frame(lo = numeric(0), hi = numeric(0)))
    data.frame(lo = z$x[starts[k]], hi = z$x[ends[k]])
  })
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
