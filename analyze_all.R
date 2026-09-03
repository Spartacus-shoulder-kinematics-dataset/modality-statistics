# =============================================================================
# GENERALIZED SWEEP — spline mixed model for every
#   joint (shoulder joint) x humeral_motion (movement) x degree_of_freedom (angle)
# across all angular data (unit == rad), then assemble PLANCHES (one per motion).
#
# Per-combo model (validated in analysis_shoulder.R):
#   compare (both conditions): Y ~ condO + s(TIME) + s(TIME,by=condO) + s(ID,re) [+AR1]
#   single  (one condition)  : Y ~ s(TIME) + s(ID, bs="re")                      [+AR1]
#
# Outputs (figures/generalized/):
#   planche_<motion>.png                one plate per movement:
#                                       rows = joints (sternoclav -> acromioclav
#                                       -> scapulothoracic -> glenohumeral),
#                                       cols = DoF 1/2/3, each panel = in-vivo vs
#                                       ex-vivo fitted curves over raw data.
#   <joint>/<motion-slug>/dof<k>/*.png  per-combo drill-down figures
#   00_master_summary.csv               one row per combo (effect sizes, p-values)
#   00_SUMMARY.md                       human-readable index
#
# Run:  python3 prepare_monolix_data.py   # writes spartacus_angles_long.csv
#       Rscript analyze_all.R
# =============================================================================

suppressMessages({ library(mgcv); library(dplyr) })

LONG    <- "spartacus_angles_long.csv"
OUT     <- file.path("figures", "generalized")
MIN_PER_COND <- 3     # >=3 units in BOTH conditions -> compare mode
MIN_TOTAL    <- 4     # >=4 units total              -> at least single mode
COL <- c("ex vivo" = "#d95f02", "in vivo" = "#1b9e77")

# joint row order (proximal -> distal) and anatomical DoF names (col headers),
# matching shoulder-kinematics/spartacus/plots/constants_plot.py
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
slug <- function(s) gsub("(^_|_$)", "", gsub("[^a-z0-9]+", "_", tolower(s)))
capj <- function(j) paste0(toupper(substring(j, 1, 1)), substring(j, 2))
stars <- function(p) if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"

# annotate one planche panel from its master-summary row `ri`
annot <- function(ri) {
  u <- par("usr"); dx <- u[2]-u[1]; dy <- u[4]-u[3]
  if (nrow(ri) != 1) return(invisible())
  if (ri$mode == "compare") {
    sig <- if ("diff_p_fdr" %in% names(ri)) ri$diff_p_fdr else ri$diff_p   # FDR-adjusted
    text(u[2]-0.03*dx, u[4]-0.04*dy, stars(sig), adj = c(1, 1), font = 2, cex = 1.2, col = "grey10")
    ln <- c(sprintf("%d/%d sh · %d/%d st", ri$n_ex, ri$n_in, ri$st_ex, ri$st_in),
            sprintf("Δ %+.0f° (|max| %.0f°)", ri$diff_mean_signed, ri$diff_max_abs),
            sprintf("shift %+.0f° %s", ri$level_shift, stars(ri$level_p)))
    for (k in seq_along(ln))
      text(u[1]+0.03*dx, u[3]+0.05*dy + (length(ln)-k)*0.08*dy, ln[k], adj = c(0, 0), cex = 0.62, col = "grey30")
  } else if (ri$mode == "single") {
    text(u[1]+0.03*dx, u[3]+0.05*dy, paste0("single group (", sub("^only: ", "", ri$note), ")"),
         adj = c(0, 0), cex = 0.62, col = "grey30")
  }
}

if (!file.exists(LONG))
  stop("Missing ", LONG, " — run  python3 prepare_monolix_data.py  first.")
cat("Reading", LONG, "...\n")
D <- read.csv(LONG, stringsAsFactors = FALSE)
D$cond <- factor(ifelse(D$in_vivo %in% c("True","TRUE","1","yes"),
                        "in vivo", "ex vivo"), levels = c("ex vivo","in vivo"))
if (is.null(D$study)) D$study <- sub("_[^_]*$", "", D$ID)   # recover article if absent

png2 <- function(path, expr, w = 1000, h = 680, res = 120) {
  png(path, width = w, height = h, res = res); on.exit(dev.off()); force(expr)
}

# ---- fit one combination; return list(row = summary, curves = fitted curves) --
analyze_one <- function(sub, joint, motion, dof) {
  sub <- sub[is.finite(sub$TIME) & is.finite(sub$Y), ]
  sub$ID <- factor(sub$ID)
  sub <- sub[order(sub$ID, sub$TIME), ]
  n_ex <- length(unique(sub$ID[sub$cond == "ex vivo"]))
  n_in <- length(unique(sub$ID[sub$cond == "in vivo"]))
  # STUDY-level replication (the real n for the condition contrast: it is
  # confounded with study — see STATISTICAL_METHODS_REVIEW.md)
  st_ex <- length(unique(sub$study[sub$cond == "ex vivo"]))
  st_in <- length(unique(sub$study[sub$cond == "in vivo"]))

  row <- data.frame(joint = joint, motion = motion, dof = dof,
                    n_rows = nrow(sub), n_units = nlevels(sub$ID),
                    n_ex = n_ex, n_in = n_in, st_ex = st_ex, st_in = st_in,
                    mode = NA_character_, rho = NA_real_,
                    level_shift = NA_real_, level_p = NA_real_, diff_p = NA_real_,
                    diff_mean_signed = NA_real_, diff_mean_abs = NA_real_,
                    diff_max_abs = NA_real_, bic = NA_real_,
                    note = "", stringsAsFactors = FALSE)
  empty <- function(r) list(row = r, curves = NULL)

  if (nlevels(sub$ID) < MIN_TOTAL) { row$mode <- "skip"; row$note <- "too few units"
    return(empty(row)) }

  outdir <- file.path(OUT, slug(joint), slug(motion), paste0("dof", dof))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  sub$ar_start <- c(TRUE, sub$ID[-1] != sub$ID[-nrow(sub)])
  compare <- (n_ex >= MIN_PER_COND && n_in >= MIN_PER_COND)

  fit_ar <- function(form) {
    m0  <- bam(form, data = sub, method = "fREML", discrete = TRUE)
    r   <- resid(m0)                        # within-shoulder lag-1 autocorrelation
    num <- den <- 0
    for (u in levels(sub$ID)) { ru <- r[sub$ID == u]
      if (length(ru) > 1) { num <- num + sum(head(ru,-1)*tail(ru,-1)); den <- den + sum(ru^2) } }
    rho <- num / den
    list(m = bam(form, data = sub, method = "fREML", discrete = TRUE,
                 rho = rho, AR.start = sub$ar_start), rho = rho)
  }
  # predict one condition's curve clipped to its own x-support
  cond_curve <- function(m, cc, condO_levels) {
    rr <- range(sub$TIME[sub$cond == cc]); xc <- seq(rr[1], rr[2], length.out = 200)
    nd <- data.frame(TIME = xc, ID = sub$ID[1],
                     condO = factor(cc, levels = condO_levels, ordered = TRUE))
    pc <- predict(m, nd, exclude = "s(ID)", newdata.guaranteed = TRUE, se.fit = TRUE)
    data.frame(cond = cc, x = xc, fit = pc$fit, lo = pc$fit-1.96*pc$se.fit, hi = pc$fit+1.96*pc$se.fit)
  }

  if (compare) {
    sub$condO <- as.ordered(sub$cond)
    fa <- fit_ar(Y ~ condO + s(TIME) + s(TIME, by = condO) + s(ID, bs = "re"))
    m  <- fa$m; row$mode <- "compare"; row$rho <- fa$rho; row$bic <- BIC(m)
    st <- summary(m)$s.table; pt <- summary(m)$p.table
    dr <- grep("condO", rownames(st))[1]; lr <- grep("condO", rownames(pt))[1]
    row$diff_p      <- if (!is.na(dr)) st[dr, "p-value"]  else NA_real_
    row$level_shift <- if (!is.na(lr)) pt[lr, "Estimate"] else NA_real_
    row$level_p     <- if (!is.na(lr)) pt[lr, "Pr(>|t|)"] else NA_real_

    curves <- do.call(rbind, lapply(levels(sub$condO), cond_curve,
                                    m = m, condO_levels = levels(sub$condO)))
    png2(file.path(outdir, "01_population_curves.png"), {
      plot(sub$TIME, sub$Y, col = adjustcolor(COL[as.character(sub$cond)], 0.15),
           pch = 16, cex = 0.4, xlab = "humerothoracic elevation (x)", ylab = "joint angle (y)",
           main = sprintf("%s | %s | DoF %d", joint, motion, dof))
      for (cc in unique(curves$cond)) { g <- curves[curves$cond == cc, ]
        polygon(c(g$x, rev(g$x)), c(g$lo, rev(g$hi)), col = adjustcolor(COL[cc], 0.25), border = NA)
        lines(g$x, g$fit, col = COL[cc], lwd = 3) }
      legend("topleft", names(COL), col = COL, lwd = 3, bty = "n", cex = 0.9)
    })
    # difference on the OVERLAP x-range only (avoids nonsense extrapolation)
    rng_ex <- range(sub$TIME[sub$cond=="ex vivo"]); rng_in <- range(sub$TIME[sub$cond=="in vivo"])
    ov_lo <- max(rng_ex[1], rng_in[1]); ov_hi <- min(rng_ex[2], rng_in[2])
    if (ov_hi > ov_lo) {
      xo <- seq(ov_lo, ov_hi, length.out = 250)
      mk <- function(cc) data.frame(TIME = xo, ID = sub$ID[1],
                     condO = factor(cc, levels = levels(sub$condO), ordered = TRUE))
      Xd  <- predict(m, mk("in vivo"), type = "lpmatrix", newdata.guaranteed = TRUE) -
             predict(m, mk("ex vivo"), type = "lpmatrix", newdata.guaranteed = TRUE)
      dif <- as.vector(Xd %*% coef(m)); se <- sqrt(rowSums((Xd %*% vcov(m)) * Xd))
      row$diff_mean_signed <- mean(dif)                # signed (shows crossings)
      row$diff_mean_abs <- mean(abs(dif)); row$diff_max_abs <- max(abs(dif))
      png2(file.path(outdir, "02_difference.png"), {
        plot(xo, dif, type = "n", ylim = range(c(dif-2*se, dif+2*se, 0)),
             xlab = "humerothoracic elevation (x)  [overlap]", ylab = "in vivo - ex vivo",
             main = sprintf("difference: mean|.| = %.1f deg (p = %.2g)", row$diff_mean_abs, row$diff_p))
        polygon(c(xo, rev(xo)), c(dif-1.96*se, rev(dif+1.96*se)), col = adjustcolor("#377eb8", 0.25), border = NA)
        lines(xo, dif, lwd = 3, col = "#377eb8"); abline(h = 0, lty = 2, col = "grey40")
      })
    } else row$note <- "no x-overlap between conditions"
    return(list(row = row, curves = curves))
  }

  # single mode: one pooled trajectory (condition not tested)
  fa <- fit_ar(Y ~ s(TIME) + s(ID, bs = "re")); m <- fa$m
  row$mode <- "single"; row$rho <- fa$rho; row$bic <- BIC(m)
  row$note <- paste0("only: ", paste(sort(unique(as.character(sub$cond))), collapse = "+"))
  xc <- seq(min(sub$TIME), max(sub$TIME), length.out = 200)
  pp <- predict(m, data.frame(TIME = xc, ID = sub$ID[1]),
                exclude = "s(ID)", newdata.guaranteed = TRUE, se.fit = TRUE)
  curves <- data.frame(cond = "(pooled)", x = xc, fit = pp$fit,
                       lo = pp$fit-1.96*pp$se.fit, hi = pp$fit+1.96*pp$se.fit)
  png2(file.path(outdir, "01_population_curve.png"), {
    plot(sub$TIME, sub$Y, col = adjustcolor(COL[as.character(sub$cond)], 0.15), pch = 16, cex = 0.4,
         xlab = "humerothoracic elevation (x)", ylab = "joint angle (y)",
         main = sprintf("%s | %s | DoF %d  (%s)", joint, motion, dof, row$note))
    polygon(c(xc, rev(xc)), c(curves$lo, rev(curves$hi)), col = adjustcolor("grey50", 0.3), border = NA)
    lines(xc, curves$fit, lwd = 3)
  })
  list(row = row, curves = curves)
}

# ---- loop over every combination --------------------------------------------
combos <- D %>% distinct(joint, humeral_motion, degree_of_freedom) %>%
  arrange(joint, humeral_motion, degree_of_freedom)
cat(sprintf("Sweeping %d combinations...\n", nrow(combos)))

rows <- vector("list", nrow(combos)); CURVES <- list()
for (i in seq_len(nrow(combos))) {
  j <- combos$joint[i]; mo <- combos$humeral_motion[i]; dof <- combos$degree_of_freedom[i]
  sub <- D[D$joint == j & D$humeral_motion == mo & D$degree_of_freedom == dof, ]
  out <- tryCatch(analyze_one(sub, j, mo, dof), error = function(e) list(
    row = data.frame(joint=j, motion=mo, dof=dof, n_rows=nrow(sub), n_units=NA, n_ex=NA,
      n_in=NA, st_ex=NA, st_in=NA, mode="error", rho=NA, level_shift=NA, level_p=NA,
      diff_p=NA, diff_mean_signed=NA, diff_mean_abs=NA, diff_max_abs=NA, bic=NA,
      note=conditionMessage(e), stringsAsFactors=FALSE), curves = NULL))
  rows[[i]] <- out$row
  if (!is.null(out$curves)) CURVES[[paste(j, mo, dof)]] <- out$curves
  cat(sprintf("  [%2d/%2d] %-16s | %-42s | dof%s -> %s\n", i, nrow(combos), j, mo, dof, out$row$mode))
}
res <- do.call(rbind, rows)
# Benjamini-Hochberg FDR across the family of compared cells (multiplicity;
# 33 tests). p-values remain EXPLORATORY: n is study-level (st_ex/st_in), not rows.
cmp_idx <- which(res$mode == "compare")
res$diff_p_fdr <- NA_real_; res$level_p_fdr <- NA_real_
res$diff_p_fdr[cmp_idx]  <- p.adjust(res$diff_p[cmp_idx],  method = "BH")
res$level_p_fdr[cmp_idx] <- p.adjust(res$level_p[cmp_idx], method = "BH")
write.csv(res, file.path(OUT, "00_master_summary.csv"), row.names = FALSE)

# ---- PLANCHES: one plate per motion (rows = joints, cols = DoF 1/2/3) --------
motions <- intersect(MOTION_ORDER, unique(res$motion))
motions <- c(motions, setdiff(unique(res$motion), motions))
for (mo in motions) {
  png(file.path(OUT, paste0("planche_", slug(mo), ".png")), width = 1300, height = 1620, res = 135)
  par(mfrow = c(4, 3), oma = c(5, 5, 9, 2), mar = c(2.6, 3.4, 3.2, 1))
  xr_mo <- range(D$TIME[D$humeral_motion == mo & is.finite(D$TIME)])
  for (jt in JOINT_ORDER) {
    # common y-range for this joint ROW: raw data + fitted CI across the 3 DoF
    yy <- numeric(0)
    for (dof in 1:3) {
      yy <- c(yy, D$Y[D$joint == jt & D$humeral_motion == mo & D$degree_of_freedom == dof])
      cv <- CURVES[[paste(jt, mo, dof)]]; if (!is.null(cv)) yy <- c(yy, cv$lo, cv$hi)
    }
    yy <- yy[is.finite(yy)]; rr <- if (length(yy)) range(yy) else c(0, 1)
    for (dof in 1:3) {
      raw <- D[D$joint == jt & D$humeral_motion == mo & D$degree_of_freedom == dof, ]
      ttl <- DOF_LEGEND[[jt]][dof]
      if (nrow(raw) == 0) {
        plot(NA, xlim = xr_mo, ylim = rr, xlab = "", ylab = "", main = ttl,   # same row y-axis
             cex.main = 0.85, font.main = 1, cex.axis = 0.8); box(col = "grey85")
        text(mean(xr_mo), mean(rr), "no data", col = "grey60")
      } else {
        plot(raw$TIME, raw$Y, col = adjustcolor(COL[as.character(raw$cond)], 0.18), pch = 16, cex = 0.3,
             xlab = "", ylab = "", main = ttl, cex.main = 0.85, font.main = 1,
             ylim = rr, cex.axis = 0.8)                                        # shared row y-axis
        cv <- CURVES[[paste(jt, mo, dof)]]
        if (!is.null(cv)) for (cc in unique(cv$cond)) { g <- cv[cv$cond == cc, ]
          colc <- if (cc %in% names(COL)) COL[cc] else "grey20"
          polygon(c(g$x, rev(g$x)), c(g$lo, rev(g$hi)), col = adjustcolor(colc, 0.22), border = NA)
          lines(g$x, g$fit, col = colc, lwd = 2.6) }
        annot(res[res$joint == jt & res$motion == mo & res$dof == dof, ])   # stars + stats
      }
      if (dof == 1) mtext(paste0(capj(jt), " (deg)"), side = 2, line = 2.4, cex = 0.72, font = 2)
    }
  }
  mtext(paste0("Shoulder kinematics  —  ", mo), outer = TRUE, side = 3, line = 6.4, font = 2, cex = 1.2)
  mtext("in vivo (green) · ex vivo (orange) · black = single-condition pooled fit    [spline fit +/-95% CI over raw data]",
        outer = TRUE, side = 3, line = 4.6, cex = 0.78)
  mtext("top-right stars = trajectory-difference test, BH-FDR adjusted (EXPLORATORY):  *** p<.001  ** p<.01  * p<.05  ns  -- condition is confounded with study",
        outer = TRUE, side = 3, line = 2.9, cex = 0.66, col = "grey25")
  mtext("panel text: shoulders ex/in · STUDIES ex/in (the real replication) · signed Δ & |max| (deg) · level shift (deg) + its significance",
        outer = TRUE, side = 3, line = 1.4, cex = 0.66, col = "grey45")
  mtext("rows = joints (proximal -> distal)   |   cols = degree of freedom 1 / 2 / 3",
        outer = TRUE, side = 3, line = 0.1, cex = 0.68, col = "grey45")
  mtext("Humerothoracic angle (deg)", outer = TRUE, side = 1, line = 2.6, cex = 0.9)
  dev.off()
  message("  planche_", slug(mo), ".png")
}

# ---- human-readable summary --------------------------------------------------
n_cmp <- sum(res$mode == "compare")
n_sig <- sum(res$mode == "compare" & is.finite(res$diff_p) & res$diff_p < 0.05, na.rm = TRUE)
md <- c("# Generalized sweep — results index", "",
  sprintf("Ran over **%d** combinations (joint x movement x DoF) of angular data.", nrow(res)),
  sprintf("- **compare** (in-vivo vs ex-vivo tested): %d", n_cmp),
  sprintf("- **single** (one condition only): %d", sum(res$mode == "single")),
  sprintf("- **skipped** (too few units): %d", sum(res$mode == "skip")),
  sprintf("- **errors**: %d", sum(res$mode == "error")),
  sprintf("- of the compared, **%d** reach BH-FDR p < 0.05 on the shape test (EXPLORATORY).",
          sum(res$mode=="compare" & is.finite(res$diff_p_fdr) & res$diff_p_fdr < 0.05, na.rm=TRUE)),
  "",
  "> **These contrasts are EXPLORATORY and NOT causal.** Condition is perfectly",
  "> confounded with source study (no study measured both), so a difference mixes",
  "> muscle activation with protocol / measurement / population. The real",
  "> replication is the number of **studies** (`st ex/in`), not rows or shoulders,",
  "> and p-values are pseudo-replicated — read the **effect size in degrees** and",
  "> the study counts, not the stars. p-values are BH-FDR adjusted across the",
  "> compared family. See `STATISTICAL_METHODS_REVIEW.md` / `_RESPONSE.md`.",
  "", "**Planches** (`planche_<motion>.png`): one plate per movement, rows = joints",
  "(sternoclavicular -> acromioclavicular -> scapulothoracic -> glenohumeral),",
  "cols = DoF 1/2/3. See `00_master_summary.csv` for every combo and each combo's",
  "folder for its drill-down curves.",
  "", "## Compared combinations (sorted by |effect size|, largest first)", "",
  "| joint | movement | DoF | shoulders ex/in | **studies ex/in** | signed Δ° | \\|max\\|° | level shift° | level p (FDR) | shape p (FDR) |",
  "|---|---|---|---|---|---|---|---|---|---|")
cmp <- res[res$mode == "compare", ]; cmp <- cmp[order(-cmp$diff_mean_abs), ]
for (r in seq_len(nrow(cmp))) md <- c(md, sprintf(
  "| %s | %s | %s | %d/%d | **%d/%d** | %+.1f | %.1f | %+.2f | %.2g | %.2g |",
  cmp$joint[r], cmp$motion[r], cmp$dof[r], cmp$n_ex[r], cmp$n_in[r], cmp$st_ex[r], cmp$st_in[r],
  cmp$diff_mean_signed[r], cmp$diff_max_abs[r], cmp$level_shift[r], cmp$level_p_fdr[r], cmp$diff_p_fdr[r]))
writeLines(md, file.path(OUT, "00_SUMMARY.md"))

cat(sprintf("\nDONE. %d combos: %d compared (%d significant), %d single, %d skipped, %d errors.\n",
    nrow(res), n_cmp, n_sig, sum(res$mode=="single"), sum(res$mode=="skip"), sum(res$mode=="error")))
cat("Planches:", length(motions), "| see", file.path(OUT, "00_SUMMARY.md"), "\n")
