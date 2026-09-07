# =============================================================================
# ITERATION 04 — RANDOM spline mixed model with lcmm::hlme
#
# Iteration 03 left one problem open: the within-shoulder residual
# autocorrelation was 0.970 and UNMODELLED, because a shoulder's deviation from
# the population curve is itself a smooth function of elevation and a random
# INTERCEPT absorbs only its level.
#
# This iteration applies the canonical lcmm remedy (Prague & Proust-Lima), namely
# putting the whole spline basis in the random effects so each shoulder gets its
# own random CURVE:
#
#     fixed  = Y ~ ns(TIME, knots) * cond
#     random = ~ ns(TIME, knots)              <- the change
#
# Two further changes adopted from the same source:
#   * EXPLICIT knots instead of df =. Iteration 03's data-driven boundary knots
#     landed at -3.46 and 187.02, which is why its gamma_0 was the difference at
#     an extrapolated point and flipped sign with df. Anchoring the lower
#     boundary at 0 deg makes gamma_0 "the difference with the arm at the side".
#   * The random-effect covariance is compared diagonal vs unstructured, since
#     44 shoulders is not many for an unstructured 5x5.
#
# Run:  Rscript analysis_random_spline.R
# Prereq: python3 prepare_monolix_data.py
# =============================================================================

suppressMessages({ library(lcmm); library(splines); library(dplyr) })

INPUT   <- "monolix_st_frontal_dof2.csv"
FIG_DIR <- "figures"
S4      <- "04_natural_spline_hlme"
CAP     <- 100
NPROC   <- max(1L, min(8L, parallel::detectCores() - 1L))
WALD_CRIT <- 1.96
# Explicit and anatomically meaningful, not data-driven quantiles.
KNOTS    <- c(40, 70, 100, 130)     # interior, for the main K = 5 model
BOUNDARY <- c(0, 160)               # 99.6% of rows lie inside; the rest extrapolate linearly
K        <- length(KNOTS) + 1       # ns() columns = interior + 1 = 5
KS       <- 2:5                     # df sweep, as in iteration 03
COL     <- c("ex vivo" = "#d95f02", "in vivo" = "#1b9e77")
BLUE    <- "#377eb8"

OUT <- file.path(FIG_DIR, S4)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fig <- function(name, expr, w = 1100, h = 750, res = 130) {
  png(file.path(OUT, name), width = w, height = h, res = res)
  on.exit(dev.off()); force(expr); message("  saved ", file.path(S4, name))
}

# Knots are PHYSICAL positions on the elevation axis, not an abstract tuning knob,
# so mark them as such rather than as rug ticks:
#   dashed vertical  = interior knot — where the piecewise cubics join
#   dotted vertical  = boundary knot — where the spline is forced linear beyond
#   filled dot       = the fitted value AT that knot, on each curve
# The dots are the informative half: they tie a knot to the trajectory it shapes.
knot_marks <- function(interior, boundary = NULL, curves = NULL, colmap = COL,
                       xvar = "x", yvar = "fit", byvar = "cond", cex = 1.1) {
  usr <- par("usr")
  ik <- interior[interior > usr[1] & interior < usr[2]]
  if (length(ik)) abline(v = ik, lty = 2, lwd = 1, col = adjustcolor("grey35", 0.55))
  if (!is.null(boundary)) {
    bk <- boundary[boundary > usr[1] & boundary < usr[2]]
    if (length(bk)) abline(v = bk, lty = 3, lwd = 1.3, col = adjustcolor("grey15", 0.6))
  }
  if (!is.null(curves) && length(ik)) {
    grps <- if (byvar %in% names(curves)) unique(curves[[byvar]]) else NA
    for (g in grps) {
      cc <- if (is.na(g)) curves else curves[curves[[byvar]] == g, ]
      cc <- cc[order(cc[[xvar]]), ]
      inside <- ik >= min(cc[[xvar]]) & ik <= max(cc[[xvar]])
      if (!any(inside)) next
      yy <- approx(cc[[xvar]], cc[[yvar]], xout = ik[inside])$y
      col <- if (is.na(g)) BLUE else colmap[[as.character(g)]]
      points(ik[inside], yy, pch = 21, bg = col, col = "white", cex = cex, lwd = 1.2)
    }
  }
}

# -----------------------------------------------------------------------------
# 1. LOAD + THIN (identical to iteration 03, so the two are comparable)
# -----------------------------------------------------------------------------
if (!file.exists(INPUT)) stop("Missing ", INPUT, " — run prepare_monolix_data.py first.")
d <- read.csv(INPUT, stringsAsFactors = FALSE)
d$cond <- factor(ifelse(d$in_vivo %in% c("True","TRUE","1","yes"), "in vivo", "ex vivo"),
                 levels = c("ex vivo","in vivo"))
d$ID <- factor(d$ID)
d <- d[is.finite(d$TIME) & is.finite(d$Y), ]
d <- d[order(d$ID, d$TIME), ]

thin_by_unit <- function(df, cap) {
  do.call(rbind, lapply(split(df, df$ID), function(u) {
    if (nrow(u) <= cap) return(u)
    u[unique(round(seq(1, nrow(u), length.out = cap))), ]
  }))
}
dt <- thin_by_unit(d, CAP); dt <- dt[order(dt$ID, dt$TIME), ]
dt$IDnum <- as.integer(dt$ID)

cat("\n================ ITERATION 04: RANDOM spline + hlme ================\n")
cat(sprintf("%d obs -> %d thinned (<=%d/shoulder), %d shoulders\n",
            nrow(d), nrow(dt), CAP, nlevels(dt$ID)))
cat(sprintf("knots: interior %s | boundary %s  (%d of %d rows outside the boundary)\n",
            paste(KNOTS, collapse=","), paste(BOUNDARY, collapse=","),
            sum(dt$TIME < BOUNDARY[1] | dt$TIME > BOUNDARY[2]), nrow(dt)))

# Basis built ONCE with the explicit knots; reused for every prediction.
B <- ns(dt$TIME, knots = KNOTS, Boundary.knots = BOUNDARY)
stopifnot(ncol(B) == K)
nsv <- paste0("ns", 1:K)
for (k in 1:K) dt[[nsv[k]]] <- B[, k]

fixed_f  <- as.formula(paste("Y ~", paste(nsv, collapse=" + "), "+ cond +",
                             paste(paste0("cond:", nsv), collapse=" + ")))
random_f <- as.formula(paste("~", paste(nsv, collapse = " + ")))
nfix <- 2*K + 2

# -----------------------------------------------------------------------------
# 2. THE COMPARISON — how much of the autocorrelation does a random CURVE absorb?
#    lag-1 residual ACF is the target metric: iteration 03 left it at 0.970.
# -----------------------------------------------------------------------------
acf_within1 <- function(r, id) {
  num <- den <- 0
  for (u in unique(id)) { ru <- r[id == u]
    if (length(ru) > 1) { num <- num + sum(head(ru,-1)*tail(ru,-1)); den <- den + sum(ru^2) } }
  num/den
}
acf_curve <- function(r, id, lag.max = 40) {
  acc <- matrix(NA_real_, 0, lag.max + 1)
  for (u in unique(id)) { ru <- r[id == u]
    if (length(ru) > lag.max + 5)
      acc <- rbind(acc, drop(acf(ru, lag.max = lag.max, plot = FALSE)$acf)) }
  colMeans(acc, na.rm = TRUE)
}
norm_stats <- function(r) {
  r <- r[is.finite(r)]; z <- (r - mean(r))/sd(r)
  sw <- { s <- if (length(r) > 5000) sample(r, 5000) else r; shapiro.test(s)$p.value }
  c(skew = mean(z^3), kurtosis = mean(z^4) - 3,
    frac_within_1.96 = mean(abs(z) <= 1.96), shapiro_p = sw)
}
wald_from <- function(tb) {
  tm <- rownames(tb); w <- tb[, "coef"]/tb[, "Se"]
  data.frame(curve = ifelse(grepl("cond", tm), "difference (in vivo - ex vivo)",
                            "reference (ex vivo)"),
             term = tm, estimate = tb[,"coef"], se = tb[,"Se"], wald = w,
             p_value = tb[,"p-value"], passed = abs(w) >= WALD_CRIT, row.names = NULL)
}

# The unstructured variant was dropped: it cost 40x the runtime (659 s vs 16 s)
# for 101 BIC units, an identical lag-1 ACF (0.8352 vs 0.8356) and a slightly
# WORSE kurtosis, and 21 covariance parameters from 44 shoulders is a stretch.
# The off-diagonal covariances of D are therefore not estimated here.
VARIANTS <- list(
  list(id = "R1_intercept",   lab = "random = ~1 (iteration 03)",  rnd = ~1,       idiag = FALSE),
  list(id = "R2_spline_diag", lab = "random = ~ns(...), diagonal", rnd = random_f, idiag = TRUE)
)

fits <- list(); rows <- list()
for (v in VARIANTS) {
  cat(sprintf("\n--- %s ---\n", v$lab)); flush.console()
  t0 <- Sys.time()
  m <- try(hlme(fixed = fixed_f, random = v$rnd, subject = "IDnum", ng = 1,
                idiag = v$idiag, data = dt, verbose = FALSE, nproc = NPROC), silent = TRUE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  okfit <- !inherits(m, "try-error") && !is.null(m$conv) && m$conv == 1 &&
           all(is.finite(m$best))
  if (!okfit) {
    cat(sprintf("  NOT usable (%.0f s)\n", secs))
    rows[[length(rows)+1]] <- data.frame(variant = v$id, spec = v$lab, converged = FALSE,
      n_param = NA, loglik = NA, BIC = NA, resid_acf1 = NA, resid_sd = NA,
      skew = NA, kurtosis = NA, wald_pass = NA, wald_total = NA, seconds = secs)
    next
  }
  invisible(capture.output(tb <- summary(m)))
  # the fixed effects must be the FIRST nfix entries of best — everything below
  # (curve prediction, VarCov slicing) depends on that ordering
  stopifnot(nrow(tb) == nfix,
            identical(names(m$best)[1:nfix], rownames(tb)))
  r  <- m$pred$resid_ss
  ns_st <- norm_stats(r); a1 <- acf_within1(r, dt$IDnum)
  w  <- abs(tb[,"coef"]/tb[,"Se"])
  cat(sprintf("  conv=1 %.0fs | npm=%d | BIC %.0f | lag-1 ACF %.3f | kurt %.2f | Wald %d/%d\n",
              secs, length(m$best), m$BIC, a1, ns_st["kurtosis"], sum(w >= WALD_CRIT), length(w)))
  fits[[v$id]] <- list(m = m, tab = tb, lab = v$lab)
  rows[[length(rows)+1]] <- data.frame(variant = v$id, spec = v$lab, converged = TRUE,
    n_param = length(m$best), loglik = m$loglik, BIC = m$BIC,
    resid_acf1 = a1, resid_sd = sd(r),
    skew = ns_st["skew"], kurtosis = ns_st["kurtosis"],
    wald_pass = sum(w >= WALD_CRIT), wald_total = length(w), seconds = secs)
}
vt <- do.call(rbind, rows); rownames(vt) <- NULL
write.csv(vt, file.path(OUT, "00_variants.csv"), row.names = FALSE)
cat("\n================ VARIANT COMPARISON ================\n"); print(vt, digits = 4)

okv <- vt[vt$converged, ]
if (!nrow(okv)) stop("no variant converged — see 00_variants.csv")
best_id <- okv$variant[which.min(okv$BIC)]
cat(sprintf("\nRETAINED: %s\n", okv$spec[okv$variant == best_id]))
m <- fits[[best_id]]$m; tab <- fits[[best_id]]$tab

# -----------------------------------------------------------------------------
# 3. TABLES
# -----------------------------------------------------------------------------
wald_tab <- wald_from(tab)
write.csv(wald_tab, file.path(OUT, "00_wald_tests.csv"), row.names = FALSE)
wc <- tapply(wald_tab$passed, wald_tab$curve, function(x) sprintf("%d/%d", sum(x), length(x)))
ref_lab  <- wc[["reference (ex vivo)"]]; diff_lab <- wc[["difference (in vivo - ex vivo)"]]
cat("\n================ WALD TESTS ================\n"); print(wald_tab, digits = 4)
cat(sprintf("\nPASSED — reference: %s | difference: %s\n", ref_lab, diff_lab))

V <- VarCov(m)
allp <- data.frame(parameter = names(m$best), estimate = as.vector(m$best),
                   se = sqrt(diag(V)), row.names = NULL)
allp$block <- ifelse(seq_len(nrow(allp)) <= nfix, "fixed_effect", "variance_component")
allp$wald  <- ifelse(allp$block == "fixed_effect", allp$estimate/allp$se, NA)
write.csv(allp, file.path(OUT, "00_coefficients.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. FIGURES
# -----------------------------------------------------------------------------
bf <- m$best[1:nfix]; Vf <- V[1:nfix, 1:nfix]
Xrow  <- function(x, inv) { b <- predict(B, x); cbind(1, b, inv, b*inv) }
Xdiff <- function(x) { b <- predict(B, x); cbind(0, matrix(0, length(x), K), 1, b) }
band  <- function(X) { fit <- as.vector(X %*% bf)
  se <- sqrt(pmax(0, rowSums((X %*% Vf) * X)))
  data.frame(fit = fit, lo = fit - 1.96*se, hi = fit + 1.96*se) }
curves <- do.call(rbind, lapply(levels(dt$cond), function(cc) {
  rr <- range(dt$TIME[dt$cond == cc]); xc <- seq(rr[1], rr[2], length.out = 300)
  cbind(cond = cc, x = xc, band(Xrow(xc, if (cc == "in vivo") 1 else 0)))
}))

fig("01_population_curves_by_condition.png", {
  plot(dt$TIME, dt$Y, col = adjustcolor(COL[as.character(dt$cond)], 0.25), pch = 16, cex = 0.5,
       xlab = "Humerothoracic elevation (x)", ylab = "Scapulothoracic angle (y)",
       main = sprintf("Random spline: %s", fits[[best_id]]$lab))
  for (cc in unique(curves$cond)) { g <- curves[curves$cond == cc, ]
    polygon(c(g$x, rev(g$x)), c(g$lo, rev(g$hi)), col = adjustcolor(COL[cc], 0.3), border = NA)
    lines(g$x, g$fit, col = COL[cc], lwd = 3) }
  knot_marks(KNOTS, BOUNDARY, curves)
  legend("topleft", names(COL), col = COL, lwd = 3, bty = "n")
  mtext(sprintf("knots: %s  (dashed) | boundary %s (dotted); dots = fitted value at each knot",
                paste(KNOTS, collapse=", "), paste(BOUNDARY, collapse="-")),
        side = 3, cex = 0.7, line = 0.2, col = "grey30")
  legend("bottomright", bty = "n", cex = 0.85,
         legend = c(sprintf("Wald |coef/se| >= %.2f", WALD_CRIT),
                    sprintf("reference : %s", ref_lab), sprintf("difference: %s", diff_lab)))
})

ov <- c(max(min(dt$TIME[dt$cond=="ex vivo"]), min(dt$TIME[dt$cond=="in vivo"])),
        min(max(dt$TIME[dt$cond=="ex vivo"]), max(dt$TIME[dt$cond=="in vivo"])))
xd  <- seq(ov[1], ov[2], length.out = 300)
dif <- band(Xdiff(xd))
fig("02_difference_invivo_minus_exvivo.png", {
  plot(xd, dif$fit, type = "n", ylim = range(dif$lo, dif$hi, 0),
       xlab = "Humerothoracic elevation (x)", ylab = "in vivo - ex vivo (ST angle)",
       main = "Estimated difference +/- 95% CI (overlap only)")
  polygon(c(xd, rev(xd)), c(dif$lo, rev(dif$hi)), col = adjustcolor(BLUE, 0.25), border = NA)
  lines(xd, dif$fit, lwd = 3, col = BLUE); abline(h = 0, lty = 2, col = "grey40")
  knot_marks(KNOTS, BOUNDARY, data.frame(x = xd, fit = dif$fit), byvar = "none")
  legend("topleft", bty = "n", cex = 0.9,
         legend = c(sprintf("mean |difference| = %.2f deg", mean(abs(dif$fit))),
                    sprintf("overlap: %.0f .. %.0f deg", ov[1], ov[2])))
})

res <- m$pred$resid_ss; z <- (res - mean(res))/sd(res); nstat <- norm_stats(res)
fig("03_diagnostics.png", {
  op <- par(mfrow = c(2,2), mar = c(4,4,3,1)); on.exit(par(op), add = TRUE)
  plot(m$pred$pred_ss, res, pch = 16, cex = 0.4, col = adjustcolor("grey20", 0.3),
       xlab = "fitted", ylab = "residual", main = "Residuals vs fitted"); abline(h = 0, col = "red")
  qqnorm(z, pch = 16, cex = 0.4, col = adjustcolor("grey20", 0.3),
         main = "Normal Q-Q (standardised)"); qqline(z, col = "red", lwd = 2)
  hist(z, breaks = 60, freq = FALSE, border = NA, col = "grey80", xlab = "standardised residual",
       main = sprintf("skew %.2f | kurtosis %.2f", nstat["skew"], nstat["kurtosis"]))
  curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
  plot(dt$TIME, res, pch = 16, cex = 0.4, col = adjustcolor(COL[as.character(dt$cond)], 0.35),
       xlab = "elevation (x)", ylab = "residual", main = "Residuals vs x")
  knot_marks(KNOTS, BOUNDARY); abline(h = 0, col = "red")
})

# every random effect, not just the intercept
RE <- m$predRE
fig("04_random_effects.png", {
  nre <- ncol(RE) - 1
  op <- par(mfrow = c(2, ceiling(nre/2)), mar = c(4,4,3,1)); on.exit(par(op), add = TRUE)
  for (j in 2:ncol(RE)) {
    qqnorm(RE[[j]], pch = 16, cex = 0.8, col = adjustcolor(BLUE, 0.7),
           main = sprintf("%s (sd %.2f)", names(RE)[j], sd(RE[[j]])), ylab = "random effect")
    qqline(RE[[j]], col = "red", lwd = 2)
  }
})

acf_res <- acf_curve(res, dt$IDnum, lag.max = 40)
acf03   <- vt$resid_acf1[vt$variant == "R1_intercept"]
fig("05_autocorrelation.png", {
  plot(0:(length(acf_res)-1), acf_res, type = "h", lwd = 2, ylim = range(0, acf_res, 1),
       xlab = "lag (points within a shoulder)", ylab = "residual ACF",
       main = sprintf("Residual ACF after a random CURVE: lag-1 = %.3f", acf_res[2]))
  abline(h = 0)
  abline(h = acf03, lty = 2, col = "red")
  legend("topright", bty = "n", cex = 0.85, lty = c(2, NA), col = c("red", NA),
         legend = c(sprintf("random intercept only: %.3f", acf03),
                    sprintf("random spline: %.3f", acf_res[2])))
})

fig("06_variant_comparison.png", {
  o <- vt[vt$converged, ]
  op <- par(mfrow = c(1,2), mar = c(7,4,3,1)); on.exit(par(op), add = TRUE)
  bp <- barplot(o$resid_acf1, names.arg = sub("^R[0-9]_", "", o$variant), las = 2,
                col = ifelse(o$variant == best_id, BLUE, "grey70"), ylim = c(0, 1),
                ylab = "lag-1 residual ACF", main = "Autocorrelation left over\n(lower is better)")
  text(bp, o$resid_acf1, sprintf("%.3f", o$resid_acf1), pos = 3, xpd = TRUE, cex = 0.85)
  bp2 <- barplot(o$BIC, names.arg = sub("^R[0-9]_", "", o$variant), las = 2,
                 col = ifelse(o$variant == best_id, BLUE, "grey70"),
                 ylab = "BIC (lower is better)", main = "BIC")
  text(bp2, o$BIC, round(o$BIC), pos = 3, xpd = TRUE, cex = 0.85)
})

# -----------------------------------------------------------------------------
# 4a. DOES TRIMMING TO THE x-OVERLAP RESCUE THE RESIDUAL NORMALITY?
#
#     Outside the range where BOTH conditions have data, one condition's curve is
#     pure extrapolation and nothing constrains it — a natural home for the
#     extreme residuals driving the kurtosis. Trimming removes those points
#     rather than merely re-describing them. (Moving the BOUNDARY KNOTS would
#     not do this: that changes the basis, not which points are fitted.)
#
#     The range is computed from the data, never hard-coded, so this carries
#     unchanged to the other 71 joint x motion x DoF cells, whose overlaps differ.
#
#     NOTE: BIC is NOT comparable between the trimmed and untrimmed fits — they
#     are different datasets. Only the residual diagnostics are comparable.
# -----------------------------------------------------------------------------
overlap_range <- function(df, xvar = "TIME", byvar = "cond") {
  rr <- tapply(df[[xvar]], droplevels(df[[byvar]]), range)
  c(max(sapply(rr, `[`, 1)), min(sapply(rr, `[`, 2)))
}
ovr <- overlap_range(dt)
dtr <- dt[dt$TIME >= ovr[1] & dt$TIME <= ovr[2], ]
dtr$ID <- droplevels(dtr$ID); dtr$IDnum <- as.integer(dtr$ID)

cat("\n================ x-OVERLAP TRIM ================\n")
cat(sprintf("overlap of the two conditions: %.1f .. %.1f deg\n", ovr[1], ovr[2]))
cat(sprintf("rows %d -> %d (%.1f%% kept) | shoulders %d -> %d\n",
            nrow(dt), nrow(dtr), 100*nrow(dtr)/nrow(dt),
            nlevels(dt$ID), nlevels(dtr$ID)))
cat(sprintf("model: UNCHANGED — K = %d, interior knots %s, boundary %s\n",
            K, paste(KNOTS, collapse=","), paste(BOUNDARY, collapse=",")))

# SAME basis as the full-range fit — same K, same interior knots, same boundary
# knots. Only the DATA changes, so the comparison below is a single-factor test.
# (Moving the boundary knots to the overlap as well was tried and is worse:
#  on the full data it takes kurtosis 12.55 -> 14.84, and on the trimmed data
#  5.25 -> 5.51. Forcing linearity over observed data hurts the tails.)
Btr <- ns(dtr$TIME, knots = KNOTS, Boundary.knots = BOUNDARY)
nstr <- paste0("ns", 1:ncol(Btr))
for (j in seq_along(nstr)) dtr[[nstr[j]]] <- Btr[, j]
ftr <- as.formula(paste("Y ~", paste(nstr, collapse=" + "), "+ cond +",
                        paste(paste0("cond:", nstr), collapse=" + ")))
rtr <- as.formula(paste("~", paste(nstr, collapse = " + ")))

t0 <- Sys.time()
m_tr <- try(hlme(fixed = ftr, random = rtr, subject = "IDnum", ng = 1, idiag = TRUE,
                 data = dtr, verbose = FALSE, nproc = NPROC), silent = TRUE)
tr_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
trim_note <- if (inherits(m_tr, "try-error") || m_tr$conv != 1) {
  "trimmed fit did not converge."
} else {
  rt <- m_tr$pred$resid_ss; nst <- norm_stats(rt)
  full <- vt[vt$variant == "R2_spline_diag", ]
  sprintf(paste0(
    "trimmed to %.1f-%.1f deg (%d rows, %.0f%% kept), converged in %.0f s\n",
    "                     full range      trimmed\n",
    "  resid sd           %8.3f     %8.3f\n",
    "  excess kurtosis    %8.2f     %8.2f\n",
    "  skew               %8.2f     %8.2f\n",
    "  lag-1 resid ACF    %8.3f     %8.3f\n",
    "  (BIC not comparable: different datasets)"),
    ovr[1], ovr[2], nrow(dtr), 100*nrow(dtr)/nrow(dt), tr_secs,
    full$resid_sd, sd(rt), full$kurtosis, nst["kurtosis"],
    full$skew, nst["skew"], full$resid_acf1, acf_within1(rt, dtr$IDnum))
}
cat(trim_note, "\n")
write.csv(data.frame(
  setting = c("full range", "trimmed to overlap"),
  x_lo = c(min(dt$TIME), ovr[1]), x_hi = c(max(dt$TIME), ovr[2]),
  n_rows = c(nrow(dt), nrow(dtr)),
  resid_sd = c(vt$resid_sd[vt$variant=="R2_spline_diag"],
               if (inherits(m_tr,"try-error")) NA else sd(m_tr$pred$resid_ss)),
  kurtosis = c(vt$kurtosis[vt$variant=="R2_spline_diag"],
               if (inherits(m_tr,"try-error")) NA else norm_stats(m_tr$pred$resid_ss)["kurtosis"]),
  resid_acf1 = c(vt$resid_acf1[vt$variant=="R2_spline_diag"],
                 if (inherits(m_tr,"try-error")) NA else acf_within1(m_tr$pred$resid_ss, dtr$IDnum))),
  file.path(OUT, "00_overlap_trim.csv"), row.names = FALSE)

if (!inherits(m_tr, "try-error") && m_tr$conv == 1) {
  fig("13_overlap_trim_diagnostics.png", {
    op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
    for (nm in c("full range", "trimmed to overlap")) {
      rr2 <- if (nm == "full range") m$pred$resid_ss else m_tr$pred$resid_ss
      zz <- (rr2 - mean(rr2))/sd(rr2); ks <- mean(zz^4) - 3
      qqnorm(zz, pch = 16, cex = 0.35, col = adjustcolor("grey20", 0.3),
             main = sprintf("%s — Q-Q", nm))
      qqline(zz, col = "red", lwd = 2)
      hist(zz, breaks = 60, freq = FALSE, border = NA, col = "grey80",
           xlab = "standardised residual",
           main = sprintf("%s — kurtosis %.2f", nm, ks))
      curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
    }
  }, w = 1200, h = 850, res = 130)
}

# -----------------------------------------------------------------------------
# 4b. df SWEEP — does a smaller basis rescue the residual normality?
#     The kurtosis blow-up is caused by the random CURVE absorbing almost all of
#     each shoulder's trajectory, leaving crumbs as eps_ij. A smaller basis gives
#     each shoulder fewer random effects, so it should absorb LESS and leave a
#     bigger, better-behaved residual. This sweep tests that directly.
#
#     Knots follow one rule for every K: quantiles of the fitted x, rounded to the
#     nearest 10 deg. At K = 5 that reproduces the main model's 40/70/100/130.
#     The diagonal structure is used throughout — the unstructured variant costs
#     40x the runtime for no measurable gain (see 00_variants.csv).
# -----------------------------------------------------------------------------
knots_for <- function(KK)
  unique(round(as.numeric(quantile(dt$TIME, seq(0, 1, length.out = KK + 1)[2:KK]))/10)*10)

sw <- list(); swfit <- list()
for (kk in KS) {
  kn <- knots_for(kk)
  Bk <- ns(dt$TIME, knots = kn, Boundary.knots = BOUNDARY)
  dk <- dt; nsk <- paste0("ns", 1:ncol(Bk))
  for (j in seq_along(nsk)) dk[[nsk[j]]] <- Bk[, j]
  fk <- as.formula(paste("Y ~", paste(nsk, collapse=" + "), "+ cond +",
                         paste(paste0("cond:", nsk), collapse=" + ")))
  rk <- as.formula(paste("~", paste(nsk, collapse = " + ")))
  cat(sprintf("\n--- df sweep: K = %d, knots %s ---\n", kk, paste(kn, collapse=",")))
  flush.console()
  t0 <- Sys.time()
  mk <- try(hlme(fixed = fk, random = rk, subject = "IDnum", ng = 1, idiag = TRUE,
                 data = dk, verbose = FALSE, nproc = NPROC), silent = TRUE)
  sk <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  okk <- !inherits(mk, "try-error") && !is.null(mk$conv) && mk$conv == 1 && all(is.finite(mk$best))
  if (!okk) {
    cat("  NOT usable\n")
    sw[[length(sw)+1]] <- data.frame(K = kk, n_random = ncol(Bk)+1, knots = paste(kn, collapse="/"),
      converged = FALSE, n_param = NA, BIC = NA, resid_acf1 = NA, resid_sd = NA,
      skew = NA, kurtosis = NA, wald_pass = NA, wald_total = NA,
      wald_pass_ref = NA, wald_total_ref = NA, wald_pass_diff = NA, wald_total_diff = NA,
      seconds = sk)
    next
  }
  invisible(capture.output(tk <- summary(mk)))
  rr <- mk$pred$resid_ss; nsx <- norm_stats(rr); wk <- abs(tk[,"coef"]/tk[,"Se"])
  wrk <- wald_from(tk); isrefk <- wrk$curve == "reference (ex vivo)"
  cat(sprintf("  conv=1 %.0fs | BIC %.0f | ACF %.3f | resid sd %.3f | kurt %.2f | Wald %d/%d (ref %d/%d, diff %d/%d)\n",
              sk, mk$BIC, acf_within1(rr, dt$IDnum), sd(rr), nsx["kurtosis"],
              sum(wk >= WALD_CRIT), length(wk),
              sum(wrk$passed[isrefk]), sum(isrefk), sum(wrk$passed[!isrefk]), sum(!isrefk)))
  swfit[[as.character(kk)]] <- list(m = mk, B = Bk, knots = kn, tab = tk)
  sw[[length(sw)+1]] <- data.frame(K = kk, n_random = ncol(Bk)+1, knots = paste(kn, collapse="/"),
    converged = TRUE, n_param = length(mk$best), BIC = mk$BIC,
    resid_acf1 = acf_within1(rr, dt$IDnum), resid_sd = sd(rr),
    skew = nsx["skew"], kurtosis = nsx["kurtosis"],
    wald_pass = sum(wk >= WALD_CRIT), wald_total = length(wk),
    wald_pass_ref  = sum(wrk$passed[isrefk]),  wald_total_ref  = sum(isrefk),
    wald_pass_diff = sum(wrk$passed[!isrefk]), wald_total_diff = sum(!isrefk),
    seconds = sk)
}
swt <- do.call(rbind, sw); rownames(swt) <- NULL
write.csv(swt, file.path(OUT, "00_df_sweep.csv"), row.names = FALSE)
cat("\n================ df SWEEP (diagonal random spline) ================\n")
print(swt, digits = 4)

# per-coefficient Wald table for EVERY df, stacked with a K column
wald_by_df <- do.call(rbind, lapply(KS, function(kk) {
  e <- swfit[[as.character(kk)]]
  if (is.null(e)) return(NULL)
  cbind(K = kk, n_fixed = 2*kk + 2, retained_df = (kk == K), wald_from(e$tab))
}))
write.csv(wald_by_df, file.path(OUT, "00_wald_tests_by_df.csv"), row.names = FALSE)
cat("\n================ WALD PASSED PER CURVE, PER df ================\n")
print(swt[, c("K","n_random","wald_pass_ref","wald_total_ref",
              "wald_pass_diff","wald_total_diff","wald_pass","wald_total")], row.names = FALSE)

# Curve builder for ANY of the swept models — same contrasts as the main model,
# but with that K's basis and its own 2K+2 fixed effects.
cm_for <- function(mm, BB, KK) {
  nf  <- 2*KK + 2
  bfk <- mm$best[1:nf]; Vfk <- VarCov(mm)[1:nf, 1:nf]
  list(
    Xrow  = function(x, inv) { b <- predict(BB, x); cbind(1, b, inv, b*inv) },
    Xdiff = function(x) { b <- predict(BB, x); cbind(0, matrix(0, length(x), KK), 1, b) },
    band  = function(X) { f <- as.vector(X %*% bfk)
      s <- sqrt(pmax(0, rowSums((X %*% Vfk) * X)))
      data.frame(fit = f, lo = f - 1.96*s, hi = f + 1.96*s) })
}

# --- the four by-df views, mirroring iteration 03 --------------------------
fig("07_by_df_population_curves.png", {
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
  for (kk in KS) {
    e <- swfit[[as.character(kk)]]
    if (is.null(e)) { plot.new(); title(main = sprintf("K=%d — not usable", kk)); next }
    cmk <- cm_for(e$m, e$B, kk); row <- swt[swt$K == kk, ]
    ck <- do.call(rbind, lapply(levels(dt$cond), function(cc) {
      rr <- range(dt$TIME[dt$cond == cc]); xc <- seq(rr[1], rr[2], length.out = 300)
      cbind(cond = cc, x = xc, cmk$band(cmk$Xrow(xc, if (cc == "in vivo") 1 else 0)))
    }))
    plot(dt$TIME, dt$Y, col = adjustcolor(COL[as.character(dt$cond)], 0.18), pch = 16, cex = 0.35,
         xlab = "elevation (x)", ylab = "ST angle (y)",
         main = sprintf("K=%d: %d rand.eff | BIC %.0f | Wald %d/%d",
                        kk, row$n_random, row$BIC, row$wald_pass, row$wald_total),
         cex.main = if (kk == K) 1.05 else 0.95, font.main = if (kk == K) 2 else 1)
    for (cn in unique(ck$cond)) { g <- ck[ck$cond == cn, ]
      polygon(c(g$x, rev(g$x)), c(g$lo, rev(g$hi)), col = adjustcolor(COL[cn], 0.30), border = NA)
      lines(g$x, g$fit, col = COL[cn], lwd = 2.6) }
    knot_marks(e$knots, BOUNDARY, ck, cex = 0.95)
    if (kk == K) { box(lwd = 2.5, col = BLUE)
      legend("bottomleft", "RETAINED df", bty = "n", cex = 0.8, text.col = BLUE, text.font = 2) }
    if (kk == KS[1]) legend("topleft", names(COL), col = COL, lwd = 2.6, bty = "n", cex = 0.8)
  }
})

# shared y-axis: the point is to compare the SIZE of the difference across df
difs4 <- lapply(KS, function(kk) { e <- swfit[[as.character(kk)]]
  if (is.null(e)) return(NULL)
  cmk <- cm_for(e$m, e$B, kk); cbind(x = xd, cmk$band(cmk$Xdiff(xd))) })
names(difs4) <- as.character(KS)
ylim4 <- range(0, unlist(lapply(difs4, function(z) if (is.null(z)) NULL else c(z$lo, z$hi))))

fig("08_by_df_difference.png", {
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
  for (kk in KS) {
    z <- difs4[[as.character(kk)]]
    if (is.null(z)) { plot.new(); title(main = sprintf("K=%d — not usable", kk)); next }
    sig <- (z$lo > 0) | (z$hi < 0)
    plot(z$x, z$fit, type = "n", ylim = ylim4, xlab = "elevation (x)",
         ylab = "in vivo - ex vivo (deg)",
         main = sprintf("K=%d: mean |diff| %.2f deg", kk, mean(abs(z$fit))),
         cex.main = if (kk == K) 1.05 else 0.95, font.main = if (kk == K) 2 else 1)
    polygon(c(z$x, rev(z$x)), c(z$lo, rev(z$hi)), col = adjustcolor(BLUE, 0.25), border = NA)
    abline(h = 0, lty = 2, col = "grey40"); lines(z$x, z$fit, lwd = 2.6, col = BLUE)
    if (any(sig)) points(z$x[sig], rep(ylim4[1], sum(sig)), pch = 15, cex = 0.35, col = BLUE)
    knot_marks(swfit[[as.character(kk)]]$knots, BOUNDARY,
               data.frame(x = z$x, fit = z$fit), byvar = "none", cex = 0.95)
    if (kk == K) { box(lwd = 2.5, col = BLUE)
      legend("bottomright", "RETAINED df", bty = "n", cex = 0.8, text.col = BLUE, text.font = 2) }
    if (kk == KS[1]) legend("topleft", "band excludes 0", bty = "n", cex = 0.75,
                            pch = 15, col = BLUE)
  }
})

# Does a smaller basis rescue normality? One row per K: Q-Q and histogram.
fig("09_by_df_diagnostics.png", {
  op <- par(mfrow = c(length(KS), 2), mar = c(4, 4, 2.6, 1)); on.exit(par(op), add = TRUE)
  for (kk in KS) {
    e <- swfit[[as.character(kk)]]
    if (is.null(e)) { plot.new(); title(main = sprintf("K=%d — not usable", kk)); plot.new(); next }
    rr <- e$m$pred$resid_ss; zz <- (rr - mean(rr))/sd(rr)
    row <- swt[swt$K == kk, ]
    tg <- sprintf("K=%d (%d random effects)", kk, row$n_random)
    qqnorm(zz, pch = 16, cex = 0.35, col = adjustcolor("grey20", 0.3),
           main = sprintf("%s — normal Q-Q", tg)); qqline(zz, col = "red", lwd = 2)
    hist(zz, breaks = 60, freq = FALSE, border = NA, col = "grey80",
         xlab = "standardised residual",
         main = sprintf("%s — skew %.2f | kurtosis %.2f", tg, row$skew, row$kurtosis))
    curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
  }
}, w = 1200, h = 1500, res = 130)

# The random INTERCEPT at every df — the one component present in all four models,
# so it is the comparable one. Its sd barely moves; the higher-order random
# effects are what change.
fig("10_by_df_random_effects.png", {
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
  lim <- range(unlist(lapply(KS, function(kk) { e <- swfit[[as.character(kk)]]
    if (is.null(e)) NULL else e$m$predRE[[2]] })), na.rm = TRUE)
  for (kk in KS) {
    e <- swfit[[as.character(kk)]]
    if (is.null(e)) { plot.new(); title(main = sprintf("K=%d — not usable", kk)); next }
    b0 <- e$m$predRE[[2]]; row <- swt[swt$K == kk, ]
    qqnorm(b0, pch = 16, cex = 0.8, col = adjustcolor(BLUE, 0.7), ylim = lim,
           ylab = "random intercept b_i0",
           main = sprintf("K=%d: sd(b_i0) %.2f | %d rand.eff", kk, sd(b0), row$n_random),
           cex.main = if (kk == K) 1.05 else 0.95, font.main = if (kk == K) 2 else 1)
    qqline(b0, col = "red", lwd = 2)
    if (kk == K) { box(lwd = 2.5, col = BLUE)
      legend("bottomright", "RETAINED df", bty = "n", cex = 0.8, text.col = BLUE, text.font = 2) }
  }
})

# Which coefficients survive, at every df. A grid of Wald statistics: one column
# per df, one row per term, so a term that only "passes" at some df is obvious.
fig("11_by_df_wald.png", {
  op <- par(mar = c(4.5, 11, 3.5, 2)); on.exit(par(op), add = TRUE)
  # order: all reference-curve terms first (intercept, ns1..nsK), then all
  # difference-curve terms — otherwise the separator line below is meaningless,
  # because unique() returns them interleaved in the order each K introduces them
  tt <- unique(wald_by_df$term)
  ord <- function(x) {
    isc <- grepl("cond", x)
    num <- suppressWarnings(as.numeric(sub("^ns([0-9]+).*$", "\\1", x)))
    num[is.na(num)] <- 0                      # intercept / bare cond term sort first
    order(isc, num)
  }
  terms_all <- tt[ord(tt)]
  ok_k <- sort(unique(wald_by_df$K))
  M <- matrix(NA_real_, length(terms_all), length(ok_k),
              dimnames = list(terms_all, paste0("K=", ok_k)))
  for (i in seq_along(ok_k)) { w <- wald_by_df[wald_by_df$K == ok_k[i], ]
    M[match(w$term, terms_all), i] <- abs(w$wald) }
  plot(NA, xlim = c(0.5, length(ok_k)+0.5), ylim = c(0.5, length(terms_all)+0.5),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "",
       main = sprintf("|Wald| per coefficient and df  (pass = >= %.2f)", WALD_CRIT))
  axis(1, at = seq_along(ok_k), labels = colnames(M))
  axis(2, at = seq_along(terms_all), labels = terms_all, las = 2, cex.axis = 0.8)
  for (i in seq_along(ok_k)) for (j in seq_along(terms_all)) {
    v <- M[j, i]; if (is.na(v)) { rect(i-.5, j-.5, i+.5, j+.5, col = "grey93", border = "white"); next }
    pass <- v >= WALD_CRIT
    rect(i-.5, j-.5, i+.5, j+.5, border = "white",
         col = if (pass) adjustcolor(BLUE, min(1, 0.18 + 0.10*log1p(v))) else "grey88")
    text(i, j, sprintf("%.1f", v), cex = 0.75,
         col = if (pass && v > 8) "white" else "black", font = if (pass) 2 else 1)
  }
  abline(h = sum(!grepl("cond", terms_all)) + 0.5, lwd = 2)   # reference | difference
  mtext("reference curve below the line, difference curve above", side = 3, cex = 0.75, line = 0.2)
}, w = 1000, h = 850, res = 130)

fig("12_by_df_tradeoff.png", {
  o <- swt[swt$converged, ]
  op <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1)); on.exit(par(op), add = TRUE)
  plot(o$K, o$kurtosis, type = "b", pch = 16, lwd = 2, col = BLUE, xlab = "spline df (K)",
       ylab = "excess kurtosis", main = "Residual normality\n(0 = normal)")
  abline(h = 0, lty = 3, col = "grey50")
  plot(o$K, o$resid_acf1, type = "b", pch = 16, lwd = 2, col = BLUE, ylim = c(0, 1),
       xlab = "spline df (K)", ylab = "lag-1 residual ACF",
       main = "Autocorrelation left\n(lower is better)")
  plot(o$K, o$resid_sd, type = "b", pch = 16, lwd = 2, col = BLUE, xlab = "spline df (K)",
       ylab = "residual sd (deg)", main = "How much is left\nfor eps_ij")
})

# -----------------------------------------------------------------------------
# 5. DIGEST
# -----------------------------------------------------------------------------
sink(file.path(OUT, "00_results_summary.txt"))
cat("Scapulothoracic rhythm — iteration 04: RANDOM spline\n")
cat("====================================================\n\n")
cat(sprintf("Data: %d obs thinned to %d (<=%d/shoulder), %d shoulders.\n",
            nrow(d), nrow(dt), CAP, nlevels(dt$ID)))
cat(sprintf("Basis: ns with EXPLICIT knots %s, boundary %s (%d columns).\n",
            paste(KNOTS, collapse=","), paste(BOUNDARY, collapse=","), K))
cat(sprintf("RETAINED: %s\n\n", fits[[best_id]]$lab))
cat("VARIANT COMPARISON (resid_acf1 is the metric this iteration exists to improve):\n")
print(vt, digits = 4)
cat(sprintf("\n  iteration 03 (random intercept) left lag-1 ACF = %.3f\n", acf03))
cat(sprintf("  this iteration                  leaves lag-1 ACF = %.3f\n", acf_res[2]))
cat("\nFIXED EFFECTS (Wald = coef/se):\n"); print(wald_tab, digits = 4)
cat(sprintf("\nWALD PASSED: reference %s | difference %s\n", ref_lab, diff_lab))
cat("\nVARIANCE COMPONENTS:\n")
print(allp[allp$block == "variance_component", c("parameter","estimate","se")], digits = 4)
cat("\nRESIDUAL NORMALITY:\n")
cat(sprintf("  skew %.3f | excess kurtosis %.3f | within +/-1.96 sd %.1f%%\n",
            nstat["skew"], nstat["kurtosis"], 100*nstat["frac_within_1.96"]))
cat("\ndf SWEEP (diagonal random spline — does a smaller basis rescue normality?):\n")
print(swt, digits = 4)
cat("\nWALD PASSED PER CURVE, PER df (see 00_wald_tests_by_df.csv for every coefficient):\n")
print(swt[, c("K","n_random","wald_pass_ref","wald_total_ref",
              "wald_pass_diff","wald_total_diff")], row.names = FALSE)
cat("\nx-OVERLAP TRIM (does removing the extrapolated tails rescue normality?):\n")
cat(trim_note, "\n")
cat("\nNOTE ON gamma_0: the boundary knot is at 0 deg, so 'condin vivo' is now the\n")
cat("  in-vivo minus ex-vivo difference AT 0 deg of elevation — an interpretable\n")
cat("  quantity, unlike iteration 03 where the basis vanished at -3.46 deg.\n")
cat("\nCAVEAT: condition remains confounded with source study — descriptive, not causal.\n")
sink()

cat(sprintf("\nDONE. See %s/\n", OUT))
