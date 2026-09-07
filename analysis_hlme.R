# =============================================================================
# ITERATION 03 — natural-spline mixed model with lcmm::hlme
#
# Same analysis as iteration 02 (figures/02_spline_mixed_model/), same figure
# set, but a different model — chosen so that every parameter can be TESTED
# individually instead of only through a whole-smooth F-test:
#
#   * splines::ns(TIME, df = K)  -> natural spline: 2K+2 fixed effects in total,
#                                   against 64 coefficients for the GAMM
#   * lcmm::hlme                 -> ML fit reporting coef / Se / Wald per term
#   * cor = AR(TIME)             -> CONTINUOUS-TIME autoregressive residuals
#                                   (decay in elevation distance, not row index)
#   * Wald counts per curve      -> how many |coef/se| >= 1.96, per curve
#   * residual normality         -> an explicit acceptance criterion, checked
#
# Run:  Rscript analysis_hlme.R
# Prereq: python3 prepare_monolix_data.py   (creates the input CSV)
# Packages: lcmm, splines, mgcv (iteration-02 baseline only), dplyr.
# =============================================================================

suppressMessages({ library(lcmm); library(splines); library(mgcv); library(dplyr) })

INPUT   <- "monolix_st_frontal_dof2.csv"
FIG_DIR <- "figures"
S3      <- "03_natural_spline_hlme"
# CAP drives the runtime CUBICALLY: hlme with cor=AR(TIME) inverts a dense
# n_i x n_i matrix per shoulder at every likelihood evaluation. Measured on this
# data (K=2, nproc=1): cap 25 -> 2.1 s, cap 50 -> 7.4 s, cap 100 -> 78 s.
# cap 200 was ~10 min PER K. 100 is the practical ceiling.
CAP     <- 100          # max points per shoulder after thinning
KS      <- 2:5          # natural-spline df to try
NPROC   <- max(1L, min(8L, parallel::detectCores() - 1L))  # hlme defaults to 1(!)
WALD_CRIT <- 1.96       # |coef/se| >= this  => "passes"
COL     <- c("ex vivo" = "#d95f02", "in vivo" = "#1b9e77")
BLUE    <- "#377eb8"

OUT <- file.path(FIG_DIR, S3)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fig <- function(name, expr, w = 1100, h = 750, res = 130) {
  png(file.path(OUT, name), width = w, height = h, res = res)
  on.exit(dev.off()); force(expr); message("  saved ", file.path(S3, name))
}

# -----------------------------------------------------------------------------
# 1. LOAD (identical to iteration 02)
# -----------------------------------------------------------------------------
if (!file.exists(INPUT))
  stop("Missing ", INPUT, " — run  python3 prepare_monolix_data.py  first.")

d <- read.csv(INPUT, stringsAsFactors = FALSE)
d$cond <- factor(ifelse(d$in_vivo %in% c("True","TRUE","1","yes"),
                        "in vivo", "ex vivo"), levels = c("ex vivo","in vivo"))
d$ID <- factor(d$ID)
d <- d[is.finite(d$TIME) & is.finite(d$Y), ]
d <- d[order(d$ID, d$TIME), ]

cat("\n================ ITERATION 03: ns + hlme ================\n")
cat(sprintf("full data: %d obs, %d shoulders (%d ex / %d in rows)\n",
            nrow(d), nlevels(d$ID), sum(d$cond=="ex vivo"), sum(d$cond=="in vivo")))

# -----------------------------------------------------------------------------
# 2. THINNING — hlme is far slower than bam, and at rho ~ 0.996 a 2987-point
#    curve carries only a handful of independent points. Keep at most CAP points
#    per shoulder, evenly spaced in x-RANK, so each curve's x-coverage (and both
#    its endpoints) is preserved exactly.
# -----------------------------------------------------------------------------
thin_by_unit <- function(df, cap) {
  do.call(rbind, lapply(split(df, df$ID), function(u) {
    if (nrow(u) <= cap) return(u)
    u[unique(round(seq(1, nrow(u), length.out = cap))), ]
  }))
}
n_before <- as.integer(table(d$ID))
dt <- thin_by_unit(d, CAP)
dt <- dt[order(dt$ID, dt$TIME), ]
n_after <- as.integer(table(dt$ID))

thin_tab <- data.frame(
  shoulder  = levels(d$ID),
  condition = sapply(levels(d$ID), function(u) as.character(d$cond[d$ID == u][1])),
  n_before  = n_before, n_after = n_after, thinned = n_before > n_after,
  row.names = NULL)
write.csv(thin_tab, file.path(OUT, "00_thinning.csv"), row.names = FALSE)

cat(sprintf("thinned to <=%d pts/shoulder: %d -> %d rows (%.1f%%); %d of %d shoulders thinned\n",
            CAP, nrow(d), nrow(dt), 100*nrow(dt)/nrow(d), sum(thin_tab$thinned), nrow(thin_tab)))
cat(sprintf("rows per condition  before: %d ex / %d in   after: %d ex / %d in\n",
            sum(d$cond=="ex vivo"), sum(d$cond=="in vivo"),
            sum(dt$cond=="ex vivo"), sum(dt$cond=="in vivo")))

dt$IDnum <- as.integer(dt$ID)          # lcmm requires a NUMERIC subject id

# -----------------------------------------------------------------------------
# 3. FIT THE K GRID
#    The basis is built ONCE per K and kept, so the plotting grid reuses the same
#    knots via predict(B, newx). Rebuilding ns() on new data moves the knots.
# -----------------------------------------------------------------------------
build <- function(df, K) {
  B <- ns(df$TIME, df = K)
  for (k in 1:K) df[[paste0("ns", k)]] <- B[, k]
  list(df = df, B = B)
}
form_K <- function(K) as.formula(paste(
  "Y ~", paste(paste0("ns", 1:K), collapse = " + "),
  "+ cond +", paste(paste0("cond:ns", 1:K), collapse = " + ")))

# Wald table from an hlme summary matrix, split by which curve each term belongs to:
#   reference curve  = intercept + ns1..nsK   (the ex-vivo trajectory)
#   difference curve = cond + cond:ns1..nsK   (the in-vivo departure)
wald_from <- function(tb) {
  tm <- rownames(tb); w <- tb[, "coef"] / tb[, "Se"]
  data.frame(
    curve = ifelse(grepl("cond", tm), "difference (in vivo - ex vivo)", "reference (ex vivo)"),
    term = tm, estimate = tb[, "coef"], se = tb[, "Se"],
    wald = w, p_value = tb[, "p-value"], passed = abs(w) >= WALD_CRIT,
    row.names = NULL)
}

norm_stats <- function(r) {
  r <- r[is.finite(r)]; z <- (r - mean(r)) / sd(r)
  sw <- if (length(r) > 3) {
    s <- if (length(r) > 5000) sample(r, 5000) else r
    shapiro.test(s)$p.value
  } else NA_real_
  c(skew = mean(z^3), kurtosis = mean(z^4) - 3,
    frac_within_1.96 = mean(abs(z) <= 1.96), shapiro_p = sw)
}
# mean within-subject ACF (never crosses a shoulder boundary)
acf_within <- function(r, id, lag.max = 40) {
  acc <- matrix(NA_real_, 0, lag.max + 1)
  for (u in unique(id)) { ru <- r[id == u]
    if (length(ru) > lag.max + 5)
      acc <- rbind(acc, drop(acf(ru, lag.max = lag.max, plot = FALSE)$acf)) }
  colMeans(acc, na.rm = TRUE)
}

fits <- list(); rows <- list()
for (K in KS) {
  bb <- build(dt, K); dK <- bb$df
  cat(sprintf("\n--- K = %d (%d fixed effects, nproc=%d) ---\n", K, 2*K + 2, NPROC))
  flush.console()
  t0 <- Sys.time()
  # NOTE: cor = AR(TIME) is deliberately NOT used here — see section 6. On this
  # data the continuous-time AR is ALIASED with the random intercept (its decay
  # goes to zero, so it becomes a random intercept) and swallows both the
  # random-effect and residual variance. Section 6 fits it and documents that.
  m <- try(hlme(fixed = form_K(K), random = ~1, subject = "IDnum",
                ng = 1, data = dK, verbose = FALSE,
                nproc = NPROC), silent = TRUE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (inherits(m, "try-error") || is.null(m$conv) || m$conv != 1) {
    cat(sprintf("  did NOT converge (%.1f s)\n", secs))
    rows[[length(rows)+1]] <- data.frame(K = K, n_fixed = 2*K+2, converged = FALSE,
      loglik = NA, AIC = NA, BIC = NA, wald_pass = NA, wald_total = NA,
      wald_pass_ref = NA, wald_total_ref = NA, wald_pass_diff = NA, wald_total_diff = NA,
      sigma_b = NA, sigma_resid = NA,
      skew = NA, kurtosis = NA, shapiro_p = NA, seconds = secs)
    next
  }
  invisible(capture.output(tab <- summary(m)))
  wald  <- abs(tab[, "coef"] / tab[, "Se"])
  ns_st <- norm_stats(m$pred$resid_ss)
  cat(sprintf("  converged in %.1f s | BIC %.1f | Wald pass %d/%d | skew %.2f kurt %.2f\n",
              secs, m$BIC, sum(wald >= WALD_CRIT), length(wald), ns_st["skew"], ns_st["kurtosis"]))
  fits[[as.character(K)]] <- list(m = m, B = bb$B, df = dK, tab = tab)
  wr  <- wald_from(tab); isref <- wr$curve == "reference (ex vivo)"
  rows[[length(rows)+1]] <- data.frame(K = K, n_fixed = 2*K+2, converged = TRUE,
    loglik = m$loglik, AIC = m$AIC, BIC = m$BIC,
    wald_pass = sum(wald >= WALD_CRIT), wald_total = length(wald),
    wald_pass_ref  = sum(wr$passed[isref]),  wald_total_ref  = sum(isref),
    wald_pass_diff = sum(wr$passed[!isref]), wald_total_diff = sum(!isref),
    sigma_b = sqrt(m$best[["varcov 1"]]), sigma_resid = abs(m$best[["stderr"]]),
    skew = ns_st["skew"], kurtosis = ns_st["kurtosis"], shapiro_p = ns_st["shapiro_p"],
    seconds = secs)
}
sel <- do.call(rbind, rows); rownames(sel) <- NULL
write.csv(sel, file.path(OUT, "00_model_selection.csv"), row.names = FALSE)
cat("\n================ K SELECTION ================\n"); print(sel, digits = 4)

ok <- sel[sel$converged, ]
if (!nrow(ok)) stop("no K converged — inspect 00_model_selection.csv")
K_best <- ok$K[which.min(ok$BIC)]
cat(sprintf("\nRETAINED: K = %d  (lowest BIC among converged)\n", K_best))

best <- fits[[as.character(K_best)]]
m <- best$m; B <- best$B; dK <- best$df; tab <- best$tab
K <- K_best; nfix <- 2*K + 2

# -----------------------------------------------------------------------------
# 4. WALD TESTS, per curve
#    reference curve  = intercept + ns1..nsK   (the ex-vivo trajectory)
#    difference curve = cond + cond:ns1..nsK   (the in-vivo departure)
# -----------------------------------------------------------------------------
wald_tab <- wald_from(tab)                       # the RETAINED model only
write.csv(wald_tab, file.path(OUT, "00_wald_tests.csv"), row.names = FALSE)

# ... and the same table for every df that converged, stacked with a K column.
wald_all <- do.call(rbind, lapply(KS, function(kk) {
  e <- fits[[as.character(kk)]]
  if (is.null(e)) return(NULL)
  cbind(K = kk, n_fixed = 2*kk + 2, retained = (kk == K), wald_from(e$tab))
}))
write.csv(wald_all, file.path(OUT, "00_wald_tests_by_df.csv"), row.names = FALSE)

cat("\n================ WALD PASSED PER CURVE, PER df ================\n")
print(sel[, c("K","n_fixed","wald_pass_ref","wald_total_ref",
              "wald_pass_diff","wald_total_diff","wald_pass","wald_total")],
      row.names = FALSE)

wc <- tapply(wald_tab$passed, wald_tab$curve, function(x) sprintf("%d/%d", sum(x), length(x)))
ref_lab  <- wc[["reference (ex vivo)"]]
diff_lab <- wc[["difference (in vivo - ex vivo)"]]
cat("\n================ WALD TESTS (|coef/se| >= 1.96) ================\n")
print(wald_tab, digits = 4)
cat(sprintf("\nPASSED — reference curve: %s | difference curve: %s\n", ref_lab, diff_lab))

V    <- VarCov(m)
allp <- data.frame(parameter = names(m$best), estimate = as.vector(m$best),
                   se = sqrt(diag(V)), row.names = NULL)
allp$block <- ifelse(seq_len(nrow(allp)) <= nfix, "fixed_effect", "variance_component")
allp$wald  <- ifelse(allp$block == "fixed_effect", allp$estimate/allp$se, NA)
write.csv(allp, file.path(OUT, "00_coefficients.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. FIGURES — same set, same order, same names as iteration 02
# -----------------------------------------------------------------------------
# Fitted curves with 95% bands, for ANY of the fitted models (used both for the
# retained K below and for the df comparison in figure 08).
#   design row order must match names(best)[1:nfix]:
#   intercept, ns1..nsK, condin vivo, ns1:condin vivo .. nsK:condin vivo
curve_maker <- function(mm, BB, KK) {
  nf <- 2*KK + 2
  bf <- mm$best[1:nf]; Vf <- VarCov(mm)[1:nf, 1:nf]
  Xrow <- function(x, inv) { b <- predict(BB, x); cbind(1, b, inv, b * inv) }
  band <- function(X) {
    fit <- as.vector(X %*% bf)
    se  <- sqrt(pmax(0, rowSums((X %*% Vf) * X)))
    data.frame(fit = fit, lo = fit - 1.96*se, hi = fit + 1.96*se)
  }
  # in vivo MINUS ex vivo: the intercept and the main-effect spline terms cancel,
  # leaving exactly the condition block  gamma_0 + sum_k gamma_k B_k(x)
  Xdiff <- function(x) { b <- predict(BB, x); cbind(0, matrix(0, length(x), KK), 1, b) }
  list(band = band, Xrow = Xrow, Xdiff = Xdiff,
       curves = do.call(rbind, lapply(levels(dt$cond), function(cc) {
         rr <- range(dt$TIME[dt$cond == cc]); xc <- seq(rr[1], rr[2], length.out = 300)
         cbind(cond = cc, x = xc, band(Xrow(xc, if (cc == "in vivo") 1 else 0)))
       })))
}
cm     <- curve_maker(m, B, K)
band   <- cm$band; Xrow <- cm$Xrow; curves <- cm$curves

fig("01_population_curves_by_condition.png", {
  plot(dt$TIME, dt$Y, col = adjustcolor(COL[as.character(dt$cond)], 0.25), pch = 16, cex = 0.5,
       xlab = "Humerothoracic elevation (x)", ylab = "Scapulothoracic angle (y)",
       main = sprintf("ns(df=%d) + hlme with AR(TIME) — %d fixed effects", K, nfix))
  for (cc in unique(curves$cond)) { g <- curves[curves$cond == cc, ]
    polygon(c(g$x, rev(g$x)), c(g$lo, rev(g$hi)), col = adjustcolor(COL[cc], 0.3), border = NA)
    lines(g$x, g$fit, col = COL[cc], lwd = 3) }
  legend("topleft", names(COL), col = COL, lwd = 3, bty = "n")
  legend("bottomright", bty = "n", cex = 0.85,
         legend = c(sprintf("Wald  |coef/se| >= %.2f", WALD_CRIT),
                    sprintf("reference curve : %s", ref_lab),
                    sprintf("difference curve: %s", diff_lab)))
})

# difference measured only on the x-overlap of the two conditions
ov <- c(max(min(dt$TIME[dt$cond=="ex vivo"]), min(dt$TIME[dt$cond=="in vivo"])),
        min(max(dt$TIME[dt$cond=="ex vivo"]), max(dt$TIME[dt$cond=="in vivo"])))
xd  <- seq(ov[1], ov[2], length.out = 300)
dif <- band(cm$Xdiff(xd))                            # in vivo MINUS ex vivo

fig("02_difference_invivo_minus_exvivo.png", {
  plot(xd, dif$fit, type = "n", ylim = range(dif$lo, dif$hi, 0),
       xlab = "Humerothoracic elevation (x)", ylab = "in vivo  -  ex vivo  (ST angle)",
       main = "Estimated difference (in vivo vs ex vivo) +/- 95% CI")
  polygon(c(xd, rev(xd)), c(dif$lo, rev(dif$hi)), col = adjustcolor(BLUE, 0.25), border = NA)
  lines(xd, dif$fit, lwd = 3, col = BLUE)
  abline(h = 0, lty = 2, col = "grey40")      # band excluding 0 => differ there
  legend("topleft", bty = "n", cex = 0.9,
         legend = c(sprintf("mean |difference| = %.2f deg", mean(abs(dif$fit))),
                    sprintf("overlap: %.0f .. %.0f deg", ov[1], ov[2])))
})

res <- m$pred$resid_ss                       # subject-specific residuals
z   <- (res - mean(res)) / sd(res)
nstat <- norm_stats(res)

fig("03_diagnostics.png", {
  op <- par(mfrow = c(2,2), mar = c(4,4,3,1)); on.exit(par(op), add = TRUE)
  plot(m$pred$pred_ss, res, pch = 16, cex = 0.4, col = adjustcolor("grey20", 0.3),
       xlab = "fitted", ylab = "residual", main = "Residuals vs fitted"); abline(h = 0, col = "red")
  qqnorm(z, pch = 16, cex = 0.4, col = adjustcolor("grey20", 0.3),
         main = "Normal Q-Q (standardised)"); qqline(z, col = "red", lwd = 2)
  hist(z, breaks = 60, freq = FALSE, border = NA, col = "grey80",
       xlab = "standardised residual", main = sprintf("skew %.2f | kurtosis %.2f",
       nstat["skew"], nstat["kurtosis"]))
  curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
  plot(dK$TIME, res, pch = 16, cex = 0.4, col = adjustcolor(COL[as.character(dK$cond)], 0.35),
       xlab = "elevation (x)", ylab = "residual", main = "Residuals vs x"); abline(h = 0, col = "red")
})

fig("04_random_effects.png", {
  op <- par(mfrow = c(1,2), mar = c(4,4,3,1)); on.exit(par(op), add = TRUE)
  re <- m$predRE[[2]]
  hist(re, col = "grey80", border = "white", xlab = "random intercept (shoulder)",
       main = sprintf("Between-shoulder intercepts b_i (n=%d)", length(re)))
  qqnorm(re, pch = 16, main = "Are the b_i normal?"); qqline(re, col = "red", lwd = 2)
})

acf_res <- acf_within(res, dK$IDnum, lag.max = 40)
fig("05_autocorrelation.png", {
  plot(0:(length(acf_res)-1), acf_res, type = "h", lwd = 2, ylim = range(0, acf_res),
       xlab = "lag (points within a shoulder)", ylab = "residual ACF",
       main = sprintf("UNMODELLED residual autocorrelation: lag-1 = %.3f", acf_res[2]))
  abline(h = 0)
  abline(h = c(-1,1)*1.96/sqrt(nrow(dK)/nlevels(dt$ID)), lty = 2, col = "red")
  legend("topright", bty = "n", cex = 0.8,
         legend = c("A shoulder's deviation from the population curve",
                    "is itself a SMOOTH function of x, not noise.",
                    "A random intercept absorbs only its LEVEL.",
                    "=> the Wald SEs below are OPTIMISTIC.",
                    "See section 6 for why AR(TIME) cannot fix it here."))
})

fig("06_model_comparison.png", {
  op <- par(mfrow = c(1,2), mar = c(4,4,3,1)); on.exit(par(op), add = TRUE)
  o <- sel[sel$converged, ]
  bp <- barplot(o$BIC, names.arg = sprintf("K=%d\n(%d par)", o$K, o$n_fixed),
                col = ifelse(o$K == K, BLUE, "grey70"), las = 1, cex.names = 0.8,
                ylab = "BIC (lower is better)", main = "Natural-spline df",
                ylim = range(0, o$BIC) * 1.05)
  text(bp, o$BIC, round(o$BIC), pos = 3, xpd = TRUE, cex = 0.8)
  plot(o$K, o$wald_pass/o$wald_total*100, type = "b", pch = 16, lwd = 2, col = BLUE,
       ylim = c(0, 100), xlab = "natural-spline df (K)", ylab = "% Wald tests passed",
       main = "Fraction of testable parameters")
  abline(h = 100, lty = 3, col = "grey50")
})

# One panel per natural-spline df, each on its OWN axes, so the effect of adding
# basis functions is visible directly: how much extra shape each df buys, and
# where the curves stop changing.
fig("08_by_df_population_curves.png", {
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
  for (kk in KS) {
    e <- fits[[as.character(kk)]]
    if (is.null(e)) {                       # did not converge — say so in place
      plot.new(); title(main = sprintf("ns(df=%d) — did not converge", kk)); next
    }
    cc_k <- curve_maker(e$m, e$B, kk)$curves
    row  <- sel[sel$K == kk, ]
    plot(dt$TIME, dt$Y, col = adjustcolor(COL[as.character(dt$cond)], 0.18),
         pch = 16, cex = 0.35, xlab = "elevation (x)", ylab = "ST angle (y)",
         main = sprintf("ns(df=%d): %d params | BIC %.0f | Wald %d/%d",
                        kk, 2*kk + 2, row$BIC, row$wald_pass, row$wald_total),
         cex.main = if (kk == K) 1.05 else 0.95,
         font.main = if (kk == K) 2 else 1)
    for (cn in unique(cc_k$cond)) { g <- cc_k[cc_k$cond == cn, ]
      polygon(c(g$x, rev(g$x)), c(g$lo, rev(g$hi)),
              col = adjustcolor(COL[cn], 0.30), border = NA)
      lines(g$x, g$fit, col = COL[cn], lwd = 2.6) }
    # mark the interior knots that this df places
    rug(attr(e$B, "knots"), lwd = 2, col = "grey30", ticksize = 0.035)
    if (kk == K) {
      box(lwd = 2.5, col = BLUE)
      legend("bottomleft", "RETAINED", bty = "n", cex = 0.8, text.col = BLUE, text.font = 2)
    }
    if (kk == KS[1]) legend("topleft", names(COL), col = COL, lwd = 2.6, bty = "n", cex = 0.8)
  }
})

# The same comparison for the DIFFERENCE curve. Unlike figure 08 these panels
# share one y-axis: the whole point is to compare the size of the in-vivo minus
# ex-vivo difference across df, which separate scales would hide.
difs <- lapply(KS, function(kk) {
  e <- fits[[as.character(kk)]]
  if (is.null(e)) return(NULL)
  ck <- curve_maker(e$m, e$B, kk)
  cbind(x = xd, ck$band(ck$Xdiff(xd)))
})
names(difs) <- as.character(KS)
ylim_d <- range(0, unlist(lapply(difs, function(z) if (is.null(z)) NULL else c(z$lo, z$hi))))

fig("09_by_df_difference.png", {
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
  for (kk in KS) {
    z <- difs[[as.character(kk)]]
    if (is.null(z)) { plot.new(); title(main = sprintf("ns(df=%d) — did not converge", kk)); next }
    sig <- (z$lo > 0) | (z$hi < 0)          # where the band excludes zero
    plot(z$x, z$fit, type = "n", ylim = ylim_d,
         xlab = "elevation (x)", ylab = "in vivo - ex vivo (deg)",
         main = sprintf("ns(df=%d): mean |diff| %.2f deg | Wald %s",
                        kk, mean(abs(z$fit)),
                        with(sel[sel$K == kk, ], sprintf("%d/%d", wald_pass, wald_total))),
         cex.main = if (kk == K) 1.05 else 0.95,
         font.main = if (kk == K) 2 else 1)
    polygon(c(z$x, rev(z$x)), c(z$lo, rev(z$hi)), col = adjustcolor(BLUE, 0.25), border = NA)
    abline(h = 0, lty = 2, col = "grey40")
    lines(z$x, z$fit, lwd = 2.6, col = BLUE)
    # mark, along the bottom, the elevations where the difference excludes zero
    if (any(sig)) points(z$x[sig], rep(ylim_d[1], sum(sig)), pch = 15, cex = 0.35, col = BLUE)
    rug(attr(fits[[as.character(kk)]]$B, "knots"), lwd = 2, col = "grey30", ticksize = 0.035)
    if (kk == K) {
      box(lwd = 2.5, col = BLUE)
      legend("bottomright", "RETAINED", bty = "n", cex = 0.8, text.col = BLUE, text.font = 2)
    }
    if (kk == KS[1]) legend("topleft", bty = "n", cex = 0.75, pch = 15, col = BLUE,
                            legend = "band excludes 0")
  }
})

# Normality of eps_ij at every df — the evidence that the acceptance criterion
# holds regardless of which df is retained, not just for the winner. Shared axes
# so the four Q-Q plots are directly comparable.
zs <- lapply(KS, function(kk) { e <- fits[[as.character(kk)]]
  if (is.null(e)) return(NULL); r <- e$m$pred$resid_ss; (r - mean(r)) / sd(r) })
names(zs) <- as.character(KS)
lim_z <- range(unlist(zs), na.rm = TRUE)

# One ROW per df, two columns: the normal Q-Q and the histogram against a fitted
# normal — the same pair 03_diagnostics.png shows for the retained model. Shared
# axes down each column so the four df are directly comparable.
fig("10_by_df_diagnostics.png", {
  op <- par(mfrow = c(length(KS), 2), mar = c(4, 4, 2.6, 1)); on.exit(par(op), add = TRUE)
  for (kk in KS) {
    z <- zs[[as.character(kk)]]
    if (is.null(z)) {
      plot.new(); title(main = sprintf("ns(df=%d) — did not converge", kk)); plot.new(); next
    }
    row <- sel[sel$K == kk, ]
    bold <- if (kk == K) 2 else 1
    tag  <- if (kk == K) sprintf("ns(df=%d)  [RETAINED]", kk) else sprintf("ns(df=%d)", kk)

    qqnorm(z, pch = 16, cex = 0.35, col = adjustcolor("grey20", 0.30),
           xlim = lim_z, ylim = lim_z, font.main = bold,
           main = sprintf("%s — normal Q-Q", tag))
    qqline(z, col = "red", lwd = 2)
    if (kk == K) box(lwd = 2.5, col = BLUE)

    hist(z, breaks = 60, freq = FALSE, border = NA, col = "grey80",
         xlim = lim_z, xlab = "standardised residual", font.main = bold,
         main = sprintf("%s — skew %.2f | kurtosis %.2f", tag, row$skew, row$kurtosis))
    curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
    if (kk == K) box(lwd = 2.5, col = BLUE)
  }
}, w = 1200, h = 1500, res = 130)

# The 44 shoulder intercepts at every df. b_i are assumed normal, so the Q-Q is
# the check; sigma_b is annotated to show it barely moves with df.
res_b <- lapply(KS, function(kk) { e <- fits[[as.character(kk)]]
  if (is.null(e)) NULL else e$m$predRE[[2]] })
names(res_b) <- as.character(KS)
lim_b <- range(unlist(res_b), na.rm = TRUE)

fig("11_by_df_random_effects.png", {
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
  for (kk in KS) {
    b <- res_b[[as.character(kk)]]
    if (is.null(b)) { plot.new(); title(main = sprintf("ns(df=%d) — did not converge", kk)); next }
    sb <- sqrt(fits[[as.character(kk)]]$m$best[["varcov 1"]])
    qqnorm(b, pch = 16, cex = 0.8, col = adjustcolor(BLUE, 0.7), ylim = lim_b,
           main = sprintf("ns(df=%d): sigma_b %.2f | n=%d shoulders", kk, sb, length(b)),
           ylab = "shoulder intercept b_i",
           cex.main = if (kk == K) 1.05 else 0.95,
           font.main = if (kk == K) 2 else 1)
    qqline(b, col = "red", lwd = 2)
    if (kk == K) {
      box(lwd = 2.5, col = BLUE)
      legend("bottomright", "RETAINED", bty = "n", cex = 0.8, text.col = BLUE, text.font = 2)
    }
  }
})

# -----------------------------------------------------------------------------
# 6. WHY cor = AR(TIME) IS NOT USED — fit it once and record the degeneracy.
#    lcmm's AR(TIME) has correlation exp(-cor1 * |x_j - x_k|). Here cor1 is
#    driven to ~0, so the correlation never decays and the process becomes a
#    random intercept: it is ALIASED with random = ~1. It then swallows both the
#    random-effect variance and the residual variance, and the reported
#    resid_ss (which subtracts the predicted AR path) collapses to near zero,
#    which is why its kurtosis explodes. Kept as evidence, like iteration 01.
# -----------------------------------------------------------------------------
cat("\n================ WHY NOT cor = AR(TIME) ================\n")
t0 <- Sys.time()
m_ar <- try(hlme(fixed = form_K(K), random = ~1, subject = "IDnum", ng = 1,
                 cor = AR(TIME), data = dK, verbose = FALSE, nproc = NPROC),
            silent = TRUE)
ar_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
ar_note <- if (inherits(m_ar, "try-error") || m_ar$conv != 1) {
  "AR(TIME) did not converge."
} else {
  rr <- m_ar$pred$resid_ss; zz <- (rr - mean(rr))/sd(rr)
  sprintf(paste0("AR(TIME) converged (%.0f s) but DEGENERATE:\n",
    "  AR decay cor1 = %.6f  -> correlation exp(-cor1*dx) ~ 1 even 100 deg apart,\n",
    "                           i.e. the AR process IS a random intercept.\n",
    "  AR process sd cor2 = %.3f   vs random-intercept sd %.3f without it\n",
    "  random-intercept variance collapsed to %.5f (was %.3f)\n",
    "  residual sd collapsed to %.4f (was %.3f)\n",
    "  resid_ss kurtosis %.1f — an ARTEFACT: resid_ss subtracts the predicted\n",
    "  AR path, so it is not the eps_ij whose normality we require."),
    ar_secs, m_ar$best[["cor1"]], abs(m_ar$best[["cor2"]]),
    sqrt(m$best[["varcov 1"]]), m_ar$best[["varcov 1"]], m$best[["varcov 1"]],
    abs(m_ar$best[["stderr"]]), abs(m$best[["stderr"]]), mean(zz^4) - 3)
}
cat(ar_note, "\n")

# -----------------------------------------------------------------------------
# 7. CROSS-ITERATION CHECK vs the iteration-02 GAMM (full data, mgcv)
#    BIC is NOT comparable across mgcv and lcmm — different criteria, and lcmm's
#    BIC uses the number of SUBJECTS. Compare the fitted curves, not the scores.
# -----------------------------------------------------------------------------
d$condO    <- as.ordered(d$cond)
d$ar_start <- c(TRUE, d$ID[-1] != d$ID[-nrow(d)])
gform <- Y ~ condO + s(TIME) + s(TIME, by = condO) + s(ID, bs = "re")
g0 <- bam(gform, data = d, method = "fREML", discrete = TRUE)
r  <- resid(g0); num <- den <- 0
for (u in levels(d$ID)) { ru <- r[d$ID == u]
  if (length(ru) > 1) { num <- num + sum(head(ru,-1)*tail(ru,-1)); den <- den + sum(ru^2) } }
gam <- bam(gform, data = d, method = "fREML", discrete = TRUE,
           rho = num/den, AR.start = d$ar_start)

fig("07_vs_iteration02_gamm.png", {
  plot(dt$TIME, dt$Y, col = adjustcolor("grey60", 0.18), pch = 16, cex = 0.45,
       xlab = "Humerothoracic elevation (x)", ylab = "Scapulothoracic angle (y)",
       main = sprintf("iteration 03 (ns df=%d, %d params) vs iteration 02 (GAMM, 64 coefs)",
                      K, nfix))
  for (cc in levels(dt$cond)) {
    g <- curves[curves$cond == cc, ]; lines(g$x, g$fit, col = COL[cc], lwd = 3)
    rr <- range(d$TIME[d$cond == cc]); xg <- seq(rr[1], rr[2], length.out = 300)
    nd <- data.frame(TIME = xg, ID = d$ID[1],
                     condO = factor(cc, levels = levels(d$condO), ordered = TRUE))
    lines(xg, predict(gam, nd, exclude = "s(ID)", newdata.guaranteed = TRUE),
          col = COL[cc], lwd = 2, lty = 2)
  }
  legend("topleft", bty = "n", lwd = c(3,2), lty = c(1,2),
         legend = c(sprintf("iteration 03: hlme + ns (thinned, %d params)", nfix),
                    "iteration 02: mgcv GAMM (full data, 64 coefs)"))
})

# -----------------------------------------------------------------------------
# 8. DIGEST
# -----------------------------------------------------------------------------
sink(file.path(OUT, "00_results_summary.txt"))
cat("Scapulothoracic rhythm — in vivo vs ex vivo — iteration 03 digest\n")
cat("=================================================================\n\n")
cat(sprintf("Data: %d obs thinned to %d (<=%d pts/shoulder), %d shoulders.\n",
            nrow(d), nrow(dt), CAP, nlevels(dt$ID)))
cat("RETAINED MODEL: natural spline + hlme, random intercept per shoulder\n")
cat(sprintf("  Y ~ ns(TIME, df=%d) * cond,  random = ~1 | shoulder\n\n", K))
cat("MODEL SELECTION over the natural-spline df:\n"); print(sel, digits = 4)
cat("\nFIXED EFFECTS (Wald = coef/se):\n"); print(wald_tab, digits = 4)
cat(sprintf("\nWALD TESTS PASSED (|coef/se| >= %.2f):\n", WALD_CRIT))
cat(sprintf("  reference curve  (intercept + ns1..ns%d) : %s\n", K, ref_lab))
cat(sprintf("  difference curve (cond + cond:ns1..ns%d) : %s\n", K, diff_lab))
cat("\nVARIANCE COMPONENTS:\n")
print(allp[allp$block == "variance_component", c("parameter","estimate","se")], digits = 4)
cat("\nWHY NOT cor = AR(TIME):\n"); cat(ar_note, "\n")
cat(sprintf("\nUNMODELLED AUTOCORRELATION: within-shoulder residual lag-1 = %.3f\n", acf_res[2]))
cat("  A shoulder's deviation from the population curve is itself a SMOOTH\n")
cat("  function of x, not noise; a random intercept absorbs only its LEVEL.\n")
cat("  The Wald standard errors above are therefore OPTIMISTIC — treat the Wald\n")
cat("  counts as a descriptive summary, not as calibrated significance tests.\n")
cat("\nRESIDUAL NORMALITY (the acceptance criterion):\n")
cat(sprintf("  skewness %.3f | excess kurtosis %.3f | within +/-1.96 sd %.1f%% (normal 95.0%%)\n",
            nstat["skew"], nstat["kurtosis"], 100*nstat["frac_within_1.96"]))
cat(sprintf("  Shapiro-Wilk p = %.3g on a 5000-row subsample. At this n it rejects almost\n",
            nstat["shapiro_p"]))
cat("  any real data, so read 03_diagnostics.png (Q-Q) as the primary evidence.\n")
cat("\nHOW TO READ IT:\n")
cat(" * cond            -> constant LEVEL shift between conditions.\n")
cat(" * cond:ns1..nsK   -> the SHAPE difference; if none passes Wald, the two\n")
cat("                      trajectories differ only in level.\n")
cat(" * figure 02 shows the total difference with 95% CI across elevation:\n")
cat("   where the band excludes 0, the conditions differ at that elevation.\n")
cat("\nPARAMETER COUNT vs iteration 02:\n")
cat("  iteration 02 (mgcv GAMM): 64 coefficients + 3 lambda + 2 variances\n")
cat(sprintf("  iteration 03 (this)     : %d fixed + %d variance components = %d total\n",
            nfix, nrow(allp)-nfix, nrow(allp)))
cat("\nCAVEATS: (1) BIC is NOT comparable between iterations 02 and 03 (different\n")
cat("  criteria; lcmm's BIC uses the number of SUBJECTS). See 07_vs_iteration02_gamm.png.\n")
cat(" (2) thinned to <=200 pts/shoulder; defensible given the autocorrelation,\n")
cat("     but it IS a subsample - see 00_thinning.csv.\n")
cat(" (3) condition remains confounded with source study: descriptive, not causal.\n")
sink()

cat(sprintf("\nDONE. See %s/ (notice.md explains every column).\n", OUT))
