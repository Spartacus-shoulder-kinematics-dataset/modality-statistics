# =============================================================================
# GENERALIZED SWEEP v3 — STEP 2: the full model at the SAME K
#
# Refits at the K step 1 chose, on ALL shoulders, adding the ex-vivo departure:
#
#   Y ~ (ns1 + ... + nsK) * cond,  random = ~ ns1 + ... + nsK,  D diagonal
#
# K is INHERITED from step 1 and never re-selected here. in vivo is the
# reference, so the departure term is ex vivo and the difference curve is
# EX VIVO - IN VIVO — mirrored against v2, which had ex vivo as the reference.
# The script asserts the coefficient came back named `condex vivo`: if the
# reference ever failed to flip, every sign downstream would be wrong.
#
# Single-condition cells stop at their reduced curve. There is nothing to
# compare them against, so they get no H1 and no difference.
#
# Run:  Rscript analyze_all_v3_step2.R      (same V3_ONLY / V3_SHARD / V3_REFIT
#                                            controls as step 1)
# Outputs:
#   cache/step2/<cell>.rds   curves, difference, bf/Vf, coefficients, predictions
# =============================================================================

source("analyze_all_v3_common.R")
suppressMessages(library(lcmm))

combos <- select_combos()
cat(sprintf("STEP 2 — full model at each cell's own K, %d cells\n", nrow(combos)))
D <- load_long()
t_start <- Sys.time()

for (i in seq_len(nrow(combos))) {
  jt <- combos$joint[i]; mo <- combos$motion[i]; dof <- combos$dof[i]
  tag <- sprintf("[%2d/%d] %-18s %-30s dof%d", i, nrow(combos), jt, substr(mo, 1, 30), dof)

  s1 <- cell_load(1, jt, mo, dof)
  if (is.null(s1)) { cat(tag, "  no step-1 cache — run step 1 first\n", sep = ""); next }
  if (is.na(s1$K_star)) { cat(sprintf("%s  SKIP (%s)\n", tag, s1$row$note)); next }
  if (!REFIT && !is.null(cell_load(2, jt, mo, dof))) {
    cat(tag, " cached, skipped\n", sep = ""); next }

  P <- prepare_cell(D, jt, mo, dof)
  if (is.null(P$st)) { cat(sprintf("%s  SKIP (%s)\n", tag, P$row$note)); next }
  K <- s1$K_star
  b <- cell_basis(P$st, K, P$ov)
  if (is.null(b)) { cat(sprintf("%s  SKIP (basis failed at K=%d)\n", tag, K)); next }
  st <- b$df; Kc <- b$Kc; compare <- P$compare; ov <- P$ov

  t0 <- Sys.time()
  m1 <- try(hlme(fixed = if (compare) f_full(b$nm) else f_red(b$nm),
                 random = f_rand(b$nm), subject = "IDnum", ng = 1, idiag = TRUE,
                 data = st, verbose = FALSE, nproc = NPROC), silent = TRUE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!fit_ok(m1)) {
    cat(sprintf("%s  DID NOT CONVERGE at K=%d (%.0fs)\n", tag, K, secs))
    P$row$note <- sprintf("H1 did not converge at K=%d", K)
    cell_save(2, jt, mo, dof, list(row = P$row, K = K, converged = FALSE)); next
  }

  nf <- if (compare) 2*Kc + 2 else Kc + 1
  cn <- names(m1$best)[1:nf]
  if (compare) {
    # The reference must be in vivo. If it ever flips, the difference curve, its
    # sign and every conclusion drawn from it invert silently — so stop here.
    if (any(grepl("condin vivo", cn)))
      stop("cell ", cell_key(jt, mo, dof), ": reference did not flip — coefficient is condin vivo")
    # And the departure block must sit exactly where the design matrix below
    # assumes: [intercept, ns1..nsK, condex vivo, ns1:condex vivo .. nsK:condex vivo].
    # Xrow() and Xd index by POSITION, so a name in the right set but the wrong
    # slot would silently contrast the wrong thing.
    want <- c("intercept", b$nm, "condex vivo", paste0(b$nm, ":condex vivo"))
    if (!identical(cn, want))
      stop("cell ", cell_key(jt, mo, dof), ": unexpected fixed-effect layout.\n  expected: ",
           paste(want, collapse = ", "), "\n  got:      ", paste(cn, collapse = ", "))
  }

  bf <- m1$best[1:nf]; Vf <- VarCov(m1)[1:nf, 1:nf, drop = FALSE]
  band <- function(X) { f <- as.vector(X %*% bf)
    s <- sqrt(pmax(0, rowSums((X %*% Vf) * X)))
    data.frame(fit = f, lo = f - 1.96*s, hi = f + 1.96*s) }
  Xrow <- function(x, exv) { z <- predict(b$B, x)
    if (compare) cbind(1, z, exv, z*exv) else cbind(1, z) }

  # Each condition's curve is drawn only where that condition is supported, and
  # never past the fitted range.
  curves <- do.call(rbind, lapply(levels(droplevels(st$cond)), function(cc) {
    k <- st$cond == cc
    r2 <- support_range(st$TIME[k], st$ID[k])
    r2 <- c(max(r2[1], ov[1]), min(r2[2], ov[2]))
    if (!(r2[2] > r2[1])) return(NULL)
    xc <- seq(r2[1], r2[2], length.out = 250)
    cbind(cond = cc, x = xc, band(Xrow(xc, if (cc == "ex vivo") 1 else 0))) }))

  dif <- NULL; dstat <- list()
  if (compare) {
    xd <- seq(ov[1], ov[2], length.out = 250)
    bd <- predict(b$B, xd)
    Xd <- cbind(0, matrix(0, length(xd), Kc), 1, bd)     # gamma_0 + gamma_k B_k(x)
    dif <- cbind(x = xd, band(Xd))
    # 95% CI of the MEAN difference: the average design row is itself a contrast,
    # so its variance is xbar' V xbar — exact, not an average of pointwise SEs
    xbar <- colMeans(Xd)
    mse  <- sqrt(max(0, as.numeric(t(xbar) %*% Vf %*% xbar)))
    dstat <- list(diff_mean_signed = mean(dif$fit), diff_mean_abs = mean(abs(dif$fit)),
                  diff_max_abs = max(abs(dif$fit)),
                  diff_lo_curve = min(dif$fit), diff_hi_curve = max(dif$fit),
                  # fraction of the fitted range where the 95% band excludes zero,
                  # exported so it is data rather than a number read off a picture
                  sig_frac_x = mean((dif$lo > 0) | (dif$hi < 0)),
                  diff_mean_lo = mean(dif$fit) - 1.96*mse,
                  diff_mean_hi = mean(dif$fit) + 1.96*mse)
  }

  rr <- m1$pred$resid_ss; ns_ <- norm_stats(rr)
  coefs <- data.frame(model = "full", parameter = names(m1$best),
                      estimate = as.vector(m1$best), se = sqrt(diag(VarCov(m1))),
                      row.names = NULL)     # estimate and se ONLY — no Wald column

  cell_save(2, jt, mo, dof, c(list(
    row = P$row, K = K, converged = TRUE, compare = compare,
    curves = curves, dif = dif, bf = bf, Vf = Vf, knots = b$knots, Kc = Kc,
    full_rng = P$full_rng, ov = ov, coefs = coefs,
    loglik_full = m1$loglik, BIC_full = m1$BIC, npm_full = length(m1$best),
    n_rows = nrow(st), n_subjects = m1$ns,
    pred_full = pred_rows(m1, st),
    raw = st[, c("ID", "cond", "TIME", "Y")],
    resid_sd = sd(rr), skew = ns_[["skew"]], kurtosis = ns_[["kurtosis"]],
    shapiro_p = ns_[["shapiro_p"]], seconds = secs), dstat))

  cat(sprintf("%s  %-7s K=%2d  %5.0fs%s\n", tag, P$row$mode, K, secs,
              if (compare) sprintf("  |Δ| %.2f° (max %.2f°)  sig over %.0f%% of range",
                                   dstat$diff_mean_abs, dstat$diff_max_abs,
                                   100*dstat$sig_frac_x) else ""))
  flush.console()
}

cat(sprintf("\nSTEP 2 done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("Next:  Rscript analyze_all_v3_step3.R\n")
