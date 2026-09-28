# =============================================================================
# GENERALIZED SWEEP v4 — STEP 2: the full model at the SELECTED K
#
# Refits at the K written in K_selection.csv (column selected_K), on ALL
# shoulders, adding the ex-vivo departure when the cell is compared:
#
#   Y ~ (ns1 + ... + nsK) * cond,  random = ~ ns1 + ... + nsK,  D diagonal
#
# K is READ, never re-selected here. Change selected_K in the CSV and re-run: a
# cell whose cached fit sits at a different K is refitted, the rest are
# skipped. in vivo is the reference, so the difference curve is EX VIVO - IN
# VIVO; the script asserts the coefficient came back named `condex vivo`.
#
# A cell that is not compared has no condition term, so it gets ONE pooled
# curve. v3 drew that same curve twice, once per condition colour, which read
# as two conditions with identical bands (v3 00_human_summary.md, OBS5).
#
# Run:  Rscript analyze_all_v4_step2.R      (same V4_ONLY / V4_SHARD / V4_REFIT
#                                            controls as step 1)
# Outputs:
#   cache/step2/<cell>.rds   curves, difference, bf/Vf, coefficients, predictions
# =============================================================================

source("analyze_all_v4_common.R")
suppressMessages(library(lcmm))

combos <- select_combos()
cat(sprintf("STEP 2 — full model at each cell's selected K, %d cells\n", nrow(combos)))
SEL <- read_selection()
if (is.null(SEL)) stop("no ", SEL_CSV, " — run step 1 first")
D <- load_long()
t_start <- Sys.time()

for (i in seq_len(nrow(combos))) {
  jt <- combos$joint[i]; mo <- combos$motion[i]; dof <- combos$dof[i]
  tag <- sprintf("[%2d/%d] %-18s %-30s dof%d", i, nrow(combos), jt, substr(mo, 1, 30), dof)

  s1 <- cell_load(1, jt, mo, dof)
  if (is.null(s1)) { cat(tag, "  no step-1 cache — run step 1 first\n", sep = ""); next }
  K <- selected_K_of(SEL, jt, mo, dof)
  if (is.na(K)) {
    cat(sprintf("%s  SKIP (%s)\n", tag,
                if (nzchar(nz(s1$row$note, ""))) s1$row$note else "no selected_K")); next }
  # A hand-edited K must be one step 1 actually fitted and that converged:
  # otherwise the H0 step 3 needs does not exist and the LRT cannot be formed.
  conv_K <- if (is.null(s1$bic)) integer(0) else s1$bic$K[s1$bic$converged]
  if (!K %in% conv_K)
    stop("cell ", cell_key(jt, mo, dof), ": selected_K = ", K,
         " is not a converged K from step 1 (converged: ",
         if (length(conv_K)) paste(conv_K, collapse = ", ") else "none",
         "). Fix it in ", SEL_CSV)
  s2c <- cell_load(2, jt, mo, dof)
  if (!REFIT && !is.null(s2c)) {
    if (identical(as.integer(s2c$K), as.integer(K))) {
      cat(tag, " cached, skipped\n", sep = ""); next }
    cat(sprintf("%s  selected_K changed %d -> %d, refitting\n", tag, s2c$K, K))
  }

  P <- prepare_cell(D, jt, mo, dof)
  if (is.null(P$st)) { cat(sprintf("%s  SKIP (%s)\n", tag, P$row$note)); next }
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

  # A compared cell gets one curve per condition, each drawn only where that
  # condition is supported and never past the fitted range. A cell that is not
  # compared has no condition term, so it gets ONE curve over all its shoulders.
  groups <- if (compare) levels(droplevels(st$cond)) else "pooled"
  curves <- do.call(rbind, lapply(groups, function(cc) {
    k <- if (cc == "pooled") rep(TRUE, nrow(st)) else st$cond == cc
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
cat("Next:  Rscript analyze_all_v4_step3.R\n")
