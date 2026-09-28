# =============================================================================
# ITERATION 06 — Melanie Prague's inference strategy
#
# Plan: figures/06_natural_spline_hlme/README.md
#
#   1  choose the spline df by BIC on the IN-VIVO data alone, K = 1..10 — in vivo is
#      the reference, so its curve's shape is chosen from its own data
#   1a check each shoulder is individually captured, at EVERY K (max/median RMSE vs K)
#   2  refit at the SAME K adding the ex-vivo departure, IN VIVO as the reference
#   3  likelihood-ratio test, df = K + 1
#   4  no Wald tests anywhere — the difference curve is the reporting object
#
# Run:  Rscript analysis_lrt.R
# Prereq: python3 prepare_monolix_data.py
# =============================================================================

suppressMessages({ library(lcmm); library(splines) })

INPUT   <- "monolix_st_frontal_dof2.csv"
FIG_DIR <- "figures"
S6      <- "06_natural_spline_hlme"
CAP     <- 100          # max points per shoulder
# Step 1 (in vivo only) is cheap — all of K = 1..10 took 147 s — so it sweeps far
# enough to see whether BIC ever turns.
KS      <- 1:20
# The full-data sweep (H0 AND H1 on all 44 shoulders) is NOT cheap and blows up:
# measured 21 s at K = 1, 98 s at K = 9, and 3425 s — 57 min — at K = 10. It feeds
# figures 02-03 and the LRT-by-K table.
#
# THIS MUST REACH THE K THE METHOD CHOOSES. Previously it stopped at 8 while BIC
# chose 10, and the cap silently moved the headline test to K = 8 — which read on
# figure 01 as "retained K = 8" printed next to K = 10's BIC. One method, one K:
# if BIC's minimum moves past max(KS_FULL), raise this rather than cap the answer.
KS_FULL <- 1:10
# K used for the headline test. NA = the BIC winner from step 1. The cap below is
# a guard against an overnight fit, not a choice: if it ever fires, the run is
# reporting a K the method did not select, and the figure now says so loudly.
K_TEST  <- NA
TRIM_TO_OVERLAP <- TRUE # see README: Melanie's "garder toute la plage" is a v2 note
NPROC   <- max(1L, min(8L, parallel::detectCores() - 1L))
FIG_SCALE <- 2.2
# in vivo FIRST => in vivo is the REFERENCE, the departure term is ex vivo
LEVELS  <- c("in vivo", "ex vivo")
COL     <- c("in vivo" = "#1b9e77", "ex vivo" = "#d95f02")
BLUE    <- "#377eb8"

OUT <- file.path(FIG_DIR, S6)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

fig <- function(name, draw, w = 1100, h = 750, res = 130) {
  base <- sub("\\.png$", "", name)
  png(file.path(OUT, paste0(base, ".png")),
      width = round(w*FIG_SCALE), height = round(h*FIG_SCALE), res = round(res*FIG_SCALE))
  op <- par(no.readonly = TRUE); tryCatch(draw(), finally = { par(op); dev.off() })
  if (capabilities("cairo")) cairo_pdf(file.path(OUT, paste0(base, ".pdf")),
                                       width = w/res, height = h/res)
  else pdf(file.path(OUT, paste0(base, ".pdf")), width = w/res, height = h/res)
  op <- par(no.readonly = TRUE); tryCatch(draw(), finally = { par(op); dev.off() })
  message("  saved ", file.path(S6, base), ".png/.pdf")
}

norm_stats <- function(r) {
  r <- r[is.finite(r)]; z <- (r - mean(r))/sd(r)
  sw <- { s <- if (length(r) > 5000) sample(r, 5000) else r; shapiro.test(s)$p.value }
  c(skew = mean(z^3), kurtosis = mean(z^4) - 3, shapiro_p = sw)
}
# concordance R^2 against the 1:1 line — NOT a squared correlation, so a shoulder
# fitted with a constant offset cannot score 0.99 while sitting off the diagonal
r2_conc <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)

# -----------------------------------------------------------------------------
# 1. LOAD — in vivo as the reference level
# -----------------------------------------------------------------------------
if (!file.exists(INPUT)) stop("Missing ", INPUT, " — run prepare_monolix_data.py first.")
d <- read.csv(INPUT, stringsAsFactors = FALSE)
d$cond <- factor(ifelse(d$in_vivo %in% c("True","TRUE","1","yes"), "in vivo", "ex vivo"),
                 levels = LEVELS)
stopifnot(levels(d$cond)[1] == "in vivo")          # the reference must be in vivo
d$ID <- factor(d$ID)
d <- d[is.finite(d$TIME) & is.finite(d$Y), ]
d <- d[order(d$ID, d$TIME), ]

cat("\n================ ITERATION 06: BIC on the reduced model + LRT ================\n")
cat(sprintf("reference condition: %s   (departure term = %s)\n", LEVELS[1], LEVELS[2]))
cat(sprintf("full data: %d obs, %d shoulders (%d in vivo / %d ex vivo rows)\n",
            nrow(d), nlevels(d$ID), sum(d$cond=="in vivo"), sum(d$cond=="ex vivo")))

thin_by_unit <- function(df, cap) do.call(rbind, lapply(split(df, df$ID), function(u)
  if (nrow(u) <= cap) u else u[unique(round(seq(1, nrow(u), length.out = cap))), ]))
dt <- thin_by_unit(d, CAP); dt <- dt[order(dt$ID, dt$TIME), ]

if (TRIM_TO_OVERLAP) {
  rr <- tapply(dt$TIME, droplevels(dt$cond), range)
  ov <- c(max(sapply(rr, `[`, 1)), min(sapply(rr, `[`, 2)))
  dt <- dt[dt$TIME >= ov[1] & dt$TIME <= ov[2], ]
} else ov <- range(dt$TIME)
dt$ID <- droplevels(dt$ID); dt$IDnum <- as.integer(dt$ID)
cat(sprintf("thinned to <=%d/shoulder%s: %d rows, %d shoulders, x in %.1f..%.1f\n",
            CAP, if (TRIM_TO_OVERLAP) " and trimmed to the overlap" else " (full range)",
            nrow(dt), nlevels(dt$ID), ov[1], ov[2]))
cat(sprintf("step 1 will use the IN-VIVO subset only: %d rows, %d shoulders\n",
            sum(dt$cond == "in vivo"), nlevels(droplevels(dt$ID[dt$cond == "in vivo"]))))

# -----------------------------------------------------------------------------
# 2. STEP 1 — BIC over K on the REDUCED model (no condition terms at all)
# -----------------------------------------------------------------------------
build <- function(K) {
  B <- ns(dt$TIME, df = K)                      # quantile knots, boundary at range(x)
  x <- dt; nm <- paste0("ns", seq_len(ncol(B)))
  for (j in seq_along(nm)) x[[nm[j]]] <- B[, j]
  list(df = x, B = B, nm = nm)
}
f_red  <- function(nm) as.formula(paste("Y ~", paste(nm, collapse = " + ")))
f_full <- function(nm) as.formula(paste("Y ~", paste(nm, collapse = " + "), "+ cond +",
                                        paste(paste0("cond:", nm), collapse = " + ")))
f_rand <- function(nm) as.formula(paste("~", paste(nm, collapse = " + ")))

# Step 1 runs on the IN-VIVO data ONLY. in vivo is the reference condition, so the
# shape of the reference trajectory is settled from the reference's own data — the
# ex-vivo shoulders get no say in how many degrees of freedom the common curve has.
# The basis is still built on the full trimmed range, so the K chosen here transfers
# unchanged to the H0/H1 models fitted on everything in steps 2-3.
per_shoulder <- function(m, dat) {
  pr <- m$pred
  do.call(rbind, lapply(split(seq_len(nrow(dat)), droplevels(dat$ID)), function(i) {
    if (!length(i)) return(NULL)
    data.frame(shoulder = as.character(dat$ID[i[1]]),
               condition = as.character(dat$cond[i[1]]), n = length(i),
               rmse = sqrt(mean(pr$resid_ss[i]^2)), max_abs = max(abs(pr$resid_ss[i])),
               r2_marg = r2_conc(dat$Y[i], pr$pred_m[i]),
               r2_ss   = r2_conc(dat$Y[i], pr$pred_ss[i]))
  }))
}

cat("\n---- step 1: BIC over K on the IN-VIVO data only ----\n")
rows <- list(); red <- list()
for (K in KS) {
  b <- build(K)
  din <- b$df[b$df$cond == "in vivo", ]           # <- the reference condition alone
  din$IDnum <- as.integer(droplevels(din$ID))
  t0 <- Sys.time()
  m <- try(hlme(fixed = f_red(b$nm), random = f_rand(b$nm), subject = "IDnum",
                ng = 1, idiag = TRUE, data = din, verbose = FALSE, nproc = NPROC),
           silent = TRUE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok <- !inherits(m, "try-error") && !is.null(m$conv) && m$conv == 1 && all(is.finite(m$best))
  ps <- if (ok) per_shoulder(m, din) else NULL
  cat(sprintf("  K=%2d  %2d fixed  %2d random  %6.1fs  %s\n", K, K + 1, K + 1, secs,
              if (ok) sprintf("BIC %9.1f  loglik %9.1f  rmse max %5.2f med %5.2f",
                              m$BIC, m$loglik, max(ps$rmse), median(ps$rmse))
              else "DID NOT CONVERGE"))
  flush.console()
  rows[[length(rows)+1]] <- data.frame(
    K = K, n_fixed = K + 1, n_random = K + 1, converged = ok,
    loglik = if (ok) m$loglik else NA, BIC = if (ok) m$BIC else NA,
    rmse_max = if (ok) max(ps$rmse) else NA,
    rmse_median = if (ok) median(ps$rmse) else NA,
    rmse_ratio = if (ok) max(ps$rmse)/median(ps$rmse) else NA,
    r2_ss = if (ok) r2_conc(din$Y, m$pred$pred_ss) else NA,
    r2_marg = if (ok) r2_conc(din$Y, m$pred$pred_m) else NA,
    seconds = secs)
  if (ok) red[[as.character(K)]] <- list(m = m, B = b$B, df = din, nm = b$nm, ps = ps)
}
bic <- do.call(rbind, rows); rownames(bic) <- NULL
write.csv(bic, file.path(OUT, "00_bic_by_k.csv"), row.names = FALSE)

okb <- bic[bic$converged, ]
if (!nrow(okb)) stop("no K converged — see 00_bic_by_k.csv")
K_best <- okb$K[which.min(okb$BIC)]
runner <- sort(okb$BIC)[2]
dBIC   <- if (is.na(runner)) NA else runner - min(okb$BIC)
cat(sprintf("\nRETAINED K = %d  (BIC %.1f; next best is %.1f higher)\n",
            K_best, min(okb$BIC), dBIC))

K <- if (!is.na(K_TEST)) K_TEST else min(K_best, max(KS_FULL))
CAPPED  <- K != K_best
CAP_WHY <- if (!is.na(K_TEST)) "K_TEST set explicitly" else "capped at max(KS_FULL)"
if (CAPPED)
  cat(sprintf("\n  NOTE: BIC's winner is K = %d, but the test uses K = %d (%s).\n",
              K_best, K, CAP_WHY))

# Everything reported downstream must describe the K actually USED, not BIC's
# winner — quoting the winner's BIC next to the used K reads as a broken figure.
BIC_K  <- okb$BIC[okb$K == K]
K_note <- if (CAPPED) {
  sprintf("K = %d used (%s); BIC's own winner is K = %d, %.0f lower",
          K, CAP_WHY, K_best, BIC_K - min(okb$BIC))
} else {
  sprintf("K = %d — BIC minimum, next best %.0f higher", K, dBIC)
}
m_in <- red[[as.character(K)]]$m        # step 1's model: IN-VIVO data only
d_in <- red[[as.character(K)]]$df
B    <- red[[as.character(K)]]$B
nm   <- red[[as.character(K)]]$nm

# -----------------------------------------------------------------------------
# 3. STEP 1a — is each in-vivo shoulder individually captured?
#    On step 1's own model, i.e. the in-vivo reference curve, before any ex-vivo
#    data or condition term is involved.
# -----------------------------------------------------------------------------
pr <- m_in$pred
stopifnot(nrow(pr) == nrow(d_in))
stopifnot(max(abs(pr$obs - d_in$Y)) < 1e-8)        # prediction rows align with the data

ind <- red[[as.character(K)]]$ps
ind$ratio_to_median <- ind$rmse / median(ind$rmse)
lg <- log(ind$rmse)
ind$robust_z <- (lg - median(lg)) / (1.4826 * mad(lg, constant = 1))
ind$flagged  <- ind$robust_z > 3
ind <- ind[order(-ind$rmse), ]
write.csv(ind, file.path(OUT, "00_individual_fit.csv"), row.names = FALSE)

R2_MARG <- r2_conc(d_in$Y, pr$pred_m)
R2_SS   <- r2_conc(d_in$Y, pr$pred_ss)
cat("\n---- step 1a: individual fit, IN-VIVO only, K =", K, "----\n")
cat(sprintf("  worst shoulder is %.2fx the typical one (max rmse %.2f / median %.2f)\n",
            max(ind$ratio_to_median), max(ind$rmse), median(ind$rmse)))
cat(sprintf("  flagged (robust z > 3): %d of %d%s\n", sum(ind$flagged), nrow(ind),
            if (any(ind$flagged)) paste0(" — ", paste(ind$shoulder[ind$flagged], collapse = ", ")) else ""))
cat(sprintf("  concordance R2:  marginal %.3f   subject-specific %.3f   (gap %.3f = what the random curves buy)\n",
            R2_MARG, R2_SS, R2_SS - R2_MARG))

# -----------------------------------------------------------------------------
# 4. STEP 2 — full model at the SAME K, ex-vivo departure
# -----------------------------------------------------------------------------
cat("\n---- step 2/3: H0 and H1 on ALL data, at EVERY K ----\n")
# Figures 02-04 show how the answer depends on K, so both models are fitted at
# every df, not only the retained one. That also gives the LRT at every K, which
# is the real question behind "is the verdict an artefact of the df?".
full <- list(); lrows <- list()
for (kk in KS_FULL) {
  if (is.null(red[[as.character(kk)]])) next
  bb <- build(kk); dkk <- bb$df; nmk <- bb$nm; Bk <- bb$B
  t0 <- Sys.time()
  a0 <- try(hlme(fixed = f_red(nmk), random = f_rand(nmk), subject = "IDnum",
                 ng = 1, idiag = TRUE, data = dkk, verbose = FALSE, nproc = NPROC), silent = TRUE)
  a1 <- try(hlme(fixed = f_full(nmk), random = f_rand(nmk), subject = "IDnum",
                 ng = 1, idiag = TRUE, data = dkk, verbose = FALSE, nproc = NPROC), silent = TRUE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok <- !inherits(a0,"try-error") && !inherits(a1,"try-error") &&
        a0$conv == 1 && a1$conv == 1 && all(is.finite(a0$best)) && all(is.finite(a1$best))
  if (!ok) { cat(sprintf("  K=%2d  NOT usable (%.0fs)\n", kk, secs))
    lrows[[length(lrows)+1]] <- data.frame(K = kk, converged = FALSE, loglik_reduced = NA,
      loglik_full = NA, LR = NA, df = kk + 1, p_value = NA, seconds = secs); next }
  nf <- 2*kk + 2
  bf <- a1$best[1:nf]; Vf <- VarCov(a1)[1:nf, 1:nf]
  bnd <- function(X) { f <- as.vector(X %*% bf)
    sdv <- sqrt(pmax(0, rowSums((X %*% Vf) * X)))
    data.frame(fit = f, lo = f - 1.96*sdv, hi = f + 1.96*sdv) }
  cv <- do.call(rbind, lapply(levels(dt$cond), function(cc) {
    r2 <- range(dt$TIME[dt$cond == cc]); xc <- seq(r2[1], r2[2], length.out = 250)
    b2 <- predict(Bk, xc); e <- if (cc == "ex vivo") 1 else 0
    cbind(cond = cc, x = xc, bnd(cbind(1, b2, e, b2*e))) }))
  xdk <- seq(ov[1], ov[2], length.out = 250); bdk <- predict(Bk, xdk)
  df_k <- cbind(x = xdk, bnd(cbind(0, matrix(0, length(xdk), kk), 1, bdk)))
  LRk <- 2*(a1$loglik - a0$loglik); Pk <- pchisq(LRk, df = kk + 1, lower.tail = FALSE)
  cat(sprintf("  K=%2d  %5.0fs  LR %8.2f  df %2d  p %9.2e  %s\n",
              kk, secs, LRk, kk + 1, Pk, if (Pk < 0.05) "reject H0" else "cannot reject"))
  flush.console()
  full[[as.character(kk)]] <- list(m0 = a0, m1 = a1, B = Bk, df = dkk, nm = nmk,
                                   curves = cv, dif = df_k)
  lrows[[length(lrows)+1]] <- data.frame(K = kk, converged = TRUE,
    loglik_reduced = a0$loglik, loglik_full = a1$loglik, LR = LRk, df = kk + 1,
    p_value = Pk, seconds = secs)
}
lrt_all <- do.call(rbind, lrows); rownames(lrt_all) <- NULL
write.csv(lrt_all, file.path(OUT, "00_lrt_by_k.csv"), row.names = FALSE)

stopifnot(!is.null(full[[as.character(K)]]))
m0 <- full[[as.character(K)]]$m0; m1 <- full[[as.character(K)]]$m1
dK <- full[[as.character(K)]]$df
curves <- full[[as.character(K)]]$curves
dif <- full[[as.character(K)]]$dif
xd <- dif$x

nfix <- 2*K + 2
cn <- names(m1$best)[1:nfix]
stopifnot(any(grepl("^condex vivo$", cn)))
if (any(grepl("condin vivo", cn))) stop("reference did not flip — coefficient is condin vivo")
cat("  departure coefficient is '", grep("^cond", cn, value = TRUE)[1],
    "' — reference is in vivo, as intended\n", sep = "")

stopifnot(nrow(m0$pred) == nrow(m1$pred), m0$ns == m1$ns)
LR <- 2 * (m1$loglik - m0$loglik)
DF <- K + 1
P  <- pchisq(LR, df = DF, lower.tail = FALSE)
lrt <- data.frame(K = K, loglik_reduced = m0$loglik, loglik_full = m1$loglik,
                  LR = LR, df = DF, p_value = P,
                  n_rows = nrow(dK), n_subjects = m0$ns,
                  BIC_reduced = m0$BIC, BIC_full = m1$BIC)
write.csv(lrt, file.path(OUT, "00_lrt.csv"), row.names = FALSE)
cat("\n---- step 3: likelihood-ratio test at the retained K ----\n")
cat(sprintf("  LR = 2(%.1f - %.1f) = %.2f   df = %d   p = %.3e   -> %s\n",
            m1$loglik, m0$loglik, LR, DF, P,
            if (P < 0.05) "reject H0: the conditions differ" else "cannot reject H0"))

# coefficients, WITHOUT any Wald column (step 4)
cf <- function(m, lab) {
  V <- VarCov(m)
  data.frame(model = lab, parameter = names(m$best), estimate = as.vector(m$best),
             se = sqrt(diag(V)), row.names = NULL)
}
write.csv(rbind(cf(m0, "reduced"), cf(m1, "full")),
          file.path(OUT, "00_coefficients.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. FIGURES — no knots drawn anywhere
# -----------------------------------------------------------------------------
bfix <- m1$best[1:nfix]; Vfix <- VarCov(m1)[1:nfix, 1:nfix]
Xrow <- function(x, exv) { b <- predict(B, x); cbind(1, b, exv, b*exv) }
band <- function(X) { f <- as.vector(X %*% bfix)
  s <- sqrt(pmax(0, rowSums((X %*% Vfix) * X)))
  data.frame(fit = f, lo = f - 1.96*s, hi = f + 1.96*s) }

curves <- do.call(rbind, lapply(levels(dt$cond), function(cc) {
  r2 <- range(dt$TIME[dt$cond == cc]); xc <- seq(r2[1], r2[2], length.out = 300)
  cbind(cond = cc, x = xc, band(Xrow(xc, if (cc == "ex vivo") 1 else 0)))
}))
xd  <- seq(ov[1], ov[2], length.out = 300)
bd  <- predict(B, xd)
dif <- band(cbind(0, matrix(0, length(xd), K), 1, bd))   # ex vivo MINUS in vivo

fig("01_choosing_k.png", function() {
  op <- par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.4, 1)); on.exit(par(op), add = TRUE)
  o <- bic; ok <- o$converged

  # -- left: BIC, the criterion Melanie's step 1 specifies --------------------
  plot(o$K, o$BIC, type = "n", xlab = "natural-spline df (K)",
       ylab = "BIC (lower is better)", xaxt = "n",
       main = sprintf("BIC on the in-vivo data — K = %d used%s", K,
                      if (CAPPED) sprintf(" (BIC's winner: %d)", K_best) else ""))
  axis(1, at = pretty(KS)); grid(col = "grey90", lty = 1)
  lines(o$K[ok], o$BIC[ok], lwd = 2, col = "grey55")
  points(o$K[ok], o$BIC[ok], pch = 21, bg = "white", col = "grey35", cex = 1.2, lwd = 1.4)
  points(K, o$BIC[o$K == K], pch = 21, bg = BLUE, col = "white", cex = 1.9, lwd = 2)
  # when the test's K is capped below BIC's winner, show BOTH so the blue dot
  # sitting off the minimum reads as deliberate rather than as an error
  if (CAPPED)
    points(K_best, min(okb$BIC), pch = 21, bg = "white", col = BLUE, cex = 1.9, lwd = 2.4)
  if (any(!ok)) { u <- par("usr")
    points(o$K[!ok], rep(u[3] + 0.05*(u[4]-u[3]), sum(!ok)), pch = 4, col = "red", cex = 1.3, lwd = 2)
    text(mean(o$K[!ok]), u[3] + 0.11*(u[4]-u[3]),
         sprintf("no convergence: K = %s", paste(o$K[!ok], collapse = ", ")), col = "red", cex = 0.75) }
  legend("topright", bty = "n", cex = 0.8, pch = if (CAPPED) c(21, 21) else NA,
         pt.bg = if (CAPPED) c(BLUE, "white") else NA,
         col = if (CAPPED) c("white", BLUE) else "black",
         legend = if (CAPPED)
           c(sprintf("K = %d used (%s), BIC %.0f", K, CAP_WHY, BIC_K),
             sprintf("BIC minimum: K = %d, BIC %.0f", K_best, min(okb$BIC)))
         else c(sprintf("retained K = %d, BIC %.0f", K, BIC_K),
                sprintf("next best %.0f higher", dBIC)))

  # -- right: the per-shoulder error, which BIC does not show -----------------
  yr <- range(o$rmse_max[ok], o$rmse_median[ok], na.rm = TRUE)
  plot(o$K, o$rmse_max, type = "n", log = "y", ylim = yr, xaxt = "n",
       xlab = "natural-spline df (K)", ylab = "per-shoulder RMSE (°, log scale)",
       main = "Worst and typical shoulder — in vivo")
  axis(1, at = pretty(KS)); grid(col = "grey90", lty = 1)
  lines(o$K[ok], o$rmse_max[ok], lwd = 2.4, col = "#d95f02")
  points(o$K[ok], o$rmse_max[ok], pch = 21, bg = "#d95f02", col = "white", cex = 1.2, lwd = 1.3)
  lines(o$K[ok], o$rmse_median[ok], lwd = 2.4, col = "#1b9e77", lty = 2)
  points(o$K[ok], o$rmse_median[ok], pch = 21, bg = "#1b9e77", col = "white", cex = 1.1, lwd = 1.3)
  abline(v = K, lty = 3, col = BLUE, lwd = 2)
  legend("topright", bty = "n", cex = 0.8, lwd = 2.4, lty = c(1, 2),
         col = c("#d95f02", "#1b9e77"),
         legend = c("MAX over shoulders", "median over shoulders"))
  legend("bottomleft", bty = "n", cex = 0.75, text.col = "grey30",
         legend = sprintf("at K = %d:  max %.2f°,  median %.2f°  (ratio %.1fx)",
                          K, o$rmse_max[o$K == K], o$rmse_median[o$K == K],
                          o$rmse_ratio[o$K == K]))
}, w = 1500, h = 700, res = 130)

KOK <- KS_FULL[sapply(KS_FULL, function(k) !is.null(full[[as.character(k)]]))]
NK  <- length(KOK)
hi_k <- function(kk) if (kk == K) BLUE else NA      # ring the retained df

# ---- 02: population curves, one ROW per K --------------------------------
# left column is STEP 1 (in vivo only — nothing has been fitted to ex vivo yet),
# right column is STEP 2 (both conditions). Stacked vertically over every K.
fig("02_population_curves.png", function() {
  op <- par(mfcol = c(NK, 2), mar = c(2.4, 3.6, 1.9, 0.6), mgp = c(2.2, 0.6, 0),
            oma = c(2.6, 0, 2.4, 0)); on.exit(par(op), add = TRUE)
  yl <- range(dt$Y)
  for (kk in KOK) {                                   # column 1 — step 1, in vivo
    e <- red[[as.character(kk)]]; mi <- e$m; Bi <- e$B; di <- e$df
    bi <- mi$best[1:(kk+1)]; Vi <- VarCov(mi)[1:(kk+1), 1:(kk+1)]
    rr <- range(di$TIME); xi <- seq(rr[1], rr[2], length.out = 250)
    Xi <- cbind(1, predict(Bi, xi)); fi <- as.vector(Xi %*% bi)
    si <- sqrt(pmax(0, rowSums((Xi %*% Vi) * Xi)))
    plot(di$TIME, di$Y, col = adjustcolor(COL[["in vivo"]], 0.18), pch = 16, cex = 0.25,
         ylim = yl, xlab = "", ylab = sprintf("K=%d", kk), cex.axis = 0.75, cex.lab = 0.95)
    polygon(c(xi, rev(xi)), c(fi-1.96*si, rev(fi+1.96*si)),
            col = adjustcolor(COL[["in vivo"]], 0.3), border = NA)
    lines(xi, fi, col = COL[["in vivo"]], lwd = 2)
    if (kk == K) box(lwd = 2.5, col = BLUE)
    if (kk == KOK[1]) mtext("STEP 1 — in vivo only", side = 3, line = 0.6, font = 2, cex = 0.8)
  }
  for (kk in KOK) {                                   # column 2 — step 2, both
    cv <- full[[as.character(kk)]]$curves
    plot(dt$TIME, dt$Y, col = adjustcolor(COL[as.character(dt$cond)], 0.15), pch = 16,
         cex = 0.25, ylim = yl, xlab = "", ylab = "", cex.axis = 0.75)
    for (cc in unique(cv$cond)) { g <- cv[cv$cond == cc, ]
      polygon(c(g$x, rev(g$x)), c(g$lo, rev(g$hi)), col = adjustcolor(COL[cc], 0.28), border = NA)
      lines(g$x, g$fit, col = COL[cc], lwd = 2) }
    if (kk == K) box(lwd = 2.5, col = BLUE)
    if (kk == KOK[1]) { mtext("STEP 2 — both conditions", side = 3, line = 0.6, font = 2, cex = 0.8)
      legend("topleft", names(COL), col = COL, lwd = 2, bty = "n", cex = 0.7) }
  }
  mtext("thoracohumeral elevation (°)", side = 1, outer = TRUE, line = 1.1, cex = 0.85)
  mtext(sprintf("Population curves at every K — blue box = retained K = %d", K),
        side = 3, outer = TRUE, line = 0.9, font = 2, cex = 0.95)
}, w = 1150, h = 190 * NK + 120, res = 130)

# ---- 03: the difference curve at every K ---------------------------------
fig("03_difference.png", function() {
  op <- par(mfrow = c(NK, 1), mar = c(2.2, 4.4, 1.7, 1), mgp = c(2.6, 0.6, 0),
            oma = c(2.8, 0, 2.4, 0)); on.exit(par(op), add = TRUE)
  yl <- range(0, unlist(lapply(KOK, function(k) { z <- full[[as.character(k)]]$dif
                                                 c(z$lo, z$hi) })))
  for (kk in KOK) {
    z <- full[[as.character(kk)]]$dif; r <- lrt_all[lrt_all$K == kk, ]
    sig <- (z$lo > 0) | (z$hi < 0)
    plot(z$x, z$fit, type = "n", ylim = yl, xlab = "", cex.axis = 0.75,
         ylab = sprintf("K=%d", kk), cex.lab = 0.95)
    polygon(c(z$x, rev(z$x)), c(z$lo, rev(z$hi)), col = adjustcolor(BLUE, 0.25), border = NA)
    abline(h = 0, lty = 2, col = "grey40"); lines(z$x, z$fit, lwd = 2.2, col = BLUE)
    if (any(sig)) points(z$x[sig], rep(yl[1], sum(sig)), pch = 15, cex = 0.3, col = BLUE)
    legend("topleft", bty = "n", cex = 0.7,
           legend = sprintf("mean |Δ| %.2f°   LR %.1f, p %.1e", mean(abs(z$fit)), r$LR, r$p_value))
    if (kk == K) box(lwd = 2.5, col = BLUE)
  }
  mtext("thoracohumeral elevation (°)", side = 1, outer = TRUE, line = 1.2, cex = 0.85)
  mtext("EX VIVO − IN VIVO at every K, 95% CI — blue box = retained K",
        side = 3, outer = TRUE, line = 0.9, font = 2, cex = 0.95)
}, w = 900, h = 175 * NK + 120, res = 130)

# ---- 04: observed vs predicted at every K (step 1a, in vivo) --------------
fig("04_observed_vs_predicted.png", function() {
  op <- par(mfcol = c(NK, 2), mar = c(2.3, 3.4, 1.8, 0.6), mgp = c(2.1, 0.6, 0),
            oma = c(2.8, 0, 2.4, 0)); on.exit(par(op), add = TRUE)
  lim <- range(d_in$Y)
  for (w in c("marg", "ss")) for (kk in KOK) {
    e <- red[[as.character(kk)]]; prk <- e$m$pred; dk2 <- e$df
    yh <- if (w == "marg") prk$pred_m else prk$pred_ss
    plot(yh, dk2$Y, pch = 16, cex = 0.22, xlim = lim, ylim = lim, asp = 1,
         col = adjustcolor(COL[["in vivo"]], 0.3), cex.axis = 0.72,
         xlab = "", ylab = if (w == "marg") sprintf("K=%d", kk) else "", cex.lab = 0.95)
    abline(0, 1, col = "red", lwd = 1.6)
    legend("topleft", bty = "n", cex = 0.68,
           legend = sprintf("R² = %.3f", r2_conc(dk2$Y, yh)))
    if (kk == K) box(lwd = 2.5, col = BLUE)
    if (kk == KOK[1]) mtext(if (w == "marg") "marginal — population curve only"
                            else "subject-specific — + BLUP",
                            side = 3, line = 0.6, font = 2, cex = 0.8)
  }
  mtext(expression(hat(y)[ij]), side = 1, outer = TRUE, line = 1.2, cex = 0.9)
  mtext(sprintf("Step 1a — observed vs predicted, IN VIVO, every K (%d shoulders)", nrow(ind)),
        side = 3, outer = TRUE, line = 0.9, font = 2, cex = 0.95)
}, w = 1000, h = 175 * NK + 120, res = 130)

# ---- 05: per-shoulder fit, with EVERY K overlaid on the worst shoulders ----
fig("05_individual_fit.png", function() {
  op <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
  ramp <- colorRampPalette(c("#c6dbef", "#08306b"))(length(KOK))
  o <- ind[order(ind$rmse), ]
  barplot(o$rmse, col = ifelse(o$flagged, "red", "grey70"), border = NA,
          names.arg = rep("", nrow(o)), ylab = "per-shoulder RMSE (°)",
          xlab = "in-vivo shoulders, sorted",
          main = sprintf("At K=%d: worst / median = %.2fx", K, max(ind$ratio_to_median)))
  abline(h = median(ind$rmse), lty = 2, col = "grey30")
  thr <- exp(median(log(ind$rmse)) + 3*1.4826*mad(log(ind$rmse), constant = 1))
  abline(h = thr, lty = 3, col = "red")
  legend("topleft", bty = "n", cex = 0.75, lty = c(2,3), col = c("grey30","red"),
         legend = c("median", "robust flag (z = 3)"))
  for (sh in head(ind$shoulder, 5)) {
    i2 <- which(as.character(d_in$ID) == sh)
    plot(d_in$TIME[i2], d_in$Y[i2], pch = 16, cex = 0.5, col = adjustcolor("grey40", 0.6),
         xlab = "elevation (°)", ylab = "angle (°)", cex.main = 0.82,
         main = sprintf("%s%s\nrmse %.2f° at K=%d", sh,
                        if (ind$flagged[ind$shoulder == sh]) "  [FLAGGED]" else "",
                        ind$rmse[ind$shoulder == sh], K))
    # every K overlaid, pale -> dark, so the retained curve is seen in context
    for (q in seq_along(KOK)) {
      kk <- KOK[q]; e <- red[[as.character(kk)]]
      ii <- which(as.character(e$df$ID) == sh); if (!length(ii)) next
      o2 <- order(e$df$TIME[ii])
      lines(e$df$TIME[ii][o2], e$m$pred$pred_ss[ii][o2],
            col = ramp[q], lwd = if (kk == K) 2.6 else 1)
    }
    if (sh == head(ind$shoulder, 1))
      legend("bottomright", bty = "n", cex = 0.62, lwd = c(1, 2.6),
             col = c(ramp[1], ramp[length(ramp)]),
             legend = c(sprintf("K=%d", min(KOK)), sprintf("K=%d (retained)", K)))
  }
}, w = 1500, h = 900, res = 130)

res <- m1$pred$resid_ss; z <- (res - mean(res))/sd(res); ns_st <- norm_stats(res)
fig("06_diagnostics.png", function() {
  op <- par(mfrow = c(2,2), mar = c(4,4,3,1)); on.exit(par(op), add = TRUE)
  plot(m1$pred$pred_ss, res, pch = 16, cex = 0.4, col = adjustcolor("grey20", 0.3),
       xlab = "fitted", ylab = "residual", main = "Residuals vs fitted"); abline(h = 0, col = "red")
  qqnorm(z, pch = 16, cex = 0.4, col = adjustcolor("grey20", 0.3),
         main = "Normal Q-Q"); qqline(z, col = "red", lwd = 2)
  hist(z, breaks = 60, freq = FALSE, border = NA, col = "grey80", xlab = "standardised residual",
       main = sprintf("skew %.2f | kurtosis %.2f", ns_st["skew"], ns_st["kurtosis"]))
  curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
  plot(dK$TIME, res, pch = 16, cex = 0.4, col = adjustcolor(COL[as.character(dK$cond)], 0.35),
       xlab = "elevation (°)", ylab = "residual", main = "Residuals vs x"); abline(h = 0, col = "red")
})

# -----------------------------------------------------------------------------
# 7. DIGEST
# -----------------------------------------------------------------------------
sink(file.path(OUT, "00_results_summary.txt"))
cat("Iteration 06 — BIC on the reduced model, then a likelihood-ratio test\n")
cat("=====================================================================\n\n")
cat(sprintf("Reference condition: %s.  The departure term is %s.\n", LEVELS[1], LEVELS[2]))
cat(sprintf("The difference curve is EX VIVO minus IN VIVO — mirrored vs iterations 03-05.\n\n"))
cat(sprintf("Data: %d obs -> %d after thinning to <=%d/shoulder%s; %d shoulders; x in %.1f..%.1f\n\n",
            nrow(d), nrow(dK), CAP,
            if (TRIM_TO_OVERLAP) " and trimming to the overlap" else " (full range kept)",
            nlevels(dK$ID), ov[1], ov[2]))
cat("STEP 1 — BIC over the reduced model (no condition terms):\n"); print(bic, digits = 6)
cat(sprintf("\n  %s\n", K_note))
cat("\nSTEP 1a — individual fit at that K:\n")
cat(sprintf("  worst shoulder / median RMSE = %.2f\n", max(ind$ratio_to_median)))
cat(sprintf("  flagged (robust z > 3): %d of %d\n", sum(ind$flagged), nrow(ind)))
cat(sprintf("  concordance R2: marginal %.3f, subject-specific %.3f (gap %.3f)\n",
            R2_MARG, R2_SS, R2_SS - R2_MARG))
cat("\n  The gap is what the random curves buy: the marginal prediction uses the\n")
cat("  population curve alone, the subject-specific one adds each shoulder's BLUP.\n")
cat("\nSTEP 3 — likelihood-ratio test (H0 = no ex-vivo departure):\n")
cat(sprintf("  loglik reduced %.2f | full %.2f\n", m0$loglik, m1$loglik))
cat(sprintf("  LR = %.2f, df = %d, p = %.4e  ->  %s\n", LR, DF, P,
            if (P < 0.05) "reject H0" else "cannot reject H0"))
cat("\nSTEP 4 — no Wald tests. The reporting object is 03_difference.png.\n")
cat(sprintf("  mean |ex vivo - in vivo| = %.2f deg over %.0f..%.0f deg\n",
            mean(abs(dif$fit)), ov[1], ov[2]))
cat("\nRESIDUALS:\n")
cat(sprintf("  skew %.3f | excess kurtosis %.3f | Shapiro p %.3g (subsampled)\n",
            ns_st["skew"], ns_st["kurtosis"], ns_st["shapiro_p"]))
cat("\nCAVEAT: condition is perfectly confounded with source study - descriptive, not causal.\n")
cat("No multiple-comparison adjustment here: one cell, one test. Bonferroni across the\n")
cat("plate belongs to the generalized rewrite.\n")
sink()

cat(sprintf("\nDONE. See %s/\n", OUT))
