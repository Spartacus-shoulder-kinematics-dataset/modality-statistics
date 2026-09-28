# =============================================================================
# WHAT IS THE SPLINE df ACTUALLY CHOOSING?
#
# Iteration 06 selects K by BIC on the reduced in-vivo model, over K = 1..20,
# with the random part MIRRORING K.  That rule gives K = 10, and with the grid
# extended to 20 it is a genuine interior minimum rather than the edge of the
# grid.  This script does not propose a rival rule -- BIC decides.  It answers a
# different question: what is the number 10 a statement ABOUT?
#
# Two checks, in order:
#
#   1. Does the answer depend on what `n` means?  lcmm penalises by the number
#      of SUBJECTS (31 in-vivo shoulders), so an extra spline df costs only
#      2*log(31) = 6.87 BIC.  Recomputed on the 2,277 in-vivo OBSERVATIONS the
#      penalty is 2*log(2277) = 15.46, more than twice as steep.  -> figure 08.
#
#   2. The random part mirrors K, so every added df buys flexibility twice --
#      once in the shared curve, once again in each shoulder's own.  Hold the
#      random part FIXED at R+1 random effects and sweep the fixed df alone: if
#      the steep BIC curve survives, K really is about the population curve; if
#      it collapses, K is mostly about the random dimension.  -> figure 09.
#
# The second check is the informative one.  Two of its cells are, by
# construction, the same model as the mirrored sweep (R = 1 at K = 1, R = 3 at
# K = 3) and are asserted against it, so a wiring error cannot pass silently.
#
# REQUIRES analysis_lrt.R to have run first -- it reads that script's
# 00_bic_by_k.csv rather than refitting the mirrored sweep.
#
# Run:  Rscript analysis_random_dim.R          (~13 min)
# =============================================================================

suppressMessages({ library(lcmm); library(splines) })

INPUT <- "monolix_st_frontal_dof2.csv"
S6    <- "06_natural_spline_hlme"
OUT   <- file.path("figures", S6)
CAP   <- 100                 # max points per shoulder      as analysis_lrt.R
KS    <- 1:20                # fixed spline df to sweep     as analysis_lrt.R
RS    <- c(1L, 3L)           # random spline cols -> R+1 random effects
TRIM_TO_OVERLAP <- TRUE
LEVELS    <- c("in vivo", "ex vivo")
NPROC     <- max(1L, min(8L, parallel::detectCores() - 1L))
FIG_SCALE <- 2.2
BLUE <- "#377eb8"; CM <- "#7570b3"; C3 <- "#1b9e77"; C1 <- "#d95f02"; ORANGE <- "#d95f02"

fig <- function(name, draw, w = 1500, h = 640, res = 130) {
  png(file.path(OUT, paste0(name, ".png")), width = round(w*FIG_SCALE),
      height = round(h*FIG_SCALE), res = round(res*FIG_SCALE))
  op <- par(no.readonly = TRUE); tryCatch(draw(), finally = { par(op); dev.off() })
  if (capabilities("cairo")) cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w/res, height = h/res)
  else pdf(file.path(OUT, paste0(name, ".pdf")), width = w/res, height = h/res)
  op <- par(no.readonly = TRUE); tryCatch(draw(), finally = { par(op); dev.off() })
  message("  saved ", file.path(S6, name), ".png/.pdf")
}
r2_conc <- function(y, yh) 1 - sum((y - yh)^2) / sum((y - mean(y))^2)

# -----------------------------------------------------------------------------
# 1. LOAD — identical to analysis_lrt.R, so every K is comparable to its output
# -----------------------------------------------------------------------------
if (!file.exists(INPUT)) stop("Missing ", INPUT, " — run prepare_monolix_data.py first.")
d <- read.csv(INPUT, stringsAsFactors = FALSE)
d$cond <- factor(ifelse(d$in_vivo %in% c("True","TRUE","1","yes"), "in vivo", "ex vivo"),
                 levels = LEVELS)
d$ID <- factor(d$ID)
d <- d[is.finite(d$TIME) & is.finite(d$Y), ]
d <- d[order(d$ID, d$TIME), ]

thin_by_unit <- function(df, cap) do.call(rbind, lapply(split(df, df$ID), function(u)
  if (nrow(u) <= cap) u else u[unique(round(seq(1, nrow(u), length.out = cap))), ]))
dt <- thin_by_unit(d, CAP); dt <- dt[order(dt$ID, dt$TIME), ]
if (TRIM_TO_OVERLAP) {
  rr <- tapply(dt$TIME, droplevels(dt$cond), range)
  ov <- c(max(sapply(rr, `[`, 1)), min(sapply(rr, `[`, 2)))
  dt <- dt[dt$TIME >= ov[1] & dt$TIME <= ov[2], ]
} else ov <- range(dt$TIME)
dt$ID <- droplevels(dt$ID)

N_OBS <- sum(dt$cond == "in vivo")
N_SUB <- nlevels(droplevels(dt$ID[dt$cond == "in vivo"]))
cat(sprintf("\nstep-1 in-vivo data: %d obs, %d shoulders, x in %.1f..%.1f\n",
            N_OBS, N_SUB, ov[1], ov[2]))
cat(sprintf("BIC penalty per extra spline df: %.2f (lcmm, %d subjects) vs %.2f (%d obs)\n\n",
            2*log(N_SUB), N_SUB, 2*log(N_OBS), N_OBS))

# -----------------------------------------------------------------------------
# 2. CHECK 1 — the mirrored sweep under an observation-based penalty
# -----------------------------------------------------------------------------
f_bic <- file.path(OUT, "00_bic_by_k.csv")
if (!file.exists(f_bic)) stop("Missing ", f_bic, " — run analysis_lrt.R first.")
mir <- read.csv(f_bic)
# recover the parameter count from lcmm's own identity, rather than assuming it:
#   BIC = npm*log(ns) - 2*loglik
mir$npm      <- round((mir$BIC + 2*mir$loglik) / log(N_SUB))
stopifnot(all(mir$npm == 2*mir$K + 3))   # (K+1) fixed + (K+1) variances + sigma
mir$deviance <- -2*mir$loglik
mir$BIC_obs  <- mir$deviance + mir$npm*log(N_OBS)
write.csv(mir, file.path(OUT, "08_bic_by_k_extended.csv"), row.names = FALSE)
message("  saved ", file.path(S6, "08_bic_by_k_extended.csv"))

o <- mir[mir$converged, ]
turn <- function(v, K, label) {
  kb <- K[which.min(v)]
  if (kb == max(K)) cat(sprintf("  %-26s minimum STILL at the last K (%d) — has not turned\n", label, kb))
  else cat(sprintf("  %-26s minimum at K = %d  (BIC rises from K = %d on)\n", label, kb, kb + 1))
  kb
}
K_bic     <- turn(o$BIC,     o$K, "BIC (lcmm, subjects)")
K_bic_obs <- turn(o$BIC_obs, o$K, "BIC (observations)")
if (K_bic != K_bic_obs)
  cat("  NOTE: the two penalties disagree — say so in the README.\n")

fig("08_bic_two_penalties", function() {
  op <- par(mfrow = c(1, 2), mar = c(4.4, 4.8, 3.6, 1)); on.exit(par(op), add = TRUE)
  plot(o$K, o$BIC, type = "n", ylim = range(o$BIC, o$BIC_obs), xaxt = "n",
       xlab = "natural-spline df (K)", ylab = "BIC (lower is better)",
       main = "BIC turns — and where does not depend\non what n means")
  axis(1, at = seq(2, 20, 2)); grid(col = "grey90", lty = 1)
  lines(o$K, o$BIC, lwd = 2.4, col = "grey40")
  points(o$K, o$BIC, pch = 21, bg = "white", col = "grey40", cex = 1.0, lwd = 1.3)
  lines(o$K, o$BIC_obs, lwd = 2.4, col = ORANGE, lty = 2)
  points(o$K, o$BIC_obs, pch = 21, bg = "white", col = ORANGE, cex = 1.0, lwd = 1.3)
  abline(v = K_bic, lty = 3, col = "grey40", lwd = 2)
  points(K_bic, o$BIC[o$K == K_bic], pch = 21, bg = "grey40", col = "grey20", cex = 1.6, lwd = 1.4)
  legend("top", bty = "n", cex = 0.78, lwd = 2.4, lty = c(1, 2), col = c("grey40", ORANGE),
         legend = c(sprintf("BIC on %d subjects (lcmm): K = %d", N_SUB, K_bic),
                    sprintf("BIC on %d observations: K = %d", N_OBS, K_bic_obs)))

  plot(o$K, o$rmse_max, type = "n", ylim = c(0, max(o$rmse_max)), xaxt = "n",
       xlab = "natural-spline df (K)", ylab = "per-shoulder RMSE (deg)",
       main = "What K = 10 buys: the worst shoulder stops\nimproving while the typical one does not")
  axis(1, at = seq(2, 20, 2)); grid(col = "grey90", lty = 1)
  lines(o$K, o$rmse_max, lwd = 2.4, col = ORANGE)
  points(o$K, o$rmse_max, pch = 21, bg = "white", col = ORANGE, cex = 1.0, lwd = 1.3)
  lines(o$K, o$rmse_median, lwd = 2.4, col = BLUE)
  points(o$K, o$rmse_median, pch = 21, bg = "white", col = BLUE, cex = 1.0, lwd = 1.3)
  abline(v = K_bic, lty = 3, col = "grey40", lwd = 2)
  legend("topright", bty = "n", cex = 0.78, lwd = 2.4, col = c(ORANGE, BLUE),
         legend = c("max over shoulders", "median over shoulders"))
})

# -----------------------------------------------------------------------------
# 3. CHECK 2 — decouple the random part from K
# -----------------------------------------------------------------------------
f_dec <- file.path(OUT, "09_decoupled_random.csv")
want  <- do.call(rbind, lapply(RS, function(R) data.frame(R = R, K = KS[KS >= R])))
# Redrawing a legend should not cost 13 minutes of refitting, so the sweep is
# cached: it re-runs only when the cache is missing or does not cover the grid.
# Set REFIT <- TRUE (or delete the CSV) to force it.
REFIT <- FALSE
have  <- if (!REFIT && file.exists(f_dec)) read.csv(f_dec) else NULL
if (!is.null(have) && !all(paste(want$R, want$K) %in% paste(have$R, have$K))) {
  cat("  cache does not cover the requested grid — refitting\n"); have <- NULL
}
if (!is.null(have)) {
  cat("\n---- reusing cached sweep from ", f_dec, " ----\n", sep = "")
  dec <- have
} else {
cat("\n---- sweeping fixed K with the random part HELD FIXED ----\n")
rows <- list()
for (R in RS) for (K in KS) {
  if (R > K) next
  BF <- ns(dt$TIME, df = K); BR <- ns(dt$TIME, df = R)
  x <- dt
  nmF <- paste0("f", seq_len(ncol(BF))); for (j in seq_along(nmF)) x[[nmF[j]]] <- BF[, j]
  nmR <- paste0("r", seq_len(ncol(BR))); for (j in seq_along(nmR)) x[[nmR[j]]] <- BR[, j]
  din <- x[x$cond == "in vivo", ]; din$IDnum <- as.integer(droplevels(din$ID))
  t0 <- Sys.time()
  m <- try(hlme(fixed  = as.formula(paste("Y ~", paste(nmF, collapse = " + "))),
                random = as.formula(paste("~", paste(nmR, collapse = " + "))),
                subject = "IDnum", ng = 1, idiag = TRUE, data = din,
                verbose = FALSE, nproc = NPROC), silent = TRUE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok <- !inherits(m, "try-error") && !is.null(m$conv) && m$conv == 1 && all(is.finite(m$best))
  rs <- if (ok) sapply(split(seq_len(nrow(din)), droplevels(din$ID)),
                       function(i) sqrt(mean(m$pred$resid_ss[i]^2))) else NA
  cat(sprintf("  R=%d K=%2d  %6.1fs  %s\n", R, K, secs,
      if (ok) sprintf("npm %2d  BIC %8.1f  -2LL %8.1f  rmse max %5.2f med %5.2f",
                      length(m$best), m$BIC, -2*m$loglik, max(rs), median(rs))
      else "DID NOT CONVERGE")); flush.console()
  rows[[length(rows)+1]] <- data.frame(R = R, K = K, converged = ok,
    npm = if (ok) length(m$best) else NA, loglik = if (ok) m$loglik else NA,
    BIC = if (ok) m$BIC else NA, rmse_max = if (ok) max(rs) else NA,
    rmse_median = if (ok) median(rs) else NA,
    r2_ss = if (ok) r2_conc(din$Y, m$pred$pred_ss) else NA, seconds = secs)
}
dec <- do.call(rbind, rows); rownames(dec) <- NULL
dec$deviance <- -2*dec$loglik
write.csv(dec, f_dec, row.names = FALSE)
message("  saved ", file.path(S6, "09_decoupled_random.csv"))
}
if (is.null(dec$deviance)) dec$deviance <- -2*dec$loglik

r1 <- dec[dec$R == 1 & dec$converged, ]; r3 <- dec[dec$R == 3 & dec$converged, ]

# R=1,K=1 and R=3,K=3 ARE cells of the mirrored sweep — assert, do not assume
chk <- function(sub, k, label) {
  a <- sub$BIC[sub$K == k]; b <- mir$BIC[mir$K == k]
  cat(sprintf("  check %-22s BIC %.2f vs mirrored %.2f  (diff %.3f)\n", label, a, b, a - b))
  if (abs(a - b) > 1) stop("decoupled sweep does not reproduce the mirrored model at K = ", k)
}
cat("\n---- consistency against the mirrored sweep ----\n")
chk(r1, 1, "R=1, K=1"); chk(r3, 3, "R=3, K=3")

# -----------------------------------------------------------------------------
# 4. FIGURE 09 — the punchline
# -----------------------------------------------------------------------------
fig("09_what_the_df_buys", function() {
  op <- par(mfrow = c(1, 2), mar = c(4.4, 4.8, 4.0, 1)); on.exit(par(op), add = TRUE)

  plot(c(1, 23), range(c(mir$deviance, r1$deviance, r3$deviance)), type = "n", log = "y",
       xaxt = "n", xlab = "FIXED natural-spline df (K)", ylab = "deviance, -2 logLik (log scale)",
       main = "Fixed df buy almost nothing.\nThe RANDOM dimension is doing all the work.")
  axis(1, at = seq(2, 20, 2)); grid(col = "grey90", lty = 1)
  for (s in list(list(r1, C1), list(r3, C3), list(mir, CM))) {
    lines(s[[1]]$K, s[[1]]$deviance, lwd = 2.6, col = s[[2]])
    points(s[[1]]$K, s[[1]]$deviance, pch = 21, bg = "white", col = s[[2]], cex = 0.95, lwd = 1.4)
  }
  # placed in the empty band between the flat orange and green curves — a corner
  # legend here would sit on top of the purple curve's minimum, the whole point
  legend(x = 4.2, y = 8300, cex = 0.78, lwd = 2.6, col = c(CM, C3, C1), bty = "n",
         legend = c("random part MIRRORS K  (K+1 random effects)",
                    "random part FIXED at 4 random effects",
                    "random part FIXED at 2 random effects"))
  for (s in list(list(r1, C1), list(r3, C3), list(mir, CM))) {
    v <- s[[1]]$deviance[nrow(s[[1]])]
    text(20.6, v, sprintf("%.0f", v), pos = 4, col = s[[2]], cex = 0.8)
  }

  plot(c(1, 20), range(c(mir$BIC, r1$BIC, r3$BIC)), type = "n", log = "y", xaxt = "n",
       xlab = "FIXED natural-spline df (K)", ylab = "BIC (log scale, lower is better)",
       main = "So BIC's 'choice of K' is really a choice\nof how many random effects to allow")
  axis(1, at = seq(2, 20, 2)); grid(col = "grey90", lty = 1)
  for (s in list(list(r1, C1), list(r3, C3), list(mir, CM))) {
    lines(s[[1]]$K, s[[1]]$BIC, lwd = 2.6, col = s[[2]])
    points(s[[1]]$K, s[[1]]$BIC, pch = 21, bg = "white", col = s[[2]], cex = 0.95, lwd = 1.4)
    kb <- s[[1]]$K[which.min(s[[1]]$BIC)]
    points(kb, min(s[[1]]$BIC), pch = 21, bg = s[[2]], col = "white", cex = 1.8, lwd = 2)
  }
  legend(x = 4.2, y = 8300, cex = 0.78, pch = 21, pt.lwd = 2, pt.cex = 1.5, bty = "n",
         pt.bg = c(CM, C3, C1), col = "white",
         legend = c(sprintf("mirrors K:  BIC min at K = %d", mir$K[which.min(mir$BIC)]),
                    sprintf("4 random effects:  K = %d (basin flat from K = 8)", r3$K[which.min(r3$BIC)]),
                    sprintf("2 random effects:  K = %d", r1$K[which.min(r1$BIC)])))
})

# -----------------------------------------------------------------------------
# 5. THE NUMBERS THE README QUOTES
# -----------------------------------------------------------------------------
g_mir <- mir$deviance[1] - mir$deviance[nrow(mir)]
g_r1  <- r1$deviance[1]  - r1$deviance[nrow(r1)]
g_r3  <- r3$deviance[1]  - r3$deviance[nrow(r3)]
cat("\n---- what an extra fixed df buys, by random dimension ----\n")
cat(sprintf("  mirrors K        : deviance %8.1f -> %8.1f   gain %8.1f\n",
            mir$deviance[1], mir$deviance[nrow(mir)], g_mir))
cat(sprintf("  4 random effects : deviance %8.1f -> %8.1f   gain %8.1f  (from K = 3)\n",
            r3$deviance[1], r3$deviance[nrow(r3)], g_r3))
cat(sprintf("  2 random effects : deviance %8.1f -> %8.1f   gain %8.1f\n",
            r1$deviance[1], r1$deviance[nrow(r1)], g_r1))
cat(sprintf("\n  fixed df alone are worth %.2f%% of the mirrored gain (at 2 random effects)\n",
            100*g_r1/g_mir))
cat(sprintf("  R=3 BIC basin over K = 8..20: %.1f .. %.1f  (spread only %.1f)\n",
            min(r3$BIC[r3$K >= 8]), max(r3$BIC[r3$K >= 8]), diff(range(r3$BIC[r3$K >= 8]))))
cat(sprintf("  marginal R2 across all K   : %.4f .. %.4f  <- the population curve never moves\n",
            min(mir$r2_marg), max(mir$r2_marg)))
cat("\n============================ done ============================\n\n")
