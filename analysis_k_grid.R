# =============================================================================
# EXTENDED K GRID — push the spline df until BIC turns around
#
# Iteration 06 chose K by BIC over K = 1..10 and BIC never turned: its minimum
# sat on the last K of the grid, i.e. it reported the edge rather than selecting.
# The reason is arithmetic, and it is worth stating before reading any result:
#
#     lcmm computes  BIC = npm * log(ns) - 2*logLik   with ns = number of
#     SUBJECTS (31 in-vivo shoulders), not observations.
#
# Each extra spline df adds 2 parameters (one fixed effect, one diagonal random
# variance), so the penalty is only 2*log(31) = 6.87 BIC per df, against
# deviance gains still worth 49.5 at K = 9->10. Under an observation-based
# penalty (2*log(2277) = 15.46) BIC is still falling at K = 10 too.
#
# This script refits step 1 (IN-VIVO data only, no condition terms) over
# K = 1..20 with the settings of analysis_lrt.R untouched, and reports:
#   - where BIC turns under BOTH penalties (subjects, as lcmm does; observations)
#   - the L-method knee of the fit curve, over the longer grid
#   - which K stop converging, which is its own answer about how far to go
#
# NOTE ON WHAT BREAKS FIRST: at K = 20 the model asks for 21 diagonal random-
# effect variances from 31 shoulders. Expect non-convergence or variances
# collapsing to 0 well before K = 20. That ceiling is a real constraint on the
# model, not a failure of the script — it is reported, not hidden.
#
# Run:  Rscript analysis_k_grid.R
# =============================================================================

suppressMessages({ library(lcmm); library(splines) })
source("utils_lmethod.R")

INPUT   <- "monolix_st_frontal_dof2.csv"
S6      <- "06_natural_spline_hlme"
OUT     <- file.path("figures", S6)
CAP     <- 100                # max points per shoulder   (as analysis_lrt.R)
KS      <- 1:20               # <- the whole point of this script
TRIM_TO_OVERLAP <- TRUE
NPROC   <- max(1L, min(8L, parallel::detectCores() - 1L))
FIG_SCALE <- 2.2
LEVELS  <- c("in vivo", "ex vivo")
BLUE    <- "#377eb8"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

fig <- function(name, draw, w = 1300, h = 620, res = 130) {
  png(file.path(OUT, paste0(name, ".png")), width = round(w*FIG_SCALE),
      height = round(h*FIG_SCALE), res = round(res*FIG_SCALE))
  op <- par(no.readonly = TRUE); tryCatch(draw(), finally = { par(op); dev.off() })
  if (capabilities("cairo")) cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w/res, height = h/res)
  else pdf(file.path(OUT, paste0(name, ".pdf")), width = w/res, height = h/res)
  op <- par(no.readonly = TRUE); tryCatch(draw(), finally = { par(op); dev.off() })
  message("  saved ", file.path(S6, name), ".png/.pdf")
}
r2_conc <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)

# -----------------------------------------------------------------------------
# 1. LOAD — identical to analysis_lrt.R so the K are comparable
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
cat("\n============ EXTENDED K GRID — step 1 refitted over K = 1..20 ============\n")
cat(sprintf("in-vivo step-1 data: %d obs, %d shoulders, x in %.1f..%.1f\n",
            N_OBS, N_SUB, ov[1], ov[2]))
cat(sprintf("BIC penalty per extra spline df:  %.2f (lcmm, on %d subjects)   %.2f (on %d obs)\n",
            2*log(N_SUB), N_SUB, 2*log(N_OBS), N_OBS))
cat(sprintf("=> BIC turns when an extra df buys less deviance than that.\n\n"))

build <- function(K) {
  B <- ns(dt$TIME, df = K)
  x <- dt; nm <- paste0("ns", seq_len(ncol(B)))
  for (j in seq_along(nm)) x[[nm[j]]] <- B[, j]
  list(df = x, B = B, nm = nm)
}
f_red  <- function(nm) as.formula(paste("Y ~", paste(nm, collapse = " + ")))
f_rand <- function(nm) as.formula(paste("~", paste(nm, collapse = " + ")))

per_shoulder_rmse <- function(m, dat) {
  s <- split(seq_len(nrow(dat)), droplevels(dat$ID))
  sapply(s, function(i) sqrt(mean(m$pred$resid_ss[i]^2)))
}

# -----------------------------------------------------------------------------
# 2. FIT
# -----------------------------------------------------------------------------
rows <- list()
for (K in KS) {
  b   <- build(K)
  din <- b$df[b$df$cond == "in vivo", ]
  din$IDnum <- as.integer(droplevels(din$ID))
  t0 <- Sys.time()
  m <- try(hlme(fixed = f_red(b$nm), random = f_rand(b$nm), subject = "IDnum",
                ng = 1, idiag = TRUE, data = din, verbose = FALSE, nproc = NPROC),
           silent = TRUE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok <- !inherits(m, "try-error") && !is.null(m$conv) && m$conv == 1 && all(is.finite(m$best))
  npm <- if (ok) length(m$best) else 2*(K+1) + 1
  # a random-effect variance driven to ~0 means that df has nothing left to
  # explain between shoulders — worth seeing next to the BIC column
  nvar0 <- if (ok) sum(abs(m$best[grep("varcov|std", names(m$best), ignore.case = TRUE)]) < 1e-4) else NA
  rs <- if (ok) per_shoulder_rmse(m, din) else NA

  cat(sprintf("  K=%2d  npm %2d  %7.1fs  %s\n", K, npm, secs,
              if (ok) sprintf("BIC %8.1f  -2LL %8.1f  rmse max %5.2f med %5.2f",
                              m$BIC, -2*m$loglik, max(rs), median(rs))
              else paste("DID NOT CONVERGE", if (!inherits(m,"try-error") && !is.null(m$conv))
                                              paste0("(conv=", m$conv, ")") else "(error)")))
  flush.console()
  rows[[length(rows)+1]] <- data.frame(
    K = K, npm = npm, converged = ok,
    loglik = if (ok) m$loglik else NA, deviance = if (ok) -2*m$loglik else NA,
    BIC = if (ok) m$BIC else NA,
    BIC_obs = if (ok) -2*m$loglik + npm*log(N_OBS) else NA,
    n_var_at_zero = nvar0,
    rmse_max = if (ok) max(rs) else NA, rmse_median = if (ok) median(rs) else NA,
    r2_ss = if (ok) r2_conc(din$Y, m$pred$pred_ss) else NA,
    seconds = secs)
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
write.csv(tab, file.path(OUT, "08_bic_by_k_extended.csv"), row.names = FALSE)
message("\n  saved ", file.path(S6, "08_bic_by_k_extended.csv"))

# -----------------------------------------------------------------------------
# 3. DID BIC TURN?
# -----------------------------------------------------------------------------
o <- tab[tab$converged, ]
if (!nrow(o)) stop("nothing converged")
cat("\n---- results ----\n")
print(data.frame(K = o$K, npm = o$npm, deviance = round(o$deviance, 1),
                 d_dev = round(c(NA, diff(o$deviance)), 1),
                 BIC = round(o$BIC, 1), dBIC = round(c(NA, diff(o$BIC)), 1),
                 BIC_obs = round(o$BIC_obs, 1), dBIC_obs = round(c(NA, diff(o$BIC_obs)), 1),
                 rmse_max = round(o$rmse_max, 3)), row.names = FALSE)

turn <- function(v, K, label) {
  kb <- K[which.min(v)]
  if (kb == max(K))
    cat(sprintf("  %-22s minimum STILL at the last K tried (%d) — has not turned\n", label, kb))
  else
    cat(sprintf("  %-22s minimum at K = %d  (turns: BIC rises from K = %d on)\n", label, kb, kb + 1))
  kb
}
cat("\n")
K_bic     <- turn(o$BIC,     o$K, "BIC (lcmm, subjects)")
K_bic_obs <- turn(o$BIC_obs, o$K, "BIC (observations)")
if (any(!tab$converged))
  cat(sprintf("  non-convergence at K = %s  <- the model's own ceiling\n",
              paste(tab$K[!tab$converged], collapse = ", ")))

cmp <- lmethod_vs_bic(o$K, metric = o$loglik, BIC = o$BIC, metric_is_loglik = TRUE)
print(cmp)

cat("\n  knee of each fit curve over the extended grid:\n")
curves <- list("deviance"          = lmethod(o$K, o$deviance),
               "log deviance"      = lmethod(o$K, o$deviance, transform = "log"),
               "median RMSE (deg)" = lmethod(o$K, o$rmse_median),
               "max RMSE (deg)"    = lmethod(o$K, o$rmse_max))
print(data.frame(curve = names(curves), knee = sapply(curves, `[[`, "knee"),
                 strength = round(sapply(curves, `[[`, "strength"), 3)), row.names = FALSE)

cat("\n  does the knee survive extending the grid?\n")
print(lmethod_sensitivity(o$K, o$deviance), digits = 4, row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. FIGURES
# -----------------------------------------------------------------------------
fig("08_lmethod_vs_bic_extended", function() plot_lmethod_vs_bic(cmp))
fig("08_sensitivity_extended", function() {
  op <- par(mar = c(4.4, 4.6, 3.4, 1)); on.exit(par(op), add = TRUE)
  plot_lmethod_sensitivity(o$K, o$deviance)
}, w = 800)
fig("08_bic_two_penalties", function() {
  op <- par(mar = c(4.4, 4.6, 3.4, 1)); on.exit(par(op), add = TRUE)
  yr <- range(o$BIC, o$BIC_obs)
  plot(o$K, o$BIC, type = "n", ylim = yr, xaxt = "n", xlab = "natural-spline df (K)",
       ylab = "BIC (lower is better)",
       main = "BIC depends on what n means — subjects (lcmm) vs observations")
  axis(1, at = o$K); grid(col = "grey90", lty = 1)
  lines(o$K, o$BIC, lwd = 2.4, col = "grey40")
  points(o$K, o$BIC, pch = 21, bg = "white", col = "grey40", cex = 1.1, lwd = 1.3)
  lines(o$K, o$BIC_obs, lwd = 2.4, col = "#d95f02", lty = 2)
  points(o$K, o$BIC_obs, pch = 21, bg = "white", col = "#d95f02", cex = 1.1, lwd = 1.3)
  abline(v = K_bic, lty = 3, col = "grey40", lwd = 2)
  abline(v = K_bic_obs, lty = 3, col = "#d95f02", lwd = 2)
  abline(v = cmp$K_knee, lty = 3, col = BLUE, lwd = 2.4)
  legend("topright", bty = "n", cex = 0.8, lwd = 2.4, lty = c(1, 2, 3),
         col = c("grey40", "#d95f02", BLUE),
         legend = c(sprintf("BIC on %d subjects (lcmm): K = %d", N_SUB, K_bic),
                    sprintf("BIC on %d observations: K = %d", N_OBS, K_bic_obs),
                    sprintf("L-method knee: K = %d", cmp$K_knee)))
}, w = 900)

write.csv(cmp$table, file.path(OUT, "08_lmethod_vs_bic_extended.csv"), row.names = FALSE)
message("  saved ", file.path(S6, "08_lmethod_vs_bic_extended.csv"))
cat("\n============================ done ============================\n\n")
