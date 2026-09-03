# =============================================================================
# Full analysis in R — Scapulothoracic rhythm: in vivo vs ex vivo
# Pedagogical, Monolix-free path. Companion to METHODOLOGY.md.
#
# Figures are kept as ITERATION SETS (see figures/README.md for the story):
#   figures/00_data_exploration/   model-agnostic "look at the data" figures
#   figures/01_sigmoid_nlme/       first attempt: parametric sigmoid NLME
#   figures/02_spline_mixed_model/ the fix: penalised spline mixed model + AR(1)
# We keep every set on purpose so the reasoning (why we moved on) is visible.
#
# Run:  Rscript analysis_shoulder.R
# Prereq: run  python3 prepare_monolix_data.py  first (creates the input CSV).
# Packages (all preinstalled): mgcv, nlme, dplyr. Base graphics for plots.
# =============================================================================

suppressMessages({ library(mgcv); library(nlme); library(dplyr) })

INPUT   <- "monolix_st_frontal_dof2.csv"
FIG_DIR <- "figures"
S0 <- "00_data_exploration"; S1 <- "01_sigmoid_nlme"; S2 <- "02_spline_mixed_model"
COL <- c("ex vivo" = "#d95f02", "in vivo" = "#1b9e77")

# helper: save a PNG into a given iteration subfolder
fig <- function(subdir, name, expr, w = 1100, h = 750, res = 130) {
  dir.create(file.path(FIG_DIR, subdir), showWarnings = FALSE, recursive = TRUE)
  png(file.path(FIG_DIR, subdir, name), width = w, height = h, res = res)
  on.exit(dev.off()); force(expr); message("  saved ", file.path(subdir, name))
}

# -----------------------------------------------------------------------------
# 1. LOAD  (long format: one row per (unit, x, y))
# -----------------------------------------------------------------------------
if (!file.exists(INPUT))
  stop("Missing ", INPUT, " — run  python3 prepare_monolix_data.py  first.")

d <- read.csv(INPUT, stringsAsFactors = FALSE)
d$cond <- factor(ifelse(d$in_vivo %in% c("True","TRUE","1","yes"),
                        "in vivo", "ex vivo"), levels = c("ex vivo","in vivo"))
d$ID <- factor(d$ID)
d <- d[is.finite(d$TIME) & is.finite(d$Y), ]
d <- d[order(d$ID, d$TIME), ]                   # order within unit (needed for AR1)

cat("\n================ DATA SUMMARY ================\n")
cat(sprintf("observations : %d\nunits (ID)   : %d\n", nrow(d), nlevels(d$ID)))
cat("rows per condition:\n"); print(table(d$cond))
cat(sprintf("x range: %.1f..%.1f | y range: %.1f..%.1f\n",
            min(d$TIME), max(d$TIME), min(d$Y), max(d$Y)))
# DATA-OWNER NOTE: values are labelled unit=="rad" but magnitudes look like
# DEGREES. Flag for provenance; modelling is unaffected.

# =============================================================================
# 2. ITERATION 00 — DATA APPRECIATION ("look before you model")
# =============================================================================
fig(S0, "01_raw_scatter.png", {
  plot(d$TIME, d$Y, col = adjustcolor(COL[as.character(d$cond)], 0.35), pch = 16,
       cex = 0.5, xlab = "Humerothoracic elevation (x)", ylab = "Scapulothoracic angle (y)",
       main = "Raw data: every observation")
  legend("topleft", names(COL), col = COL, pch = 16, bty = "n")
})
fig(S0, "02_spaghetti_by_unit.png", {
  plot(NA, xlim = range(d$TIME), ylim = range(d$Y), xlab = "elevation (x)",
       ylab = "ST angle (y)", main = "One curve per unit (article x shoulder)")
  for (u in levels(d$ID)) { du <- d[d$ID == u, ]
    lines(du$TIME, du$Y, col = adjustcolor(COL[as.character(du$cond[1])], 0.5)) }
  legend("topleft", names(COL), col = COL, lwd = 2, bty = "n")
})
fig(S0, "03_x_sampling_density.png", {
  br <- seq(floor(min(d$TIME)), ceiling(max(d$TIME)), length.out = 40)
  h1 <- hist(d$TIME[d$cond=="ex vivo"], breaks = br, plot = FALSE)
  h2 <- hist(d$TIME[d$cond=="in vivo"], breaks = br, plot = FALSE)
  plot(h1, col = adjustcolor(COL["ex vivo"],0.5), border = NA, xlim = range(br),
       ylim = c(0, max(h1$counts,h2$counts)), xlab = "elevation (x)",
       main = "Where is x sampled? (irregular => not wide-format LGCM)")
  plot(h2, col = adjustcolor(COL["in vivo"],0.5), border = NA, add = TRUE)
  legend("topright", names(COL), fill = adjustcolor(COL,0.5), border = NA, bty = "n")
})
pts <- as.integer(table(d$ID))
cnd <- sapply(levels(d$ID), function(u) as.character(d$cond[d$ID==u][1]))
fig(S0, "04_points_per_unit.png", {
  o <- order(pts)
  barplot(pts[o], col = COL[cnd[o]], border = NA, names.arg = rep("",length(o)),
          ylab = "points in unit", xlab = "units (sorted)",
          main = "Points per unit (11..2987: extreme imbalance)")
  legend("topleft", names(COL), fill = COL, border = NA, bty = "n")
})
fig(S0, "05_binned_mean_sd.png", {
  d$bin <- cut(d$TIME, breaks = seq(min(d$TIME), max(d$TIME), length.out = 25))
  agg <- d %>% group_by(cond, bin) %>%
    summarise(x = mean(TIME), m = mean(Y), s = sd(Y), .groups = "drop") %>% filter(is.finite(x))
  plot(NA, xlim = range(d$TIME), ylim = range(d$Y), xlab = "elevation (x)",
       ylab = "ST angle (y)", main = "Binned mean +/- SD by condition")
  for (cc in levels(d$cond)) { a <- agg[agg$cond==cc, ]; a <- a[order(a$x), ]
    arrows(a$x, a$m-a$s, a$x, a$m+a$s, length=0.02, angle=90, code=3, col = adjustcolor(COL[cc],0.4))
    lines(a$x, a$m, col = COL[cc], lwd = 2); points(a$x, a$m, col = COL[cc], pch = 16) }
  legend("topleft", names(COL), col = COL, lwd = 2, bty = "n")
})

# =============================================================================
# 3. ITERATION 01 — PARAMETRIC SIGMOID NLME (first attempt; we move on from it)
#    A 4-parameter logistic (SSfpl). On THIS data the angle rises monotonically
#    and never plateaus, so the asymptotes fall outside the observed range.
#    -> random effects on all 4 params are non-identifiable (SDs explode to
#       ~1e5, all fixed-effect p-values ~1; see figures/README.md).
#    The salvageable version uses a random INTERCEPT only; it still fits WORST.
# =============================================================================
cat("\n============ ITERATION 01: parametric sigmoid NLME ============\n")
ns <- nls(Y ~ SSfpl(TIME, A, B, xmid, scal), data = d)
cat("pooled sigmoid params (note xmid at the data edge, B far below data):\n")
print(round(coef(ns), 2))

sig_ok <- FALSE
try({
  m_sig <- nlme(Y ~ SSfpl(TIME, A, B, xmid, scal), data = d,
                fixed = A + B + xmid + scal ~ 1, random = A ~ 1 | ID, start = coef(ns),
                control = nlmeControl(maxIter = 200, msMaxIter = 200, returnObject = TRUE))
  cat(sprintf("random-intercept sigmoid BIC = %.0f\n", BIC(m_sig))); sig_ok <- TRUE
}, silent = TRUE)

fig(S1, "01_sigmoid_pooled_fit.png", {
  plot(d$TIME, d$Y, col = adjustcolor("grey60", 0.20), pch = 16, cex = 0.45,
       xlab = "elevation (x)", ylab = "ST angle (y)",
       main = "Attempt 1: parametric 4-parameter logistic")
  xx <- seq(min(d$TIME), max(d$TIME), length.out = 300)
  lines(xx, predict(ns, newdata = data.frame(TIME = xx)), lwd = 3, col = "purple")
  co <- coef(ns)
  legend("topleft", bty = "n",
         legend = c(sprintf("xmid = %.0f  (at data edge)", co["xmid"]),
                    sprintf("B = %.0f  (asymptote below data)", co["B"]),
                    "=> asymptotes out of range => hard to identify"))
})

# =============================================================================
# 4. ITERATION 02 — PRIMARY: spline mixed model via mgcv, WITH AR(1)
#    Y ~ condO + s(TIME) + s(TIME, by=condO) + s(ID, bs="re")
#      s(TIME)             : ex-vivo (reference) nonlinear trajectory
#      s(TIME, by=condO)   : in-vivo DIFFERENCE smooth  (condO is an ordered factor)
#      s(ID, bs="re")      : per-unit random intercept (the mixed part)
# =============================================================================
cat("\n============ ITERATION 02: spline mixed model (mgcv) ============\n")
d$condO    <- as.ordered(d$cond)
d$ar_start <- c(TRUE, d$ID[-1] != d$ID[-nrow(d)])

form    <- Y ~ condO + s(TIME) + s(TIME, by = condO) + s(ID, bs = "re")
m_naive <- bam(form, data = d, method = "fREML", discrete = TRUE)
r   <- resid(m_naive)
# lag-1 autocorrelation computed WITHIN each shoulder (do not cross boundaries)
num <- den <- 0
for (u in levels(d$ID)) { ru <- r[d$ID == u]
  if (length(ru) > 1) { num <- num + sum(head(ru,-1)*tail(ru,-1)); den <- den + sum(ru^2) } }
rho <- num / den
cat(sprintf("lag-1 residual autocorrelation rho = %.3f  (=> AR(1) needed)\n", rho))
m_ar <- bam(form, data = d, method = "fREML", discrete = TRUE,
            rho = rho, AR.start = d$ar_start)

st_naive <- summary(m_naive)$s.table
st_ar    <- summary(m_ar)$s.table
pt_ar    <- summary(m_ar)$p.table
cat("\nSmooth terms WITHOUT AR(1) (over-optimistic):\n"); print(round(st_naive, 3))
cat("\nSmooth terms WITH AR(1) (honest inference):\n");    print(round(st_ar, 3))
cat("\nParametric terms WITH AR(1) (condO.L = in-vivo level shift):\n"); print(round(pt_ar, 3))

# population-curve predictions (exclude the random effect)
grid <- expand.grid(TIME = seq(min(d$TIME), max(d$TIME), length.out = 300),
                    condO = factor(levels(d$condO), levels = levels(d$condO), ordered = TRUE))
grid$ID <- d$ID[1]
pp <- predict(m_ar, newdata = grid, exclude = "s(ID)", newdata.guaranteed = TRUE, se.fit = TRUE)
grid$fit <- pp$fit; grid$lo <- pp$fit - 1.96*pp$se.fit; grid$hi <- pp$fit + 1.96*pp$se.fit

# ---- iteration-01 vs iteration-02 overlay: why the spline wins -----------
fig(S1, "02_parametric_vs_spline.png", {
  plot(d$TIME, d$Y, col = adjustcolor("grey60", 0.15), pch = 16, cex = 0.4,
       xlab = "elevation (x)", ylab = "ST angle (y)",
       main = "Why we moved on: sigmoid (dashed) vs spline (solid)")
  xx  <- seq(min(d$TIME), max(d$TIME), length.out = 300)
  gg  <- grid %>% group_by(TIME) %>% summarise(fit = mean(fit), .groups = "drop")
  lines(gg$TIME, gg$fit, lwd = 3, col = "black")
  lines(xx, predict(ns, newdata = data.frame(TIME = xx)), lwd = 3, lty = 2, col = "purple")
  legend("topleft", c("spline (mgcv)", "sigmoid (NLME)"),
         col = c("black","purple"), lty = c(1,2), lwd = 3, bty = "n")
})

# ---- HEADLINE: population trajectory per condition -------------------------
fig(S2, "01_population_curves_by_condition.png", {
  plot(d$TIME, d$Y, col = adjustcolor(COL[as.character(d$cond)], 0.15), pch = 16,
       cex = 0.45, xlab = "Humerothoracic elevation (x)", ylab = "Scapulothoracic angle (y)",
       main = "Population trajectory per condition (spline mixed model, AR1)")
  for (cc in levels(d$condO)) { g <- grid[grid$condO == cc, ]
    polygon(c(g$TIME, rev(g$TIME)), c(g$lo, rev(g$hi)), col = adjustcolor(COL[cc], 0.25), border = NA)
    lines(g$TIME, g$fit, col = COL[cc], lwd = 3) }
  legend("topleft", paste("population:", names(COL)), col = COL, lwd = 3, bty = "n")
})

# ---- the hypothesis test, visualised: in vivo - ex vivo difference + 95% CI
# Correct contrast via the linear-predictor matrix: taking (in vivo - ex vivo)
# at the same grid cancels the shared intercept, s(TIME) and s(ID) columns, and
# leaves exactly the condition effect (level shift + difference smooth).
fig(S2, "02_difference_invivo_minus_exvivo.png", {
  gx <- seq(min(d$TIME), max(d$TIME), length.out = 300)
  nd_in <- data.frame(TIME = gx, ID = d$ID[1],
                      condO = factor("in vivo", levels = levels(d$condO), ordered = TRUE))
  nd_ex <- data.frame(TIME = gx, ID = d$ID[1],
                      condO = factor("ex vivo", levels = levels(d$condO), ordered = TRUE))
  Xd   <- predict(m_ar, nd_in, type = "lpmatrix", newdata.guaranteed = TRUE) -
          predict(m_ar, nd_ex, type = "lpmatrix", newdata.guaranteed = TRUE)
  dif  <- as.vector(Xd %*% coef(m_ar))
  sef  <- sqrt(rowSums((Xd %*% vcov(m_ar)) * Xd))
  plot(gx, dif, type = "n", ylim = range(c(dif-2*sef, dif+2*sef, 0)),
       xlab = "Humerothoracic elevation (x)", ylab = "in vivo  -  ex vivo  (ST angle)",
       main = "Estimated difference (in vivo vs ex vivo) +/- 95% CI")
  polygon(c(gx, rev(gx)), c(dif-1.96*sef, rev(dif+1.96*sef)),
          col = adjustcolor("#377eb8", 0.25), border = NA)
  lines(gx, dif, lwd = 3, col = "#377eb8")
  abline(h = 0, lty = 2, col = "grey40")     # where the band excludes 0 => differ
})

fig(S2, "03_diagnostics.png", {
  op <- par(mfrow = c(1,2), mar = c(4,4,2,1)); on.exit(par(op), add = TRUE)
  fv <- fitted(m_ar); rs <- resid(m_ar)
  plot(fv, rs, pch = 16, cex = 0.4, col = adjustcolor("grey20",0.25),
       xlab = "fitted", ylab = "residual", main = "Residuals vs fitted"); abline(h = 0, col = "red")
  qqnorm(rs, pch = 16, cex = 0.4, col = adjustcolor("grey20",0.25), main = "Normal Q-Q"); qqline(rs, col = "red")
})
fig(S2, "04_random_effects.png", {
  re <- coef(m_ar)[grep("s\\(ID\\)", names(coef(m_ar)))]
  hist(re, col = "grey80", border = "white", xlab = "random intercept (unit)",
       main = "Between-unit random intercepts s(ID)")
})
fig(S2, "05_autocorrelation.png", {
  ac <- acf(resid(m_naive), lag.max = 40, plot = FALSE)
  plot(ac$lag, ac$acf, type = "h", lwd = 2, xlab = "lag", ylab = "residual ACF",
       main = sprintf("Why AR(1): naive residual autocorrelation, lag-1 = %.2f", rho)); abline(h = 0)
})
fig(S2, "06_model_comparison.png", {
  bics <- c(spline_AR1 = BIC(m_ar), spline_naive = BIC(m_naive),
            if (sig_ok) sigmoid_NLME = BIC(m_sig))
  bp <- barplot(bics, col = c("#377eb8","grey70","purple")[seq_along(bics)], las = 1,
                cex.names = 0.85, ylab = "BIC (lower is better)",
                main = "Model comparison (BIC)")
  text(bp, pmax(bics, 0), round(bics), pos = 3, xpd = TRUE)
})

# =============================================================================
# 5. RESULTS DIGEST (in the winning iteration's folder)
# =============================================================================
sink(file.path(FIG_DIR, S2, "00_results_summary.txt"))
cat("Scapulothoracic rhythm — in vivo vs ex vivo — R analysis digest\n")
cat("================================================================\n\n")
cat(sprintf("Data: %d obs, %d units. lag-1 residual autocorrelation = %.3f.\n\n",
            nrow(d), nlevels(d$ID), rho))
cat("WINNING MODEL: spline mixed model (mgcv) with AR(1) residuals\n")
cat("  Y ~ condO + s(TIME) + s(TIME, by=condO) + s(ID, bs='re')\n\n")
cat("Smooth terms (AR(1)-corrected, honest inference):\n"); print(round(st_ar, 4))
cat("\nParametric terms (AR(1)):\n"); print(round(pt_ar, 4))
cat("\nHOW TO READ IT:\n")
cat(" * s(TIME):condOin vivo -> DIFFERENCE in trajectory SHAPE (small p = differ).\n")
cat(" * condO.L              -> constant LEVEL shift between conditions.\n")
cat(" * figure 02 shows the total difference (with 95% CI) across elevation:\n")
cat("   where the band excludes 0, the conditions differ at that elevation.\n\n")
cat(sprintf("BIC: spline+AR1 = %.0f | spline naive = %.0f%s\n",
            BIC(m_ar), BIC(m_naive),
            if (sig_ok) sprintf(" | sigmoid NLME = %.0f", BIC(m_sig)) else ""))
cat("\nCAVEATS: (1) values labelled 'rad' look like degrees - verify provenance.\n")
cat(" (2) 44 units, heavily imbalanced dense ex-vivo curves - between-study\n")
cat("     variance modest to estimate; interpret cautiously.\n")
cat(" (3) AR(1) handles within-curve autocorrelation approximately.\n")
sink()

cat("\nDONE. See figures/ (00_data_exploration, 01_sigmoid_nlme, 02_spline_mixed_model)\n")
cat("and figures/README.md for the iteration story.\n")
