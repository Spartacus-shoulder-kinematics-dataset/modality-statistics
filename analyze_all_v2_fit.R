# =============================================================================
# GENERALIZED SWEEP v2 — THE FIT HALF
#
# Fits one hlme per joint x humeral_motion x degree_of_freedom cell and writes
# everything the figures need to figures/generalized_v2/cache/v2_fit.rds.
# Draws nothing. See analyze_all_v2.R for the model and the per-cell recipe.
#
# Run (from the repository root):  Rscript analyze_all_v2_fit.R      # ~12 min
# Then, as often as you like:      Rscript analyze_all_v2_figures.R  # seconds
#
# Outputs:
#   cache/v2_fit.rds        curves, differences, coefficients and covariances
#   00_master_summary.csv   one row per cell, incl. Wald counts and joint test
#   00_wald_tests.csv       one row per COEFFICIENT, every cell
# =============================================================================

source("analyze_all_v2_common.R")
suppressMessages(library(lcmm))

NPROC <- max(1L, min(8L, parallel::detectCores() - 1L))

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
  key <- cell_key(jt, mo, df_)
  sub <- D[D$joint == jt & D$humeral_motion == mo & D$degree_of_freedom == df_, ]
  row <- data.frame(joint = jt, motion = mo, dof = df_, n_rows = nrow(sub),
                    n_units = 0, n_ex = 0, n_in = 0, st_ex = 0, st_in = 0,
                    mode = "skip", x_lo = NA, x_hi = NA,
                    x_ex_lo = NA, x_ex_hi = NA, x_in_lo = NA, x_in_hi = NA,
                    knots = NA, n_fitted = 0,
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
  # How far each condition reaches ON ITS OWN, before any trimming. x_lo/x_hi below
  # is the OVERLAP, which is what the model is fitted on; these two spans are what
  # the overlap is an intersection OF, and the significance map draws them in the
  # condition colours. Computed here, above the early return, so that skipped and
  # single-condition cells carry them too — that is what lets a "not compared" row
  # still show why (one condition covering ground the other never reaches).
  # Safe on `sub` rather than the thinned `st`: thin_by_unit keeps each unit's first
  # and last row, so per-unit ranges — all support_range looks at — are unchanged.
  for (cc in c("ex vivo", "in vivo")) {
    k <- sub$cond == cc
    r <- if (any(k)) support_range(sub$TIME[k], sub$ID[k]) else c(NA_real_, NA_real_)
    if (cc == "ex vivo") { row$x_ex_lo <- r[1]; row$x_ex_hi <- r[2] }
    else                 { row$x_in_lo <- r[1]; row$x_in_hi <- r[2] }
  }
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
  # curves/dif/knots/raw are what the figures draw today. bf/Vf plus the basis
  # recipe (knots + full_rng + Kc) are the sufficient statistics for anything
  # they might draw tomorrow: a different prediction grid or a different
  # contrast can be rebuilt from these without refitting, via
  #   ns(x, knots = knots, Boundary.knots = full_rng)
  # The hlme object itself is NOT kept — it carries per-observation predictions
  # and is far larger than the parts actually needed.
  CUR[[key]] <- list(curves = curves, dif = dif, knots = kn, raw = st,
                     bf = bf, Vf = Vf, full_rng = full_rng, Kc = Kc,
                     compare = compare, ov = ov)
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
wald <- if (length(wald_rows)) do.call(rbind, wald_rows) else NULL
if (!is.null(wald))
  write.csv(wald, file.path(OUT, "00_wald_tests.csv"), row.names = FALSE)

# Written AFTER the BH adjustment above: the planche annotations read
# joint_p_fdr to draw their stars, so a cache saved inside the loop would be
# missing exactly the column the figures need.
saveRDS(list(cells = CUR, master = master, wald = wald, meta = fit_meta()), CACHE)

cat(sprintf("\nswept %d combos in %.1f min | compare %d | single %d | skip %d | failed %d\n",
            nrow(master), as.numeric(difftime(Sys.time(), t_start, units="mins")),
            sum(master$mode=="compare"), sum(master$mode=="single"),
            sum(master$mode=="skip"), sum(!is.na(master$converged) & !master$converged)))
cat(sprintf("cached %d fitted cells -> %s (%.1f MB)\n",
            length(CUR), CACHE, file.size(CACHE)/1024^2))
cat("Now run:  Rscript analyze_all_v2_figures.R\n")
