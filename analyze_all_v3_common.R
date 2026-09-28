# =============================================================================
# GENERALIZED SWEEP v3 — shared configuration, helpers and the per-cell recipe
#
# Sourced by all four halves: analyze_all_v3_step{1,2,3}.R and
# analyze_all_v3_figures.R. Everything the steps have to agree on lives here —
# where the files are, what defines a cell, how a cell's rows and basis are
# built, the palette and the panel ordering.
#
# THE PREPARATION IS SHARED ON PURPOSE. Step 1 chooses K on the reference
# condition, step 2 fits H1 at that K on everything, step 3 fits H0 at that K on
# everything and differences the log-likelihoods. That difference is meaningless
# unless all three built IDENTICAL rows and an IDENTICAL basis, so prepare_cell()
# and cell_basis() are the single source for both and every step calls them.
#
# lcmm is deliberately NOT loaded here — only the fitting steps need it, and the
# point of the split is that redrawing a figure costs nothing.
# =============================================================================

suppressMessages({ library(splines) })

LONG      <- "spartacus_angles_long.csv"
OUT       <- file.path("figures", "generalized_v3")
CACHE_DIR <- file.path(OUT, "cache")
STEP_DIR  <- function(s) file.path(CACHE_DIR, paste0("step", s))
CACHE     <- file.path(CACHE_DIR, "v3_fit.rds")
CELL_DIR  <- file.path(OUT, "cells")

# ---- model constants --------------------------------------------------------
# KS is the whole selection rule: fit the reduced model at each of these, take
# the BIC minimum. No cap, no adaptive stopping — a cell whose minimum lands on
# an edge of this grid is recorded as bic_edge and drawn as such rather than
# quietly reported as a choice.
KS           <- 1:20
# V3_KS='1:3' truncates the grid for a wiring smoke test. It is part of the
# fingerprint, so a cache built under a short grid reads as stale rather than
# being mistaken for a real selection.
if (nzchar(Sys.getenv("V3_KS")))
  KS <- eval(parse(text = Sys.getenv("V3_KS")))
CAP          <- 100    # max points per shoulder
MIN_PER_COND <- 3      # >=3 shoulders in BOTH conditions -> compare
MIN_TOTAL    <- 4      # >=4 shoulders total              -> at least single
MIN_ROWS     <- 20     # below this a cell is not modelled at all
# How many shoulders must span an x before it counts as supported. v2 used 2.
# At 1, support_range() degenerates to the absolute min/max, which is Melanie's
# TODO — "garder toute la plage des abscisses meme si que deux individus": keep
# the whole abscissa even where only one shoulder reaches. Cells then fit over a
# wider x, and K can change with it, since the sparse tails are exactly where
# extra spline df get spent.
MIN_UNITS_X  <- 1
TRIM_TO_OVERLAP <- TRUE
# in vivo FIRST => in vivo is the REFERENCE, the departure term is ex vivo, and
# the difference curve is EX VIVO - IN VIVO. Mirrored against v2, which had ex
# vivo as the reference.
LEVELS       <- c("in vivo", "ex vivo")
# Threads INSIDE one hlme fit. When sharding across processes this must come
# down, or N shards x 8 threads oversubscribes the machine and every shard runs
# slower: V3_NPROC=4 with V3_NSHARD=5 keeps 20 of 22 cores busy and no more.
# Not part of fit_params(): it changes how long the sweep takes, not what comes
# out of it.
NPROC        <- max(1L, min(8L, parallel::detectCores() - 1L))
if (nzchar(Sys.getenv("V3_NPROC")))
  NPROC <- max(1L, as.integer(Sys.getenv("V3_NPROC")))

COL   <- c("ex vivo" = "#d95f02", "in vivo" = "#1b9e77")
BLUE  <- "#377eb8"

JOINT_ORDER <- c("sternoclavicular", "acromioclavicular", "scapulothoracic", "glenohumeral")
DOF_LEGEND <- list(
  glenohumeral      = c("Plane of elevation", "Elevation(-)/Depression(+)", "Int(+)/ext(-) rotation"),
  scapulothoracic   = c("Protraction(+)/retraction(-)", "Medial(+)/lateral(-) rot", "Posterior(+)/anterior(-) tilt"),
  acromioclavicular = c("Protraction(+)/retraction(-)", "Medial(+)/lateral(-) rot", "Posterior(+)/anterior(-) tilt"),
  sternoclavicular  = c("Protraction(+)/retraction(-)", "Depression(+)/elevation(-)", "Backwards(+)/forward(-) rot"))
MOTION_ORDER <- c("frontal plane elevation", "scapular plane elevation", "sagittal plane elevation",
                  "horizontal flexion", "internal-external rotation 0 degree-abducted",
                  "internal-external rotation 90 degree-abducted")

# The two families the DoF-grouped significance maps are built from. MOTION_ORDER
# above is the reading order of the planches (one plate per motion); these are the
# INNER axis of a plate whose rows are joint x DoF, so that the effect of the same
# DoF can be compared across the planes the arm moves in.
ELEV_MOTIONS <- MOTION_ORDER[1:3]
ROT_MOTIONS  <- MOTION_ORDER[5:6]
MOTION_SHORT <- c("frontal plane elevation"  = "frontal",
                  "scapular plane elevation" = "scapular",
                  "sagittal plane elevation" = "sagittal",
                  "horizontal flexion"       = "horiz. flexion",
                  "internal-external rotation 0 degree-abducted"  = "IER 0°",
                  "internal-external rotation 90 degree-abducted" = "IER 90°")

for (d in c(OUT, CACHE_DIR, CELL_DIR, STEP_DIR(1), STEP_DIR(2), STEP_DIR(3)))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- small helpers ----------------------------------------------------------
slug  <- function(s) gsub("(^_|_$)", "", gsub("[^a-z0-9]+", "_", tolower(s)))
capj  <- function(j) paste0(toupper(substring(j, 1, 1)), substring(j, 2))
stars <- function(p) if (is.na(p)) "" else if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else "ns"
# How a cell is addressed in memory; cell_slug is how it is addressed on disk.
cell_key  <- function(jt, mo, dof) paste(jt, mo, dof, sep = "|")
cell_slug <- function(jt, mo, dof) sprintf("%s_%s_dof%d", slug(jt), slug(mo), dof)

# concordance R^2 against the 1:1 line — NOT a squared correlation, so a shoulder
# fitted with a constant offset cannot score 0.99 while sitting off the diagonal
r2_conc <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)

norm_stats <- function(r) {
  r <- r[is.finite(r)]
  if (length(r) < 4 || sd(r) == 0) return(c(skew = NA, kurtosis = NA, shapiro_p = NA))
  z <- (r - mean(r))/sd(r)
  sw <- tryCatch({ s <- if (length(r) > 5000) sample(r, 5000) else r
                   shapiro.test(s)$p.value }, error = function(e) NA_real_)
  c(skew = mean(z^3), kurtosis = mean(z^4) - 3, shapiro_p = sw)
}

# ---- the house renderer -----------------------------------------------------
# Render every figure TWICE: a high-resolution PNG and a vector PDF.
# The draw argument is a FUNCTION, not a braced expression, so it can be called
# once per device — a promise would only evaluate on the first device and the
# second file would come out blank.
FIG_SCALE <- 2.2          # pixel multiplier; res scales with it so text stays proportional
figout <- function(base, draw, w = 1100, h = 750, res = 130, dir = OUT) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  png(file.path(dir, paste0(base, ".png")),
      width = round(w*FIG_SCALE), height = round(h*FIG_SCALE), res = round(res*FIG_SCALE))
  op <- par(no.readonly = TRUE); tryCatch(draw(), finally = { par(op); dev.off() })
  # cairo_pdf, NOT pdf(): the legacy device uses Type 1 fonts with no glyph for
  # the circles/arrows/Greek used here, so those characters silently vanish.
  if (capabilities("cairo")) cairo_pdf(file.path(dir, paste0(base, ".pdf")),
                                       width = w/res, height = h/res)
  else pdf(file.path(dir, paste0(base, ".pdf")), width = w/res, height = h/res)
  op <- par(no.readonly = TRUE); tryCatch(draw(), finally = { par(op); dev.off() })
}

# ---- cache metadata ---------------------------------------------------------
# The constants that actually change what a fit means. NPROC is not among them:
# it changes how long the sweep takes, not what comes out of it.
fit_params <- function()
  list(KS = KS, CAP = CAP, MIN_PER_COND = MIN_PER_COND, MIN_TOTAL = MIN_TOTAL,
       MIN_ROWS = MIN_ROWS, MIN_UNITS_X = MIN_UNITS_X,
       TRIM_TO_OVERLAP = TRIM_TO_OVERLAP, LEVELS = LEVELS)

fit_meta <- function() {
  fi <- file.info(LONG)
  list(params = fit_params(),
       input = list(file = LONG, size = fi$size, mtime = fi$mtime),
       fitted_at = Sys.time(),
       R_version = getRversion(),
       lcmm_version = tryCatch(as.character(utils::packageVersion("lcmm")),
                               error = function(e) NA_character_))
}

# Once fitting and drawing are separate commands, nothing stops you editing a
# constant and redrawing from a cache built under the old value. Warn — loudly,
# by name — rather than stop: a deliberate mismatch is sometimes what you want.
check_meta <- function(meta) {
  now <- fit_params(); was <- meta$params
  if (!length(was)) {
    warning("cache carries no settings fingerprint; cannot check it against the current ",
            "constants. Re-run the steps if in doubt.", call. = FALSE)
    return(invisible(NULL))
  }
  bad <- names(now)[!mapply(identical, now[names(now)], was[names(now)])]
  if (length(bad))
    warning("cache was built with different settings: ",
            paste(sprintf("%s = %s (cache) vs %s (now)", bad,
                          sapply(bad, function(k) paste(format(was[[k]]), collapse = ",")),
                          sapply(bad, function(k) paste(format(now[[k]]), collapse = ","))),
                   collapse = "; "),
            "\n  re-run the affected step to rebuild.", call. = FALSE)
  fi <- file.info(LONG)
  if (!is.na(fi$size) && !identical(fi$size, meta$input$size))
    warning(LONG, " has changed since the cache was built (", meta$input$size,
            " -> ", fi$size, " bytes).", call. = FALSE)
  invisible(NULL)
}

# ---- per-cell cache ---------------------------------------------------------
# One RDS per cell per step, written the moment that cell's step finishes, so a
# killed run loses at most the cell it was inside. Each file carries the
# fingerprint it was built under; a mismatch makes it stale and it is refitted.
step_file <- function(s, jt, mo, dof) file.path(STEP_DIR(s), paste0(cell_slug(jt, mo, dof), ".rds"))

cell_save <- function(s, jt, mo, dof, obj) {
  obj$meta <- fit_meta()
  saveRDS(obj, step_file(s, jt, mo, dof))
  invisible(obj)
}

# Returns NULL when there is nothing usable on disk: absent, unreadable, or built
# under different constants. Callers treat NULL as "this cell still has to run".
cell_load <- function(s, jt, mo, dof, quiet = TRUE) {
  f <- step_file(s, jt, mo, dof)
  if (!file.exists(f)) return(NULL)
  o <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(o)) return(NULL)
  if (!identical(o$meta$params, fit_params())) {
    if (!quiet) message("  stale cache (settings changed): ", basename(f))
    return(NULL)
  }
  o
}

# ---- which cells does THIS process run? -------------------------------------
# Sharding across separate R processes rather than mclapply: hlme(nproc=) already
# forks internally, and nesting the two makes a mess. Four terminals with
# V3_NSHARD=4 V3_SHARD=0..3 each write their own per-cell files into the shared
# cache dir — one writer per file, so no locking is needed.
#   V3_ONLY='<joint>|<motion>|<dof>'  runs exactly one cell (verification)
#   V3_REFIT=1                        ignores existing per-cell caches
all_combos <- function()
  expand.grid(joint = JOINT_ORDER, motion = MOTION_ORDER, dof = 1:3,
              stringsAsFactors = FALSE)[, c("joint", "motion", "dof")]

REFIT <- nzchar(Sys.getenv("V3_REFIT"))

# Log every K as it finishes, not just every cell. A cell is 20 fits and the
# slow ones run for tens of minutes, so without this a shard looks identical
# whether it is working hard or wedged. V3_QUIET_K=1 turns it off.
VERBOSE_K <- !nzchar(Sys.getenv("V3_QUIET_K"))

select_combos <- function() {
  cb <- all_combos()
  only <- Sys.getenv("V3_ONLY")
  if (nzchar(only)) {
    p <- strsplit(only, "|", fixed = TRUE)[[1]]
    if (length(p) != 3) stop("V3_ONLY must be '<joint>|<motion>|<dof>', got: ", only)
    k <- cb$joint == p[1] & cb$motion == p[2] & cb$dof == as.integer(p[3])
    if (!any(k)) stop("V3_ONLY matched no cell: ", only)
    return(cb[k, , drop = FALSE])
  }
  ns <- suppressWarnings(as.integer(Sys.getenv("V3_NSHARD", "1")))
  sh <- suppressWarnings(as.integer(Sys.getenv("V3_SHARD",  "0")))
  if (is.na(ns) || ns < 1) ns <- 1L
  if (is.na(sh) || sh < 0) sh <- 0L
  if (ns == 1L) return(cb)
  cb[(seq_len(nrow(cb)) %% ns) == (sh %% ns), , drop = FALSE]
}

# ---- data -------------------------------------------------------------------
load_long <- function() {
  if (!file.exists(LONG)) stop("Missing ", LONG, " — run prepare_monolix_data.py first.")
  D <- read.csv(LONG, stringsAsFactors = FALSE)
  D$cond <- factor(ifelse(D$in_vivo %in% c("True","TRUE","1","yes"), "in vivo", "ex vivo"),
                   levels = LEVELS)
  stopifnot(levels(D$cond)[1] == "in vivo")     # the reference must be in vivo
  D$study <- sub("_[^_]*$", "", D$ID)
  D[is.finite(D$TIME) & is.finite(D$Y), ]
}

# ---- the per-cell recipe ----------------------------------------------------
thin_by_unit <- function(df, cap) do.call(rbind, lapply(split(df, df$ID), function(u)
  if (nrow(u) <= cap) u else u[unique(round(seq(1, nrow(u), length.out = cap))), ]))

# The x-range SUPPORTED by the data: the widest interval spanned by at least
# MIN_UNITS_X shoulders.
#
# At MIN_UNITS_X = 1 — the current setting — every x any shoulder reaches counts,
# so this returns the absolute min/max and the trim keeps the whole abscissa.
# That is deliberate (Melanie's TODO), and it is a real trade-off rather than a
# free win: v2 used 2 precisely because a lone unit reaching far past the others
# drags the fitted curve across empty space — e.g. sternoclavicular horizontal
# flexion has one unit (17 of 13,064 rows) reaching 100 deg while the other nine
# stop near 0. Those stretches are now fitted, and the 95% band is what shows
# how little is holding them up.
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
  if (kdf < 2) return(NULL)
  q <- as.numeric(quantile(x, seq(0, 1, length.out = kdf + 1)[2:kdf]))
  q <- unique(round(q, 1)); q <- q[q > min(x) & q < max(x)]
  if (length(q) < 1) NULL else q
}

# Everything from "subset the plate" to "here are the rows that get fitted".
# Returns a one-row info frame in $row (populated even when the cell is skipped,
# so a skipped cell still shows WHY), and NULL $st when there is nothing to fit.
prepare_cell <- function(D, jt, mo, dof) {
  sub <- D[D$joint == jt & D$humeral_motion == mo & D$degree_of_freedom == dof, ]
  row <- data.frame(joint = jt, motion = mo, dof = dof, n_rows = nrow(sub),
                    n_units = 0L, n_ex = 0L, n_in = 0L, st_ex = 0L, st_in = 0L,
                    mode = "skip", x_lo = NA_real_, x_hi = NA_real_,
                    x_ex_lo = NA_real_, x_ex_hi = NA_real_,
                    x_in_lo = NA_real_, x_in_hi = NA_real_,
                    n_fitted = 0L, note = "", stringsAsFactors = FALSE)
  none <- function(why) { row$note <- why; list(row = row, st = NULL) }

  if (nrow(sub) < MIN_ROWS) return(none("too few rows"))
  sub$ID <- factor(sub$ID); sub <- sub[order(sub$ID, sub$TIME), ]
  units_by <- tapply(as.character(sub$cond), sub$ID, function(z) z[1])
  row$n_units <- length(units_by)
  row$n_ex <- sum(units_by == "ex vivo"); row$n_in <- sum(units_by == "in vivo")
  row$st_ex <- length(unique(sub$study[sub$cond == "ex vivo"]))
  row$st_in <- length(unique(sub$study[sub$cond == "in vivo"]))
  # How far each condition reaches ON ITS OWN, before any trimming. x_lo/x_hi
  # below is the OVERLAP, which is what the model is fitted on; these two spans
  # are what the overlap is an intersection OF. Computed above the early return
  # so skipped and single-condition cells carry them too — that is what lets a
  # "not compared" row still show why (one condition covering ground the other
  # never reaches). Safe on `sub` rather than the thinned rows: thin_by_unit
  # keeps each unit's first and last row, so per-unit ranges are unchanged.
  for (cc in c("ex vivo", "in vivo")) {
    k <- sub$cond == cc
    r <- if (any(k)) support_range(sub$TIME[k], sub$ID[k]) else c(NA_real_, NA_real_)
    if (cc == "ex vivo") { row$x_ex_lo <- r[1]; row$x_ex_hi <- r[2] }
    else                 { row$x_in_lo <- r[1]; row$x_in_hi <- r[2] }
  }
  compare <- row$n_ex >= MIN_PER_COND && row$n_in >= MIN_PER_COND
  if (!compare && row$n_units < MIN_TOTAL)
    return(none(sprintf("too few units (%d ex / %d in)", row$n_ex, row$n_in)))

  st <- thin_by_unit(sub, CAP); st <- st[order(st$ID, st$TIME), ]
  full_rng <- range(st$TIME)                    # boundary knots from the UNTRIMMED range
  ov <- if (compare) overlap_range(st) else support_range(st$TIME, st$ID)
  if (TRIM_TO_OVERLAP || !compare)
    st <- st[st$TIME >= ov[1] & st$TIME <= ov[2], ]
  else ov <- full_rng
  st$ID <- droplevels(st$ID); st$IDnum <- as.integer(st$ID)
  st$cond <- factor(as.character(st$cond), levels = LEVELS)
  if (nrow(st) < MIN_ROWS || nlevels(st$ID) < MIN_TOTAL) return(none("empty after trim"))

  row$x_lo <- ov[1]; row$x_hi <- ov[2]
  row$n_fitted <- nrow(st); row$mode <- if (compare) "compare" else "single"
  list(row = row, st = st, ov = ov, full_rng = full_rng, compare = compare)
}

# The basis, and the only place it is built. Interior knots at the cell's own
# quantiles (Melanie's step 1b); boundary knots at the FITTED range.
#
# v2 put the boundary knots at the untrimmed range instead, on the grounds that
# tying them to the overlap made the residual tails worse. That was measured at
# a fixed df = 4, and it does not survive per-cell K. A natural spline is linear
# BEYOND its boundary knots, so when the fitted data occupies only part of the
# basis's support the outer columns are near-collinear over the observed range;
# their random-effect variances then blow up and drag the fixed-effect
# covariance with them. Measured on scapulothoracic / frontal / DoF2 at K = 7,
# same 3378 rows, boundary (-3.5, 187) against (14, 150):
#
#   random-effect params   63 40 46 39 43 59 2e+02 5e+02   vs   48 13 20 23 35 46 58 1.1e+02
#   95% CI half-width at the left edge      7.1 deg        vs   4.6 deg
#
# The difference curve is the reporting object, so a 1.6x inflated edge band is
# not a detail. Anchoring to the fitted range also makes this match iteration 06
# exactly. The overlap trim, not the boundary knots, is what keeps the curve off
# empty space.
cell_basis <- function(st, K, rng) {
  kn <- knots_for(st$TIME, K)
  B  <- try(ns(st$TIME, knots = kn, Boundary.knots = rng), silent = TRUE)
  if (inherits(B, "try-error")) return(NULL)
  nm <- paste0("ns", seq_len(ncol(B)))
  for (j in seq_along(nm)) st[[nm[j]]] <- B[, j]
  list(df = st, B = B, nm = nm, knots = kn, Kc = ncol(B))
}

f_red  <- function(nm) as.formula(paste("Y ~", paste(nm, collapse = " + ")))
f_full <- function(nm) as.formula(paste("Y ~", paste(nm, collapse = " + "), "+ cond +",
                                        paste(paste0("cond:", nm), collapse = " + ")))
f_rand <- function(nm) as.formula(paste("~", paste(nm, collapse = " + ")))

# hlme converged AND produced usable numbers. m$conv == 1 alone is not enough:
# a run can report convergence with non-finite entries in best.
fit_ok <- function(m)
  !inherits(m, "try-error") && !is.null(m$conv) && m$conv == 1 && all(is.finite(m$best))

# per-shoulder fit quality from a fitted model's own prediction rows
per_shoulder <- function(m, dat) {
  pr <- m$pred
  do.call(rbind, lapply(split(seq_len(nrow(dat)), droplevels(dat$ID)), function(i) {
    if (!length(i)) return(NULL)
    data.frame(shoulder = as.character(dat$ID[i[1]]),
               condition = as.character(dat$cond[i[1]]), n = length(i),
               rmse = sqrt(mean(pr$resid_ss[i]^2)), max_abs = max(abs(pr$resid_ss[i])),
               r2_marg = r2_conc(dat$Y[i], pr$pred_m[i]),
               r2_ss   = r2_conc(dat$Y[i], pr$pred_ss[i]), row.names = NULL)
  }))
}

# The prediction rows figures 04-06 need, kept instead of the hlme object — which
# is far larger than the parts actually used.
pred_rows <- function(m, dat)
  data.frame(ID = as.character(dat$ID), cond = as.character(dat$cond), TIME = dat$TIME,
             obs = dat$Y, pred_m = m$pred$pred_m, pred_ss = m$pred$pred_ss,
             resid_ss = m$pred$resid_ss, row.names = NULL)

# RMSE is right-skewed, so a mean +/- sd rule would be dragged by the very
# outliers it is meant to catch. MAD-based z on the log scale instead.
flag_shoulders <- function(ind) {
  ind$ratio_to_median <- ind$rmse / median(ind$rmse)
  lg <- log(ind$rmse); m <- mad(lg, constant = 1)
  ind$robust_z <- if (m > 0) (lg - median(lg)) / (1.4826 * m) else 0
  ind$flagged  <- ind$robust_z > 3
  ind[order(-ind$rmse), ]
}

# ---- gather the per-cell RDS files into the plate ---------------------------
# One master row per cell, and one entry per cell that reached a K. Step 3 calls
# this to build the cache; the figures call it when there is no cache yet, which
# is what lets a HALF-FINISHED sweep be drawn — a cell with only step 1 carries
# its BIC-vs-K table and its individual fit, and its step-2 fields stay NULL
# until step 2 runs. Keeping it here rather than inside step 3 means the two
# cannot drift.
collect_cells <- function() {
  cb <- all_combos(); rows <- list(); CUR <- list(); cfl <- list()
  for (i in seq_len(nrow(cb))) {
    jt <- cb$joint[i]; mo <- cb$motion[i]; dof <- cb$dof[i]
    s1 <- cell_load(1, jt, mo, dof); if (is.null(s1)) next
    s2 <- cell_load(2, jt, mo, dof); s3 <- cell_load(3, jt, mo, dof)

    r <- s1$row
    r$K_star      <- s1$K_star
    r$k_basis     <- if (is.null(s1$k_basis)) NA_character_ else s1$k_basis
    r$bic_edge    <- if (is.null(s1$bic_edge)) NA else s1$bic_edge
    r$n_converged <- if (is.null(s1$n_converged)) NA_integer_ else s1$n_converged
    r$dBIC_next   <- if (is.null(s1$dBIC_next)) NA_real_ else s1$dBIC_next
    r$knots       <- if (is.null(s1$knots)) NA_character_ else paste(s1$knots, collapse = "/")
    r$rmse_max    <- if (is.null(s1$ind)) NA_real_ else max(s1$ind$rmse)
    r$rmse_median <- if (is.null(s1$ind)) NA_real_ else median(s1$ind$rmse)
    r$rmse_ratio  <- if (is.null(s1$ind)) NA_real_ else max(s1$ind$ratio_to_median)
    r$n_flagged   <- if (is.null(s1$ind)) NA_integer_ else sum(s1$ind$flagged)
    r$n_shoulders_ref <- if (is.null(s1$ind)) NA_integer_ else nrow(s1$ind)
    # The random part mirrors K, so the fit estimates K+1 diagonal variances from
    # n_shoulders_ref subjects. When that ratio reaches 1 the model is asking for
    # as many variance parameters as it has curves to estimate them from, and the
    # retained K stops being a statement about the trajectory. lcmm penalises BIC
    # on the SUBJECT count, so in a thin cell an extra parameter costs only
    # 2*log(n) and BIC barely resists — it is not a safeguard here. Recorded so
    # the flag is data, not something a reader has to notice on a plot.
    r$var_per_shoulder <- if (is.na(r$K_star) || is.na(r$n_shoulders_ref)) NA_real_
                          else (r$K_star + 1) / r$n_shoulders_ref
    r$overparam <- !is.na(r$var_per_shoulder) & r$var_per_shoulder >= 1
    r$r2_marg   <- if (is.null(s1$r2_marg)) NA_real_ else s1$r2_marg
    r$r2_ss     <- if (is.null(s1$r2_ss)) NA_real_ else s1$r2_ss

    for (nmc in c("diff_mean_signed","diff_mean_abs","diff_max_abs","diff_lo_curve",
                  "diff_hi_curve","sig_frac_x","diff_mean_lo","diff_mean_hi",
                  "resid_sd","skew","kurtosis","shapiro_p"))
      r[[nmc]] <- if (!is.null(s2) && !is.null(s2[[nmc]])) s2[[nmc]] else NA_real_
    r$h1_converged <- if (is.null(s2)) NA else isTRUE(s2$converged)

    r$loglik_reduced <- if (!is.null(s3) && isTRUE(s3$converged)) s3$loglik_reduced else NA_real_
    r$loglik_full    <- if (!is.null(s3) && isTRUE(s3$converged)) s3$loglik_full else NA_real_
    r$LR             <- if (!is.null(s3) && isTRUE(s3$converged)) s3$LR else NA_real_
    r$lrt_df         <- if (!is.null(s3) && isTRUE(s3$converged)) s3$df else NA_integer_
    r$p_lrt          <- if (!is.null(s3) && isTRUE(s3$converged)) s3$p_lrt else NA_real_
    r$seconds        <- sum(c(s1$seconds, s2$seconds, s3$seconds), na.rm = TRUE)
    # which stages this cell actually has, so a figure can skip what is absent
    r$has_step2      <- !is.null(s2) && isTRUE(s2$converged)
    r$has_step3      <- !is.null(s3) && isTRUE(s3$converged)

    rows[[length(rows)+1]] <- r
    if (!is.na(s1$K_star))
      CUR[[cell_key(jt, mo, dof)]] <- list(
        # step 1 — always present once a K was reached
        bic = s1$bic, K_star = s1$K_star, bic_edge = s1$bic_edge,
        dBIC_next = s1$dBIC_next, n_converged = s1$n_converged,
        k_basis = s1$k_basis, ind = s1$ind, pred_ref = s1$pred_ref,
        r2_marg = s1$r2_marg, r2_ss = s1$r2_ss,
        K = if (!is.null(s2)) s2$K else s1$K_star,
        compare = identical(s1$row$mode, "compare"),
        # step 2 — NULL until it has run
        curves = s2$curves, dif = s2$dif, raw = s2$raw, knots = s2$knots,
        bf = s2$bf, Vf = s2$Vf, full_rng = s2$full_rng, Kc = s2$Kc,
        ov = s2$ov, pred_full = s2$pred_full)
    idc <- data.frame(joint = jt, motion = mo, dof = dof)
    if (!is.null(s2) && !is.null(s2$coefs)) cfl[[length(cfl)+1]] <- cbind(idc, K = s2$K, s2$coefs)
    if (!is.null(s3) && !is.null(s3$coefs0)) cfl[[length(cfl)+1]] <- cbind(idc, K = s3$K, s3$coefs0)
  }
  master <- do.call(rbind, rows); rownames(master) <- NULL

  # Melanie: "ajustement multiple a faire ici, sur la planche au complet. Prendre
  # bonferroni." Bonferroni is the headline; BH stays beside it so v3 can be read
  # against v2, which used BH.
  master$p_bonf <- NA_real_; master$p_bh <- NA_real_
  cmp <- master$mode == "compare" & !is.na(master$p_lrt)
  if (any(cmp)) {
    master$p_bonf[cmp] <- p.adjust(master$p_lrt[cmp], "bonferroni")
    master$p_bh[cmp]   <- p.adjust(master$p_lrt[cmp], "BH")
  }
  list(cells = CUR, master = master,
       coefs = if (length(cfl)) do.call(rbind, cfl) else NULL)
}

# ---- the x-intervals where a difference band excludes zero ------------------
# The solid bars of every significance map. A run-length encoding of the
# pointwise flag, so a cell whose significance is broken into several stretches
# draws as several bars rather than one bar spanning the gaps.
sig_segments <- function(e) {
  if (is.null(e) || is.null(e$dif)) return(NULL)
  z <- e$dif; sig <- (z$lo > 0) | (z$hi < 0)
  r <- rle(sig); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  k <- which(r$values)
  if (!length(k)) return(data.frame(lo = numeric(0), hi = numeric(0)))
  data.frame(lo = z$x[starts[k]], hi = z$x[ends[k]])
}
