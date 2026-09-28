# =============================================================================
# GENERALIZED SWEEP v4 — shared configuration, helpers and the per-cell recipe
#
# Sourced by all four halves: analyze_all_v4_step{1,2,3}.R and
# analyze_all_v4_figures.R. Everything the steps have to agree on lives here —
# where the files are, what defines a cell, how a cell's rows and basis are
# built, the palette and the panel ordering.
#
# THE PREPARATION IS SHARED ON PURPOSE. Step 1 fits the reduced model on every
# row at K = 1..10, step 2 fits H1 at the selected K on everything, step 3 takes
# H0 at that K (reused from step 1) and differences the log-likelihoods. That
# difference is meaningless unless every step built IDENTICAL rows and an
# IDENTICAL basis, so prepare_cell() and cell_basis() are the single source.
#
# K IS A DECISION, NOT ONLY A COMPUTATION. Step 1 proposes one per cell and
# writes it to K_selection.csv; step 2 reads selected_K from that file, so a K
# changed by hand after looking at 01_choosing_k is what the rest of the
# pipeline uses.
#
# lcmm is deliberately NOT loaded here — only the fitting steps need it, and the
# point of the split is that redrawing a figure costs nothing.
# =============================================================================

suppressMessages({ library(splines) })

LONG      <- "spartacus_angles_long.csv"
OUT       <- file.path("figures", "generalized_v4")
CACHE_DIR <- file.path(OUT, "cache")
STEP_DIR  <- function(s) file.path(CACHE_DIR, paste0("step", s))
CACHE     <- file.path(CACHE_DIR, "v4_fit.rds")
CELL_DIR  <- file.path(OUT, "cells")
# The hand-editable K decisions. Deliberately outside cache/, and tracked by git.
SEL_CSV   <- file.path(OUT, "K_selection.csv")

# ---- model constants --------------------------------------------------------
# The K grid step 1 fits. Capped at 10: past that the extra df describe how far
# each shoulder's own curve can bend rather than the trajectory, and thin cells
# ended up asking for more random variances than they had shoulders.
KS           <- 1:10
# BIC gain thresholds drawn on 01_choosing_k and used for K_5pct / K_1pct:
# rel(K) = (BIC(K-1) - BIC(K)) / (BIC(first K) - min BIC), the gain of step K
# as a fraction of the total BIC drop over the grid.
REL_THRESH   <- c(pct5 = 0.05, pct1 = 0.01)
# V4_KS='1:3' truncates the grid for a wiring smoke test. It is part of the
# fingerprint, so a cache built under a short grid reads as stale rather than
# being mistaken for a real selection.
if (nzchar(Sys.getenv("V4_KS")))
  KS <- eval(parse(text = Sys.getenv("V4_KS")))
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
# slower: V4_NPROC=4 with V4_NSHARD=5 keeps 20 of 22 cores busy and no more.
# Not part of fit_params(): it changes how long the sweep takes, not what comes
# out of it.
NPROC        <- max(1L, min(8L, parallel::detectCores() - 1L))
if (nzchar(Sys.getenv("V4_NPROC")))
  NPROC <- max(1L, as.integer(Sys.getenv("V4_NPROC")))

COL   <- c("ex vivo" = "#d95f02", "in vivo" = "#1b9e77")
BLUE  <- "#377eb8"
# The single curve of a cell that is not compared: neither condition's colour,
# since the model behind it has no condition term.
COL_POOLED <- "grey35"
cond_col <- function(cc) unname(ifelse(cc %in% names(COL), COL[cc], COL_POOLED))

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
# V4_NSHARD=4 V4_SHARD=0..3 each write their own per-cell files into the shared
# cache dir — one writer per file, so no locking is needed.
#   V4_ONLY='<joint>|<motion>|<dof>'  runs exactly one cell (verification)
#   V4_REFIT=1                        ignores existing per-cell caches
all_combos <- function()
  expand.grid(joint = JOINT_ORDER, motion = MOTION_ORDER, dof = 1:3,
              stringsAsFactors = FALSE)[, c("joint", "motion", "dof")]

REFIT <- nzchar(Sys.getenv("V4_REFIT"))

# Log every K as it finishes, not just every cell. A cell is up to 10 fits and
# the slow ones run for minutes, so without this a shard looks identical
# whether it is working hard or wedged. V4_QUIET_K=1 turns it off.
VERBOSE_K <- !nzchar(Sys.getenv("V4_QUIET_K"))

select_combos <- function() {
  cb <- all_combos()
  only <- Sys.getenv("V4_ONLY")
  if (nzchar(only)) {
    p <- strsplit(only, "|", fixed = TRUE)[[1]]
    if (length(p) != 3) stop("V4_ONLY must be '<joint>|<motion>|<dof>', got: ", only)
    k <- cb$joint == p[1] & cb$motion == p[2] & cb$dof == as.integer(p[3])
    if (!any(k)) stop("V4_ONLY matched no cell: ", only)
    return(cb[k, , drop = FALSE])
  }
  ns <- suppressWarnings(as.integer(Sys.getenv("V4_NSHARD", "1")))
  sh <- suppressWarnings(as.integer(Sys.getenv("V4_SHARD",  "0")))
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

# ---- K selection: candidate K, and the hand-editable CSV ---------------------
nz <- function(x, d = NA) if (is.null(x) || !length(x)) d else x
fmtK <- function(k) if (is.null(k) || !length(k) || is.na(k)) " —" else sprintf("%2d", as.integer(k))
isTRUE_vec <- function(x) {
  if (is.null(x)) return(logical(0))
  x <- as.logical(x); !is.na(x) & x
}

# Candidate K from one cell's BIC table — recomputed from the cached table every
# time it is needed, never frozen at fit time, so changing REL_THRESH or this
# rule costs no refit. rel(K) is the BIC gain from the previous CONVERGED K
# (so one K that failed does not break the chain) as a fraction of the total
# BIC drop over the grid. Unlike a ratio over |BIC(K-1)|, it stays defined when
# the BIC sits near zero or turns negative, which is common in the thin cells.
# K_5pct / K_1pct are the first K whose gain falls below 5 % / 1 % of that drop;
# the default K is K_5pct, or the BIC minimum when there is no drop at all.
k_thresholds <- function(bic) {
  out <- list(K_bic = NA_integer_, K_5pct = NA_integer_, K_1pct = NA_integer_,
              default_K = NA_integer_, dBIC_next = NA_real_)
  if (is.null(bic)) return(out)
  ok <- bic[bic$converged & is.finite(bic$BIC), , drop = FALSE]
  ok <- ok[order(ok$K), , drop = FALSE]
  if (!nrow(ok)) return(out)
  out$K_bic <- as.integer(ok$K[which.min(ok$BIC)])
  if (nrow(ok) > 1) {
    out$dBIC_next <- sort(ok$BIC)[2] - min(ok$BIC)
    total <- ok$BIC[1] - min(ok$BIC)
    rel <- if (total > 0) (head(ok$BIC, -1) - ok$BIC[-1]) / total
           else rep(NA_real_, nrow(ok) - 1)
    Kr  <- ok$K[-1]
    first_below <- function(t) {
      w <- which(rel < t)
      if (length(w)) as.integer(Kr[w[1]]) else NA_integer_
    }
    out$K_5pct <- first_below(REL_THRESH[["pct5"]])
    out$K_1pct <- first_below(REL_THRESH[["pct1"]])
  }
  out$default_K <- if (!is.na(out$K_5pct)) out$K_5pct else out$K_bic
  out
}

# Write-then-rename, so a reader — another shard, or a person with the file open
# — never sees a half-written file. Sharded step-1 runs all write these tables.
write_atomic <- function(df, path) {
  tmp <- paste0(path, ".tmp", Sys.getpid())
  write.csv(df, tmp, row.names = FALSE, na = "")
  invisible(file.rename(tmp, path))
}

# The CSV as edited by hand. A spreadsheet saved in a French locale uses ';' as
# the separator, so that is detected rather than failing on "missing columns".
# selected_K must be a whole number: "8.5" is refused, not truncated.
read_selection <- function() {
  if (!file.exists(SEL_CSV)) return(NULL)
  first <- readLines(SEL_CSV, n = 1, warn = FALSE)
  rd <- if (grepl(";", first) && !grepl(",", first)) utils::read.csv2 else utils::read.csv
  s <- rd(SEL_CSV, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  need <- c("joint", "motion", "dof", "selected_K")
  if (!all(need %in% names(s)))
    stop(SEL_CSV, " is missing column(s): ", paste(setdiff(need, names(s)), collapse = ", "))
  k <- suppressWarnings(as.numeric(s$selected_K))
  bad <- !is.na(s$selected_K) & (is.na(k) | k != round(k))
  if (any(bad))
    stop("selected_K must be a whole number — check line(s) ",
         paste(which(bad) + 1, collapse = ", "), " of ", SEL_CSV)
  s$selected_K <- as.integer(k)
  if (!"note" %in% names(s)) s$note <- NA_character_
  s
}

selected_K_of <- function(sel, jt, mo, dof) {
  if (is.null(sel)) return(NA_integer_)
  i <- which(sel$joint == jt & sel$motion == mo & sel$dof == dof)
  if (length(i) > 1)
    stop(SEL_CSV, " has ", length(i), " rows for ", cell_key(jt, mo, dof), " — keep one")
  if (!length(i)) NA_integer_ else sel$selected_K[i]
}

# Rebuild K_selection.csv from the step-1 caches WITHOUT losing a decision.
# Computed columns are refreshed; a selected_K or note already in the file is
# kept; a row for a cell with no cache right now (cleared, still running) is
# kept too. An empty selected_K is filled with the default, so blanking a cell
# is how it goes back to the default after the data changed. `old` can be passed
# in — e.g. with every selected_K blanked — to reset the defaults in ONE atomic
# write, with no moment where the file on disk holds empty K.
write_selection <- function(old = read_selection()) {
  cb <- all_combos(); rows <- list()
  for (i in seq_len(nrow(cb))) {
    jt <- cb$joint[i]; mo <- cb$motion[i]; dof <- cb$dof[i]
    o <- cell_load(1, jt, mo, dof); if (is.null(o)) next
    r <- o$row; th <- k_thresholds(o$bic)
    rows[[length(rows)+1]] <- data.frame(
      joint = jt, motion = mo, dof = dof, n_in = r$n_in, n_ex = r$n_ex, mode = r$mode,
      n_converged = nz(o$n_converged, 0L),
      K_bic = th$K_bic, K_5pct = th$K_5pct, K_1pct = th$K_1pct,
      default_K = th$default_K, selected_K = NA_integer_,
      note = if (nzchar(nz(r$note, ""))) r$note else NA_character_,
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(invisible(NULL))
  new <- do.call(rbind, rows)
  key <- function(d) paste(d$joint, d$motion, d$dof, sep = "|")
  if (!is.null(old)) {
    m <- match(key(new), key(old))
    keepK <- !is.na(m) & !is.na(old$selected_K[m])
    new$selected_K[keepK] <- old$selected_K[m[keepK]]
    keepN <- !is.na(m) & !is.na(old$note[m])
    new$note[keepN] <- old$note[m[keepN]]
    orphan <- old[!key(old) %in% key(new), intersect(names(new), names(old)), drop = FALSE]
    if (nrow(orphan)) {
      for (nm in setdiff(names(new), names(orphan))) orphan[[nm]] <- NA
      new <- rbind(new, orphan[, names(new)])
    }
  }
  fill <- is.na(new$selected_K)
  new$selected_K[fill] <- new$default_K[fill]
  new <- new[order(match(new$joint, JOINT_ORDER), match(new$motion, MOTION_ORDER), new$dof), ]
  write_atomic(new, SEL_CSV)
  invisible(new)
}

# Everything step 1 writes to disk besides its caches: the selection CSV, the
# BIC-by-K table and the individual fit at each cell's selected K. Called at the
# end of every step-1 shard, once more by run_v4_sweep.sh after all shards, and
# by step 3, so the individual-fit table follows a K changed by hand.
write_step1_tables <- function() {
  write_selection()
  sel <- read_selection()
  cb <- all_combos(); bl <- list(); il <- list()
  for (i in seq_len(nrow(cb))) {
    jt <- cb$joint[i]; mo <- cb$motion[i]; dof <- cb$dof[i]
    o <- cell_load(1, jt, mo, dof); if (is.null(o) || is.null(o$bic)) next
    Ks <- selected_K_of(sel, jt, mo, dof); th <- k_thresholds(o$bic)
    idc <- data.frame(joint = jt, motion = mo, dof = dof)
    bl[[length(bl)+1]] <- cbind(idc, K_bic = th$K_bic, K_5pct = th$K_5pct,
                                K_1pct = th$K_1pct, selected_K = Ks, o$bic)
    pk <- if (is.na(Ks)) NULL else o$perK[[as.character(Ks)]]
    if (!is.null(pk)) il[[length(il)+1]] <- cbind(idc, K = Ks, flag_shoulders(pk$ps))
  }
  if (length(bl)) write_atomic(do.call(rbind, bl), file.path(OUT, "00_bic_by_k.csv"))
  if (length(il)) write_atomic(do.call(rbind, il), file.path(OUT, "00_individual_fit.csv"))
  invisible(sel)
}

# ---- gather the per-cell RDS files into the plate ---------------------------
# One master row per cell, and one entry per cell with a usable K. Step 3 calls
# this to build the cache; the figures call it when there is no cache or it is
# stale, which is what lets a HALF-FINISHED sweep be drawn. Keeping it here
# rather than inside step 3 means the two cannot drift.
#
# The K in use is the CSV's selected_K (step 1's default for a cell missing from
# the CSV). Step-2 and step-3 results only count when they were fitted at that
# K: after a hand edit they describe another model, so they are set aside and
# the cell is marked refit_pending rather than drawn as if current.
collect_cells <- function() {
  sel <- read_selection()
  cb <- all_combos(); rows <- list(); CUR <- list(); cfl <- list()
  for (i in seq_len(nrow(cb))) {
    jt <- cb$joint[i]; mo <- cb$motion[i]; dof <- cb$dof[i]
    s1 <- cell_load(1, jt, mo, dof); if (is.null(s1)) next
    s2 <- cell_load(2, jt, mo, dof); s3 <- cell_load(3, jt, mo, dof)

    th <- k_thresholds(s1$bic)
    Ks <- selected_K_of(sel, jt, mo, dof)
    if (is.na(Ks)) Ks <- th$default_K
    at_K <- function(s)
      if (!is.null(s) && !is.na(Ks) && identical(as.integer(s$K), as.integer(Ks))) s else NULL
    s2k <- at_K(s2); s3k <- at_K(s3)
    pk  <- if (is.na(Ks)) NULL else s1$perK[[as.character(Ks)]]
    ind <- if (is.null(pk)) NULL else flag_shoulders(pk$ps)

    r <- s1$row
    r$K_star        <- Ks
    r$K_bic         <- th$K_bic
    r$K_5pct        <- th$K_5pct
    r$K_1pct        <- th$K_1pct
    r$default_K     <- th$default_K
    r$K_manual      <- !is.na(Ks) && !is.na(r$default_K) && Ks != r$default_K
    r$k_basis       <- nz(s1$k_basis, NA_character_)
    r$bic_edge      <- !is.na(th$K_bic) && th$K_bic %in% range(KS)
    r$n_converged   <- nz(s1$n_converged, NA_integer_)
    r$dBIC_next     <- th$dBIC_next
    r$knots         <- if (is.null(pk)) NA_character_ else paste(pk$knots, collapse = "/")
    r$rmse_max      <- if (is.null(ind)) NA_real_ else max(ind$rmse)
    r$rmse_median   <- if (is.null(ind)) NA_real_ else median(ind$rmse)
    r$rmse_ratio    <- if (is.null(ind)) NA_real_ else max(ind$ratio_to_median)
    r$n_flagged     <- if (is.null(ind)) NA_integer_ else sum(ind$flagged)
    r$n_shoulders_ref <- if (is.null(ind)) NA_integer_ else nrow(ind)
    # The random part mirrors K, so the fit estimates K+1 diagonal variances from
    # n_shoulders_ref subjects. When that ratio reaches 1 the model is asking for
    # as many variance parameters as it has curves to estimate them from, and the
    # K in use stops being a statement about the trajectory. lcmm penalises BIC
    # on the SUBJECT count, so in a thin cell an extra parameter costs only
    # 2*log(n) and BIC barely resists — it is not a safeguard here.
    r$var_per_shoulder <- if (is.na(r$K_star) || is.na(r$n_shoulders_ref)) NA_real_
                          else (r$K_star + 1) / r$n_shoulders_ref
    r$overparam <- !is.na(r$var_per_shoulder) & r$var_per_shoulder >= 1
    r$r2_marg   <- nz(pk$r2_marg, NA_real_)
    r$r2_ss     <- nz(pk$r2_ss, NA_real_)

    for (nmc in c("diff_mean_signed","diff_mean_abs","diff_max_abs","diff_lo_curve",
                  "diff_hi_curve","sig_frac_x","diff_mean_lo","diff_mean_hi",
                  "resid_sd","skew","kurtosis","shapiro_p"))
      r[[nmc]] <- if (!is.null(s2k) && !is.null(s2k[[nmc]])) s2k[[nmc]] else NA_real_
    r$h1_converged <- if (is.null(s2k)) NA else isTRUE(s2k$converged)

    r$loglik_reduced <- if (!is.null(s3k) && isTRUE(s3k$converged)) s3k$loglik_reduced else NA_real_
    r$loglik_full    <- if (!is.null(s3k) && isTRUE(s3k$converged)) s3k$loglik_full else NA_real_
    r$LR             <- if (!is.null(s3k) && isTRUE(s3k$converged)) s3k$LR else NA_real_
    r$lrt_df         <- if (!is.null(s3k) && isTRUE(s3k$converged)) s3k$df else NA_integer_
    r$p_lrt          <- if (!is.null(s3k) && isTRUE(s3k$converged)) s3k$p_lrt else NA_real_
    r$h0_source      <- if (!is.null(s3k)) nz(s3k$h0_source, NA_character_) else NA_character_
    r$seconds        <- sum(c(s1$seconds, s2k$seconds, s3k$seconds), na.rm = TRUE)
    # which stages this cell actually has at its K, so a figure can skip what is absent
    r$has_step2      <- !is.null(s2k) && isTRUE(s2k$converged)
    r$has_step3      <- !is.null(s3k) && isTRUE(s3k$converged)
    r$refit_pending  <- (!is.null(s2) && is.null(s2k)) ||
                        (!is.null(s3) && is.null(s3k) && !is.null(s2k))

    rows[[length(rows)+1]] <- r
    if (!is.null(pk))
      CUR[[cell_key(jt, mo, dof)]] <- list(
        # step 1 — present once the K in use has a converged fit
        bic = s1$bic, K_star = Ks, K_bic = r$K_bic, K_5pct = r$K_5pct,
        K_1pct = r$K_1pct, default_K = r$default_K, K_manual = r$K_manual,
        bic_edge = r$bic_edge,
        dBIC_next = r$dBIC_next, n_converged = r$n_converged,
        k_basis = r$k_basis, ind = ind, pred_ref = pk$pred,
        r2_marg = r$r2_marg, r2_ss = r$r2_ss, K = Ks,
        compare = identical(s1$row$mode, "compare"),
        # step 2 at the K in use — NULL until it has run
        curves = s2k$curves, dif = s2k$dif, raw = s2k$raw, knots = s2k$knots,
        bf = s2k$bf, Vf = s2k$Vf, full_rng = s2k$full_rng, Kc = s2k$Kc,
        ov = s2k$ov, pred_full = s2k$pred_full)
    idc <- data.frame(joint = jt, motion = mo, dof = dof)
    if (!is.null(s2k$coefs)) cfl[[length(cfl)+1]] <- cbind(idc, K = s2k$K, s2k$coefs)
    if (!is.null(s3k$coefs0)) cfl[[length(cfl)+1]] <- cbind(idc, K = s3k$K, s3k$coefs0)
  }
  master <- do.call(rbind, rows); rownames(master) <- NULL

  # Melanie: "ajustement multiple a faire ici, sur la planche au complet. Prendre
  # bonferroni." Bonferroni is the headline; BH stays beside it for comparison
  # with the earlier sweeps.
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
