# =============================================================================
# GENERALIZED SWEEP v2 — shared configuration and helpers
#
# Sourced by BOTH analyze_all_v2_fit.R and analyze_all_v2_figures.R. Everything
# here is what the two halves have to agree on: where the files live, the model
# constants that define a cell, the palette and the ordering of the panels.
#
# lcmm is deliberately NOT loaded here — only the fit needs it, and the whole
# point of the split is that redrawing a figure costs nothing.
# =============================================================================

suppressMessages({ library(splines); library(dplyr) })

LONG      <- "spartacus_angles_long.csv"
OUT       <- file.path("figures", "generalized_v2")
CACHE_DIR <- file.path(OUT, "cache")
CACHE     <- file.path(CACHE_DIR, "v2_fit.rds")

KDF          <- 4      # ns() columns; 3 interior knots
CAP          <- 100    # max points per shoulder
MIN_PER_COND <- 3      # >=3 shoulders in BOTH conditions -> compare
MIN_TOTAL    <- 4      # >=4 shoulders total              -> at least single
WALD_CRIT    <- 1.96
MIN_UNITS_X  <- 2      # a curve is only drawn where >=2 shoulders have data

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

dir.create(OUT,       showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

slug  <- function(s) gsub("(^_|_$)", "", gsub("[^a-z0-9]+", "_", tolower(s)))
capj  <- function(j) paste0(toupper(substring(j, 1, 1)), substring(j, 2))
stars <- function(p) if (is.na(p)) "" else if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else "ns"
# How a cell is addressed in the cache. The fit writes these keys, the figures
# read them; keeping the idiom in one place is what stops the two halves drifting.
cell_key <- function(jt, mo, dof) paste(jt, mo, dof, sep = "|")

# Render every figure TWICE: a high-resolution PNG and a vector PDF.
# The draw argument is a FUNCTION, not a braced expression, so it can be called
# once per device — a promise would only evaluate on the first device and the
# second file would come out blank.
FIG_SCALE <- 2.2          # pixel multiplier; res scales with it so text stays proportional
figout <- function(base, draw, w = 1100, h = 750, res = 130) {
  png(file.path(OUT, paste0(base, ".png")),
      width = round(w*FIG_SCALE), height = round(h*FIG_SCALE), res = round(res*FIG_SCALE))
  op <- par(no.readonly = TRUE); tryCatch(draw(), finally = { par(op); dev.off() })
  # cairo_pdf, NOT pdf(): the legacy device uses Type 1 fonts with no glyph for
  # the circles/arrows/Greek used here, so those characters silently vanish.
  if (capabilities("cairo")) cairo_pdf(file.path(OUT, paste0(base, ".pdf")),
                                       width = w/res, height = h/res)
  else pdf(file.path(OUT, paste0(base, ".pdf")), width = w/res, height = h/res)
  op <- par(no.readonly = TRUE); tryCatch(draw(), finally = { par(op); dev.off() })
}

# ---- cache metadata ---------------------------------------------------------
# The constants that actually change what a fit means. NPROC is not among them:
# it changes how long the sweep takes, not what comes out of it.
fit_params <- function()
  list(KDF = KDF, CAP = CAP, MIN_PER_COND = MIN_PER_COND, MIN_TOTAL = MIN_TOTAL,
       WALD_CRIT = WALD_CRIT, MIN_UNITS_X = MIN_UNITS_X)

fit_meta <- function() {
  fi <- file.info(LONG)
  list(params = fit_params(),
       input = list(file = LONG, size = fi$size, mtime = fi$mtime),
       fitted_at = Sys.time(),
       R_version = getRversion(),
       lcmm_version = tryCatch(as.character(utils::packageVersion("lcmm")),
                               error = function(e) NA_character_))
}

# Once fitting and drawing are separate commands, nothing stops you editing KDF
# and redrawing from a cache fitted under the old value. Warn — loudly, by name —
# rather than stop: a deliberate mismatch is sometimes what you want.
check_meta <- function(meta) {
  now <- fit_params(); was <- meta$params
  if (!length(was)) {
    warning("cache carries no settings fingerprint; cannot check it against the current ",
            "constants. Refit with  Rscript analyze_all_v2_fit.R  if in doubt.", call. = FALSE)
    return(invisible(NULL))
  }
  bad <- names(now)[!mapply(identical, now[names(now)], was[names(now)])]
  if (length(bad))
    warning("cache was fitted with different settings: ",
            paste(sprintf("%s = %s (cache) vs %s (now)", bad,
                          sapply(bad, function(k) format(was[[k]])),
                          sapply(bad, function(k) format(now[[k]]))), collapse = "; "),
            "\n  re-run  Rscript analyze_all_v2_fit.R  to refit.", call. = FALSE)
  fi <- file.info(LONG)
  if (!is.na(fi$size) && !identical(fi$size, meta$input$size))
    warning(LONG, " has changed since the cache was fitted (", meta$input$size,
            " -> ", fi$size, " bytes).\n  re-run  Rscript analyze_all_v2_fit.R  to refit.",
            call. = FALSE)
  invisible(NULL)
}
