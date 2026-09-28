# =============================================================================
# utils_lmethod.R — the L method for finding the knee of a fit-vs-complexity
# curve, as a cross-check on BIC when choosing the spline df K.
#
# Salvador, S. & Chan, P. (2004). Determining the Number of Clusters/Segments
# in Hierarchical Clustering/Segmentation Algorithms. ICTAI.
#
# WHY THIS EXISTS
#   BIC picks K by trading likelihood against a k*log(n) penalty. With ~10^4
#   observations that penalty is tiny next to what an extra spline df buys, so
#   BIC can keep improving all the way to the largest K on the grid: it then
#   reports the edge of the grid, not a choice. The L method asks a different
#   and penalty-free question — where does the error-vs-K curve stop being
#   steep and start being flat? — by fitting TWO straight lines to the curve
#   and taking their meeting point as the knee.
#
#   The two answers are complementary, not redundant:
#     BIC       "is the extra parameter paid for by the likelihood?"  (absolute)
#     L method  "has the curve already flattened out?"                (shape)
#   Report both. When they agree, K is not controversial. When BIC runs to the
#   edge of the grid and the knee sits well inside it, the knee is the
#   defensible K and the disagreement is itself the finding.
#
# WHAT THE CURVE MUST LOOK LIKE
#   The L method wants an evaluation graph that is steep on one side, flat on
#   the other, i.e. an "L". Feed it an UNPENALIZED goodness-of-fit measure
#   against K — deviance (-2*loglik), residual RMSE, per-shoulder RMSE. Feeding
#   it BIC itself is allowed and sometimes informative, but then the knee is a
#   knee of an already-penalized curve and no longer an independent check.
#
# SCALE INVARIANCE
#   The knee is unchanged by any affine rescaling of x or of y (least-squares
#   residuals are invariant to affine maps of the predictor, and scale with y).
#   So there is no need to standardize the axes, and no hidden dependence on
#   the units of Y. A log transform of y is NOT affine — it genuinely changes
#   what counts as "linear", so `transform = "log"` is a deliberate modelling
#   choice, not a normalization. Use it when the raw curve is so convex that
#   the knee is pinned at the second point.
#
# KNOWN LIMITS (stated in the paper, and they bite here)
#   - It cannot return the first or last point of the grid: each line needs two
#     points, so the knee lives in x[2] .. x[n-2]. With K = 1..10 the reachable
#     answers are K = 2..8.
#   - The answer depends on how far the grid extends. A long flat tail drags
#     the knee right. `lmethod_sensitivity()` shows exactly how much, and you
#     should look at it before quoting a knee.
#   - It cannot tell you that a straight line (K = 1) was enough.
#
# Base R only, no dependencies. Source it:
#     source("utils_lmethod.R")
# or run the self-test and the demo on this project's own BIC table:
#     Rscript utils_lmethod.R
# =============================================================================


# -----------------------------------------------------------------------------
# 1. CORE — one pass of the L method
# -----------------------------------------------------------------------------

# RMSE of the least-squares line through (x, y). Two points define a line
# exactly, hence RMSE 0, which is what makes the extreme partitions cheap and
# why the weighting in `lmethod_scan()` is needed to keep them from winning.
.rmse_line <- function(x, y) {
  n <- length(x)
  if (n < 2L) return(NA_real_)
  if (n == 2L) return(0)
  if (length(unique(x)) < 2L) return(sqrt(mean((y - mean(y))^2)))
  r <- .lm.fit(cbind(1, x), y)$residuals
  sqrt(mean(r^2))
}

#' Every admissible two-line partition of an evaluation curve
#'
#' Left  = points 1..j, right = points (j+1)..n, so the split point j belongs to
#' the left line and the two lines share no point. j runs 2..(n-2) so that each
#' side keeps at least two points. The total error is the length-weighted mean
#' of the two RMSEs (paper, eq. 1) — without the weights the partitions that
#' put two points on one side would always win, since a line through two points
#' has zero error.
#'
#' @param x complexity axis (e.g. spline df K), strictly increasing
#' @param y evaluation metric at each x (lower = better fit, or higher, either way)
#' @return data.frame, one row per candidate split, ordered by x_knee
lmethod_scan <- function(x, y) {
  n <- length(x)
  if (n < 4L) stop("the L method needs at least 4 points; got ", n)
  js <- 2L:(n - 2L)
  do.call(rbind, lapply(js, function(j) {
    li <- seq_len(j); ri <- (j + 1L):n
    rl <- .rmse_line(x[li], y[li]); rr <- .rmse_line(x[ri], y[ri])
    data.frame(split_index = j, x_knee = x[j],
               n_left = length(li), n_right = length(ri),
               rmse_left = rl, rmse_right = rr,
               rmse = (length(li) / n) * rl + (length(ri) / n) * rr)
  }))
}

#' One pass of the L method on the whole curve given
#'
#' @return list with the knee, its index, the winning partition's weighted RMSE,
#'   the single-line RMSE it is compared against, and the full scan table.
lmethod_once <- function(x, y) {
  sc <- lmethod_scan(x, y)
  # Ties are real, not just floating-point: if the curve is exactly two straight
  # lines, the vertex point lies on BOTH of them, so the partitions on either
  # side of it are equally good and the winner would otherwise be decided by
  # rounding noise. Resolve ties toward the SMALLER split — when two values of K
  # describe the curve equally well, the parsimonious one is the answer.
  tol <- 1e-9 * max(sc$rmse, na.rm = TRUE)
  b   <- which(sc$rmse <= min(sc$rmse, na.rm = TRUE) + tol)[1]
  j   <- sc$split_index[b]
  # How much better is "two lines" than "one line"? If this is near 0 the curve
  # has no elbow and the reported knee is an artefact of where the grid stops.
  # When a single line already fits essentially exactly, the ratio is 0/0 and no
  # elbow can be claimed either way — report NA rather than a noise-driven number.
  rmse_one <- .rmse_line(x, y)
  degenerate <- !is.finite(rmse_one) || rmse_one <= 1e-10 * max(abs(y), na.rm = TRUE)
  list(knee = sc$x_knee[b], index = j, rmse = sc$rmse[b],
       rmse_single_line = rmse_one,
       strength = if (degenerate) NA_real_ else 1 - sc$rmse[b] / rmse_one,
       coef_left  = stats::coef(stats::lm(y[seq_len(j)] ~ x[seq_len(j)])),
       coef_right = stats::coef(stats::lm(y[(j + 1L):length(x)] ~ x[(j + 1L):length(x)])),
       scan = sc)
}


# -----------------------------------------------------------------------------
# 2. THE FULL METHOD — trimming, iterative refinement, transforms
# -----------------------------------------------------------------------------

# Paper 3.3, second refinement: a curve that rises to a maximum before falling
# (reading from the flat end towards the steep end) has no "L" in it. Drop the
# leading points up to that turning point. Orientation is read off the overall
# trend, so this is a no-op for the monotone curves we normally feed it.
.trim_leading_extremum <- function(x, y) {
  n <- length(y)
  if (n < 5L) return(seq_len(n))
  i <- if (y[n] < y[1]) which.max(y) else which.min(y)   # decreasing vs increasing
  if (i >= n - 3L) return(seq_len(n))                    # refuse to gut the curve
  i:n
}

#' The L method, with the paper's iterative refinement
#'
#' Refinement (paper 3.3, first refinement) re-runs the method on the first
#' `2 * knee` points and repeats while the knee keeps moving left. It matters
#' when the curve has hundreds of points and the flat tail drowns out the steep
#' head. The paper never lets the focus region drop below ~20 points, so for a
#' short grid such as K = 1..10 refinement is a no-op by construction — that is
#' correct behaviour, not a silent failure, and `$n_refine_iters` reports it.
#'
#' @param x complexity axis, strictly increasing (e.g. K = 1:10)
#' @param y evaluation metric — an UNPENALIZED fit measure, see the header
#' @param refine run the iterative refinement (default TRUE)
#' @param min_cutoff smallest focus region the refinement may use (paper: 20)
#' @param trim drop leading points before a turning point (default FALSE)
#' @param transform "identity" (default) or "log" — see the header on scale
#' @return object of class "lmethod"; print() and plot_lmethod() know about it
lmethod <- function(x, y, refine = TRUE, min_cutoff = 20L,
                    trim = FALSE, transform = c("identity", "log")) {
  transform <- match.arg(transform)
  stopifnot(length(x) == length(y))

  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 4L)
    stop("need at least 4 finite (x, y) pairs; got ", sum(keep))
  x0 <- x[keep]; y0 <- y[keep]
  o  <- order(x0); x0 <- x0[o]; y0 <- y0[o]
  dropped <- sum(!keep)

  yt <- y0
  if (transform == "log") {
    if (any(yt <= 0)) stop("transform = \"log\" needs strictly positive y ",
                           "(shift a deviance curve first, or use identity)")
    yt <- log(yt)
  }

  idx <- if (trim) .trim_leading_extremum(x0, yt) else seq_along(x0)
  xs <- x0[idx]; ys <- yt[idx]
  if (length(xs) < 4L) stop("fewer than 4 points left after trimming")

  # --- iterative refinement -------------------------------------------------
  n <- length(xs)
  cutoff <- n; last_knee <- n + 1L; iters <- 0L
  fit <- NULL
  repeat {
    m   <- max(min_cutoff, 4L)
    lim <- min(n, max(m, cutoff))
    f   <- lmethod_once(xs[seq_len(lim)], ys[seq_len(lim)])
    iters <- iters + 1L
    if (!refine) { fit <- f; break }
    if (f$index >= last_knee || lim <= m) { fit <- f; break }
    last_knee <- f$index
    fit <- f
    cutoff <- f$index * 2L
    if (iters > 50L) break                      # cannot happen; cheap insurance
  }

  structure(list(
    knee = fit$knee, index = fit$index,
    rmse = fit$rmse, rmse_single_line = fit$rmse_single_line,
    strength = fit$strength,
    coef_left = fit$coef_left, coef_right = fit$coef_right,
    scan = fit$scan,
    x = xs, y = ys, y_raw = y0[idx], transform = transform,
    reachable = range(xs[2:(length(xs) - 2L)]),
    n_refine_iters = iters, refined = refine && iters > 1L,
    n_trimmed = length(x0) - length(idx), n_dropped = dropped
  ), class = "lmethod")
}

#' @export
print.lmethod <- function(x, ...) {
  cat("L method (Salvador & Chan 2004)\n")
  cat(sprintf("  curve            : %d points, x in %g..%g%s\n",
              length(x$x), min(x$x), max(x$x),
              if (x$transform != "identity") sprintf(", y on the %s scale", x$transform) else ""))
  cat(sprintf("  KNEE             : x = %g\n", x$knee))
  cat(sprintf("  two-line RMSE    : %.4g   (one line: %.4g)\n", x$rmse, x$rmse_single_line))
  cat(sprintf("  knee strength    : %.1f%% error reduction over a single line%s\n",
              100 * x$strength,
              if (is.finite(x$strength) && x$strength < 0.5) "   <- WEAK: barely an elbow" else ""))
  cat(sprintf("  reachable knees  : x = %g..%g  (each line needs 2 points)\n",
              x$reachable[1], x$reachable[2]))
  if (x$knee <= x$reachable[1])
    cat("  NOTE: the knee is at the smallest value the method can return — the true\n",
        "        knee may be further left, i.e. off this grid.\n", sep = "")
  if (x$knee >= x$reachable[2])
    cat("  NOTE: the knee is at the largest value the method can return — extend the\n",
        "        grid before quoting it.\n", sep = "")
  if (x$n_trimmed) cat(sprintf("  trimmed          : %d leading point(s)\n", x$n_trimmed))
  if (x$n_dropped) cat(sprintf("  dropped          : %d non-finite point(s)\n", x$n_dropped))
  cat(sprintf("  refinement       : %d iteration(s)%s\n", x$n_refine_iters,
              if (x$n_refine_iters == 1L) " (grid too short to refine — expected)" else ""))
  invisible(x)
}

#' How much does the knee depend on where the grid stops?
#'
#' Re-runs the L method on the first m points for every m, which is the paper's
#' main weakness made visible: if the knee only appears once the grid is long
#' enough, or drifts right with every extra point, say so instead of quoting one
#' number. A knee that is stable over a run of truncations is one to trust.
#'
#' @return data.frame(m, x_max, knee, strength)
lmethod_sensitivity <- function(x, y, refine = FALSE, ...) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]; o <- order(x); x <- x[o]; y <- y[o]
  do.call(rbind, lapply(4:length(x), function(m) {
    f <- try(lmethod(x[seq_len(m)], y[seq_len(m)], refine = refine, ...), silent = TRUE)
    if (inherits(f, "try-error")) return(NULL)
    data.frame(m = m, x_max = x[m], knee = f$knee, strength = f$strength)
  }))
}


# -----------------------------------------------------------------------------
# 3. THE COMPARISON THIS PROJECT ACTUALLY WANTS — knee vs BIC
# -----------------------------------------------------------------------------

#' Choose the spline df K by the L method and by BIC, and compare
#'
#' @param K       vector of spline df tried
#' @param metric  UNPENALIZED fit metric at each K, lower = better fit
#'                (deviance = -2*loglik, or rmse_median, ...). If you pass
#'                loglik, set `metric_is_loglik = TRUE` and it is converted.
#' @param BIC     BIC at each K
#' @param converged optional logical; non-converged fits are excluded from both
#' @param metric_name label for plots and printing
#' @param ... passed to lmethod() (refine, trim, transform)
#' @return object of class "lmethod_bic"
lmethod_vs_bic <- function(K, metric, BIC, converged = NULL,
                           metric_name = "deviance (-2 logLik)",
                           metric_is_loglik = FALSE, ...) {
  if (metric_is_loglik) { metric <- -2 * metric; metric_name <- "deviance (-2 logLik)" }
  ok <- is.finite(K) & is.finite(metric) & is.finite(BIC)
  if (!is.null(converged)) ok <- ok & as.logical(converged)
  if (sum(ok) < 4L) stop("need at least 4 converged K values; got ", sum(ok))
  K <- K[ok]; metric <- metric[ok]; BIC <- BIC[ok]
  o <- order(K); K <- K[o]; metric <- metric[o]; BIC <- BIC[o]

  fit <- lmethod(K, metric, ...)
  K_knee <- fit$knee
  K_bic  <- K[which.min(BIC)]
  # BIC that stops at the edge of the grid has not selected anything — it has
  # run out of candidates. This is the case worth flagging out loud.
  bic_at_edge <- K_bic == max(K)

  tab <- data.frame(
    K = K, metric = metric, BIC = BIC,
    dBIC = BIC - min(BIC),
    d_metric = c(NA, diff(metric)),
    pick = ifelse(K == K_knee & K == K_bic, "knee+BIC",
           ifelse(K == K_knee, "knee",
           ifelse(K == K_bic, "BIC", ""))))

  structure(list(fit = fit, table = tab, K = K, metric = metric, BIC = BIC,
                 K_knee = K_knee, K_bic = K_bic,
                 agree = K_knee == K_bic,
                 bic_at_edge = bic_at_edge,
                 dBIC_at_knee = BIC[K == K_knee] - min(BIC),
                 metric_name = metric_name),
            class = "lmethod_bic")
}

#' @export
print.lmethod_bic <- function(x, ...) {
  cat("\n---- choosing K: the L method vs BIC ----\n\n")
  print(x$fit)
  cat(sprintf("\n  BIC minimum      : K = %d  (BIC %.1f)\n", x$K_bic, min(x$BIC)))
  cat(sprintf("  L-method knee    : K = %d  on %s\n", x$K_knee, x$metric_name))
  cat(sprintf("  BIC cost of the knee: %+.1f BIC vs the BIC minimum\n", x$dBIC_at_knee))
  cat("\n")
  print(x$table, digits = 6, row.names = FALSE)
  cat("\n")
  if (x$agree) {
    cat(sprintf("  VERDICT: both criteria say K = %d. Nothing to arbitrate.\n", x$K_knee))
  } else if (x$bic_at_edge) {
    cat(sprintf("  VERDICT: BIC's minimum is at K = %d, the LARGEST K tried — BIC has not\n", x$K_bic))
    cat( "           turned around, so it is reporting the edge of the grid rather than\n")
    cat( "           selecting. Either extend the grid until BIC turns, or take the knee\n")
    cat(sprintf("           at K = %d, where the fit curve flattens, as the reported df.\n", x$K_knee))
  } else {
    cat(sprintf("  VERDICT: BIC says K = %d, the knee says K = %d. BIC has an interior\n",
                x$K_bic, x$K_knee))
    cat(sprintf("           minimum, so it IS selecting; the knee's %d costs %.1f BIC.\n",
                x$K_knee, x$dBIC_at_knee))
    cat( "           Prefer BIC unless the extra df buy nothing you can see in the fit.\n")
  }
  if (is.finite(x$fit$strength) && x$fit$strength < 0.5)
    cat("  CAUTION: weak elbow — the two-line fit barely beats one line. Treat the\n",
        "           knee as indicative only.\n", sep = "")
  invisible(x)
}


# -----------------------------------------------------------------------------
# 4. PLOTS — base graphics, same palette as the analysis scripts
# -----------------------------------------------------------------------------
.LM_BLUE <- "#377eb8"; .LM_ORANGE <- "#d95f02"; .LM_GREEN <- "#1b9e77"

#' The evaluation curve with the two fitted lines and the knee
#'
#' @param fit an "lmethod" object
#' @param vlines optional named numeric, extra reference lines (e.g. c("BIC" = 10))
plot_lmethod <- function(fit, xlab = "number of parameters (K)",
                         ylab = "evaluation metric", main = NULL,
                         vlines = NULL, legend_pos = "topright") {
  x <- fit$x; y <- fit$y
  if (is.null(main)) main <- sprintf("L method — knee at K = %g", fit$knee)
  plot(x, y, type = "n", xlab = xlab, ylab = ylab, main = main, xaxt = "n")
  axis(1, at = x); grid(col = "grey90", lty = 1)

  j <- fit$index
  # left line over its own points; right line extended back to the knee so the
  # two segments meet on the plot the way the method conceives them
  cl <- fit$coef_left; cr <- fit$coef_right
  segments(x[1], cl[1] + cl[2] * x[1], x[j], cl[1] + cl[2] * x[j],
           col = .LM_ORANGE, lwd = 2.6)
  segments(x[j], cr[1] + cr[2] * x[j], x[length(x)], cr[1] + cr[2] * x[length(x)],
           col = .LM_GREEN, lwd = 2.6)

  lines(x, y, lwd = 2, col = "grey55")
  points(x, y, pch = 21, bg = "white", col = "grey35", cex = 1.2, lwd = 1.4)
  abline(v = fit$knee, lty = 3, col = .LM_BLUE, lwd = 2)
  points(fit$knee, y[j], pch = 21, bg = .LM_BLUE, col = "white", cex = 1.9, lwd = 2)

  lg <- c(sprintf("knee: K = %g", fit$knee), "steep side", "flat side")
  lc <- c(.LM_BLUE, .LM_ORANGE, .LM_GREEN); ll <- c(3, 1, 1)
  if (!is.null(vlines)) {
    for (i in seq_along(vlines)) abline(v = vlines[i], lty = 2, col = "grey30", lwd = 1.6)
    lg <- c(lg, sprintf("%s: K = %g", names(vlines), vlines))
    lc <- c(lc, rep("grey30", length(vlines))); ll <- c(ll, rep(2, length(vlines)))
  }
  legend(legend_pos, bty = "n", cex = 0.8, lwd = 2.2, lty = ll, col = lc, legend = lg)
  legend("bottomleft", bty = "n", cex = 0.72, text.col = "grey30",
         legend = sprintf("two-line RMSE %.3g vs %.3g for one line (%.0f%% better)",
                          fit$rmse, fit$rmse_single_line, 100 * fit$strength))
  invisible(fit)
}

#' Two panels: the knee on the fit curve, and BIC, marked with both answers
plot_lmethod_vs_bic <- function(cmp, xlab = "natural-spline df (K)") {
  op <- par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.4, 1)); on.exit(par(op), add = TRUE)

  plot_lmethod(cmp$fit, xlab = xlab, ylab = cmp$metric_name,
               main = sprintf("L method on the fit curve — knee K = %d", cmp$K_knee),
               vlines = stats::setNames(cmp$K_bic, "BIC min"))

  plot(cmp$K, cmp$BIC, type = "n", xlab = xlab, ylab = "BIC (lower is better)",
       xaxt = "n", main = sprintf("BIC — minimum at K = %d%s", cmp$K_bic,
                                  if (cmp$bic_at_edge) " (grid edge)" else ""))
  axis(1, at = cmp$K); grid(col = "grey90", lty = 1)
  lines(cmp$K, cmp$BIC, lwd = 2, col = "grey55")
  points(cmp$K, cmp$BIC, pch = 21, bg = "white", col = "grey35", cex = 1.2, lwd = 1.4)
  points(cmp$K_bic, cmp$BIC[cmp$K == cmp$K_bic], pch = 21, bg = "grey30",
         col = "white", cex = 1.8, lwd = 2)
  abline(v = cmp$K_knee, lty = 3, col = .LM_BLUE, lwd = 2)
  points(cmp$K_knee, cmp$BIC[cmp$K == cmp$K_knee], pch = 21, bg = .LM_BLUE,
         col = "white", cex = 1.8, lwd = 2)
  legend("topright", bty = "n", cex = 0.8, lwd = 2.2, lty = c(2, 3),
         col = c("grey30", .LM_BLUE),
         legend = c(sprintf("BIC min: K = %d", cmp$K_bic),
                    sprintf("knee: K = %d (+%.0f BIC)", cmp$K_knee, cmp$dBIC_at_knee)))
  if (cmp$bic_at_edge)
    legend("bottomleft", bty = "n", cex = 0.72, text.col = .LM_ORANGE,
           legend = "BIC still falling at the last K — it is not selecting")
  invisible(cmp)
}

#' How the knee moves as the K grid is extended (the paper's main caveat)
plot_lmethod_sensitivity <- function(x, y, xlab = "largest K on the grid", ...) {
  s <- lmethod_sensitivity(x, y, ...)
  # keep the y = x diagonal on screen: a knee that simply tracks the grid edge
  # is not a knee, and that is only visible next to the diagonal
  plot(s$x_max, s$knee, type = "n", xlab = xlab, ylab = "knee returned",
       ylim = range(s$knee, s$x_max), main = "Does the knee survive extending the grid?",
       xaxt = "n")
  axis(1, at = s$x_max); grid(col = "grey90", lty = 1)
  abline(0, 1, col = "grey75", lty = 2, lwd = 1.6)
  lines(s$x_max, s$knee, lwd = 2.4, col = .LM_BLUE)
  points(s$x_max, s$knee, pch = 21, bg = .LM_BLUE, col = "white", cex = 1.4, lwd = 1.6)
  legend("topleft", bty = "n", cex = 0.75, text.col = "grey30",
         legend = c("flat run = a stable knee",
                    "tracking the diagonal = no elbow, just the grid edge"))
  invisible(s)
}


# -----------------------------------------------------------------------------
# 5. SELF-TEST + DEMO — only when run as a script
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L && !interactive()) {

  cat("\n================ utils_lmethod.R — self-test ================\n")

  # A curve with a knee planted between K = 5 and K = 6: steep (slope -100) up
  # to K = 5, flat (slope -2) from K = 6 on. The vertex is deliberately placed
  # BETWEEN two grid points, so exactly one partition is right — putting it on a
  # grid point makes that point lie on both lines and the answer ambiguous by 1.
  set.seed(1)
  K  <- 1:12
  yy <- ifelse(K <= 5, 1000 - 100 * K, 420 - 2 * (K - 6))
  f  <- lmethod(K, yy)
  cat(sprintf("  planted knee at K = 5, recovered K = %g   %s\n", f$knee,
              if (f$knee == 5) "OK" else "*** FAIL ***"))
  stopifnot(f$knee == 5)

  # affine invariance: rescaling either axis must not move the knee
  stopifnot(lmethod(K, yy * 1e6 + 7)$knee == 5)
  stopifnot(lmethod(K * 3 + 11, yy)$index == f$index)
  cat("  affine invariance in x and y                     OK\n")

  # a noisy straight line has no elbow: two lines must barely beat one.
  # (A NOISE-FREE straight line is 0/0 — one line already fits exactly — and is
  # reported as strength NA rather than as a noise-driven number.)
  fl <- lmethod(K, 500 - 10 * K + rnorm(length(K), sd = 2))
  cat(sprintf("  noisy straight line -> strength %.2f (want < 0.5)  %s\n",
              fl$strength, if (fl$strength < 0.5) "OK" else "*** FAIL ***"))
  stopifnot(fl$strength < 0.5)
  cat(sprintf("  exact straight line -> strength %s (0/0, no elbow)  %s\n",
              format(lmethod(K, 500 - 10 * K)$strength),
              if (is.na(lmethod(K, 500 - 10 * K)$strength)) "OK" else "*** FAIL ***"))

  # too few points must error, not guess
  stopifnot(inherits(try(lmethod(1:3, c(3, 2, 1)), silent = TRUE), "try-error"))
  cat("  refuses fewer than 4 points                      OK\n")

  # ---- run on iteration 06's real BIC table --------------------------------
  S6  <- "06_natural_spline_hlme"
  OUT <- file.path("figures", S6)
  csv <- file.path(OUT, "00_bic_by_k.csv")
  if (!file.exists(csv)) {
    cat("\n  (skipping: ", csv, " not found — run analysis_lrt.R first)\n\n", sep = "")
    quit(save = "no")
  }

  cat("\n================ iteration 06: the knee of the K curve ================\n", sep = "")
  b <- read.csv(csv)

  # same png/pdf pair the analysis scripts produce, same scale factor
  FIG_SCALE <- 2.2
  fig <- function(name, draw, w = 1100, h = 620, res = 130) {
    png(file.path(OUT, paste0(name, ".png")), width = round(w*FIG_SCALE),
        height = round(h*FIG_SCALE), res = round(res*FIG_SCALE))
    op <- par(no.readonly = TRUE); tryCatch(draw(), finally = { par(op); dev.off() })
    if (capabilities("cairo")) cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w/res, height = h/res)
    else pdf(file.path(OUT, paste0(name, ".pdf")), width = w/res, height = h/res)
    op <- par(no.readonly = TRUE); tryCatch(draw(), finally = { par(op); dev.off() })
    message("  saved ", file.path(S6, name), ".png/.pdf")
  }

  # (a) the headline comparison, on the deviance — the likelihood's own scale
  cmp <- lmethod_vs_bic(b$K, metric = b$loglik, BIC = b$BIC,
                        converged = b$converged, metric_is_loglik = TRUE)
  print(cmp)

  # (b) the same question of three other curves. The deviance falls by a factor
  #     of six across the grid, so its knee lands early; the log deviance and the
  #     per-shoulder RMSE are the fairer "has it flattened?" curves to look at.
  dev <- -2 * b$loglik
  curves <- list("deviance"          = lmethod(b$K, dev),
                 "log deviance"      = lmethod(b$K, dev, transform = "log"),
                 "median RMSE (deg)" = lmethod(b$K, b$rmse_median),
                 "max RMSE (deg)"    = lmethod(b$K, b$rmse_max))
  alt <- data.frame(curve = names(curves),
                    knee = sapply(curves, `[[`, "knee"),
                    strength = sapply(curves, `[[`, "strength"))
  cat("\n  the knee of four different fit curves (BIC's answer is K = ",
      cmp$K_bic, " for all of them):\n", sep = "")
  print(alt, digits = 3, row.names = FALSE)
  cat("\n  Agreement across curves is the thing to look for: they disagree by a\n",
      "  degree or two here, so quote a range, not a single K.\n", sep = "")

  # (c) does the knee survive extending the grid?
  sens <- lmethod_sensitivity(b$K[b$converged], dev[b$converged])
  cat("\n  sensitivity of the deviance knee to where the K grid stops:\n")
  print(sens, digits = 4, row.names = FALSE)

  fig("07_lmethod_vs_bic",   function() plot_lmethod_vs_bic(cmp), w = 1300, h = 620)
  fig("07_lmethod_sensitivity", function() {
    op <- par(mar = c(4.4, 4.6, 3.4, 1)); on.exit(par(op), add = TRUE)
    plot_lmethod_sensitivity(b$K[b$converged], dev[b$converged])
  }, w = 800, h = 620)

  write.csv(cmp$table, file.path(OUT, "07_lmethod_vs_bic.csv"), row.names = FALSE)
  message("  saved ", file.path(S6, "07_lmethod_vs_bic.csv"))
  cat("\n")
}
