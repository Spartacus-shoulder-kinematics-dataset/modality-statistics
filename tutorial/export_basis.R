# =============================================================================
# Export the REAL natural-spline basis and fitted coefficients, so the Python
# tutorial shows this project's actual model rather than a look-alike.
# Writes one CSV per df:  tutorial/basis_K<k>.csv
# and the fitted gammas:  tutorial/coefs.csv
#
# Run once:  Rscript tutorial/export_basis.R
# =============================================================================
suppressMessages({ library(splines) })
OUT <- "tutorial"
BOUNDARY <- c(0, 160)
GRID <- seq(BOUNDARY[1], BOUNDARY[2], length.out = 600)

# the knot rule used by analyze_all_v2.R: quantiles of x, rounded to nearest 10
knots_for <- function(K) switch(as.character(K),
  "2" = c(90), "3" = c(60, 110), "4" = c(50, 90, 120), "5" = c(40, 70, 100, 130))

for (K in 2:5) {
  kn <- knots_for(K)
  B  <- ns(GRID, knots = kn, Boundary.knots = BOUNDARY)
  stopifnot(ncol(B) == K)
  df <- data.frame(x = GRID, B)
  names(df) <- c("x", paste0("B", 1:K))
  write.csv(df, file.path(OUT, sprintf("basis_K%d.csv", K)), row.names = FALSE)
  write.csv(data.frame(K = K, knot = kn), file.path(OUT, sprintf("knots_K%d.csv", K)),
            row.names = FALSE)
  cat(sprintf("K=%d: %d basis columns, %d interior knots (%s)\n",
              K, ncol(B), length(kn), paste(kn, collapse = ", ")))
}

# fitted gammas from the worked example, so the sliders can start at the real fit
w <- file.path("figures", "04_natural_spline_hlme", "00_wald_tests_by_df.csv")
if (file.exists(w)) {
  d <- read.csv(w, stringsAsFactors = FALSE)
  d <- d[grepl("cond", d$term), c("K", "term", "estimate")]
  d$gamma <- ifelse(grepl("^cond", d$term), 0L,
                    suppressWarnings(as.integer(sub("^ns([0-9]+).*$", "\\1", d$term))))
  write.csv(d[order(d$K, d$gamma), c("K","gamma","term","estimate")],
            file.path(OUT, "coefs.csv"), row.names = FALSE)
  cat("wrote coefs.csv (fitted difference coefficients, per df)\n")
} else cat("NOTE: no fitted coefficients found; sliders will start at zero\n")
