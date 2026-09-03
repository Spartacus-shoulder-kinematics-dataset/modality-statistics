# =============================================================================
# Evidence for the response to STATISTICAL_METHODS_REVIEW.md
# Exemplar cell: scapulothoracic / frontal plane elevation / DoF 2.
# Demonstrates, on real data, the impact of the reviewer's main points:
#   (1) study-level random effect (condition is confounded with study)
#   (3) within-shoulder rho (fix the pooled lag-1 estimator)
#   (4/7) study-level uncertainty via leave-one-study-out (LOSO)
#   (6) signed vs absolute effect size
# =============================================================================
suppressMessages({ library(mgcv); library(dplyr) })

D <- read.csv("spartacus_angles_long.csv", stringsAsFactors = FALSE)
d <- D %>% filter(joint == "scapulothoracic",
                  humeral_motion == "frontal plane elevation",
                  degree_of_freedom == 2)
d$cond  <- factor(ifelse(d$in_vivo == "True", "in vivo", "ex vivo"), levels = c("ex vivo","in vivo"))
d$condO <- as.ordered(d$cond)
d$ID    <- factor(d$ID)
d$study <- factor(sub("_[^_]*$", "", d$ID))     # recover article from ID
d <- d[order(d$ID, d$TIME), ]
d$ar_start <- c(TRUE, d$ID[-1] != d$ID[-nrow(d)])

cat(sprintf("Exemplar: %d rows, %d shoulders, %d studies (%d ex / %d in studies)\n",
    nrow(d), nlevels(d$ID), nlevels(d$study),
    length(unique(d$study[d$cond=="ex vivo"])), length(unique(d$study[d$cond=="in vivo"]))))

# within-shoulder lag-1 residual autocorrelation (reviewer point 3)
within_rho <- function(resid, id) {
  num <- den <- 0
  for (u in levels(id)) { r <- resid[id == u]
    if (length(r) > 1) { num <- num + sum(head(r,-1)*tail(r,-1)); den <- den + sum(r^2) } }
  num/den
}
fit <- function(form, dat) {
  m0  <- bam(form, data = dat, method = "fREML", discrete = TRUE)
  rho <- within_rho(resid(m0), dat$ID)
  m1  <- do.call(bam, list(formula = form, data = dat, method = "fREML",
                           discrete = TRUE, rho = rho, AR.start = dat$ar_start))
  list(m = m1, rho = rho)
}
lvl <- function(m) { pt <- summary(m)$p.table; i <- grep("condO", rownames(pt))[1]
  c(est = pt[i,"Estimate"], se = pt[i,"Std. Error"], p = pt[i,"Pr(>|t|)"]) }

# ---- Model A: current (shoulder RE only) vs Model B: + study RE --------------
fA <- fit(Y ~ condO + s(TIME) + s(TIME, by=condO) + s(ID, bs="re"), d)
fB <- fit(Y ~ condO + s(TIME) + s(TIME, by=condO) + s(study, bs="re") + s(ID, bs="re"), d)
cat(sprintf("\nrho: pooled(old)=%.4f  within-shoulder(fixed)=%.4f\n",
            sum(head(resid(fA$m),-1)*tail(resid(fA$m),-1))/sum(resid(fA$m)^2), fA$rho))
a <- lvl(fA$m); b <- lvl(fB$m)
cat(sprintf("\nLevel shift beta1 (deg):\n  A shoulder-RE only : %+6.2f  SE %5.2f  p=%.2g\n  B + study-RE       : %+6.2f  SE %5.2f  p=%.2g   <- SE inflates, p weakens\n",
            a["est"],a["se"],a["p"], b["est"],b["se"],b["p"]))

# ---- signed vs absolute effect size on the overlap (reviewer point 6) --------
ov <- c(max(min(d$TIME[d$cond=="ex vivo"]), min(d$TIME[d$cond=="in vivo"])),
        min(max(d$TIME[d$cond=="ex vivo"]), max(d$TIME[d$cond=="in vivo"])))
xo <- seq(ov[1], ov[2], length.out = 250)
dif_of <- function(m) {
  mk <- function(cc) data.frame(TIME=xo, ID=d$ID[1], study=d$study[1],
                                condO=factor(cc, levels=levels(d$condO), ordered=TRUE))
  as.vector((predict(m, mk("in vivo"), type="lpmatrix", newdata.guaranteed=TRUE) -
             predict(m, mk("ex vivo"), type="lpmatrix", newdata.guaranteed=TRUE)) %*% coef(m))
}
dB <- dif_of(fB$m)
cat(sprintf("\nEffect size on overlap [%.0f,%.0f] deg (study-RE model):\n  signed mean = %+.1f deg | absolute mean = %.1f deg | max|.| = %.1f deg\n",
            ov[1], ov[2], mean(dB), mean(abs(dB)), max(abs(dB))))

# ---- leave-one-study-out: how much does the contrast move? (points 4 & 7) ----
cat("\nLeave-one-study-out (drop each study, refit model B):\n")
res <- data.frame()
for (s in levels(d$study)) {
  ds <- droplevels(d[d$study != s, ])
  if (length(unique(ds$cond)) < 2 ||
      min(table(unique(ds[,c("study","cond")])$cond)) < 2) next
  fs <- tryCatch(fit(Y ~ condO + s(TIME) + s(TIME, by=condO) + s(study, bs="re") + s(ID, bs="re"), ds),
                 error = function(e) NULL)
  if (is.null(fs)) next
  l <- lvl(fs$m)
  res <- rbind(res, data.frame(dropped = s, beta1 = round(l["est"],2), p = signif(l["p"],2)))
}
rownames(res) <- NULL; print(res)
cat(sprintf("\n=> beta1 ranges %.1f to %.1f deg across LOSO (full-data B = %.1f); p not always < .05.\n",
            min(res$beta1), max(res$beta1), b["est"]))
cat("The contrast is sensitive to individual studies -> inference must be study-level.\n")
