# ============================================================
# Polynomial Fits on Group Mean Curves (ST DoF=2)
# Robust model selection: AICc + weighted LOOCV RMSE
# ============================================================

library(dplyr)
library(ggplot2)

# --- Lecture & filtrage ---
df <- read.csv("corrected_confident_data.csv", stringsAsFactors = FALSE)

st <- df %>%
  filter(
    joint == "scapulothoracic",
    humeral_motion == "frontal plane elevation",
    degree_of_freedom == 2,
    unit == "rad"
  ) %>%
  mutate(
    in_vivo_chr = tolower(as.character(in_vivo)),
    in_vivo_num = ifelse(in_vivo_chr %in% c("true","1","yes","y","t"), 1, 0),
    vivo_label = ifelse(in_vivo_num == 1, "in vivo", "ex vivo"),
    angle = value,
    x = humerothoracic_angle
  ) %>%
  filter(!is.na(x), !is.na(angle))

# --- Agrégation en courbes moyennes par groupe et par x ---
st_grouped <- st %>%
  group_by(vivo_label, x) %>%
  summarise(
    mean_angle = mean(angle, na.rm = TRUE),
    sd_angle   = sd(angle, na.rm = TRUE),
    n          = n(),
    .groups = "drop"
  ) %>%
  mutate(
    se = sd_angle / sqrt(pmax(n, 1)),
    # Choix des poids — par défaut: poids = n (stable).
    # Alternative (plus théorique) = 1/(se^2 + eps).
    w  = pmax(n, 1),
    # w = 1 / (se^2 + 1e-8)
  )

cat("\nPoints moyens par groupe:\n"); print(table(st_grouped$vivo_label))

# ---------- Helpers: AICc et RMSE-LOOCV pondéré ----------
compute_aicc <- function(mod, n_eff) {
  aic <- AIC(mod)
  k <- length(coef(mod)) + 1  # paramètres (coeffs + sigma)
  if (n_eff - k - 1 > 0) aic + (2 * k * (k + 1)) / (n_eff - k - 1) else NA_real_
}

cv_rmse_weighted <- function(dat, order) {
  # LOOCV par point moyen (laisse 1 point-x moyen hors fit)
  n <- nrow(dat)
  if (n <= order + 1) return(NA_real_)  # pas assez de points pour CV
  se_sum <- 0; w_sum <- 0; used <- 0
  for (i in seq_len(n)) {
    train <- dat[-i, , drop = FALSE]
    if (nrow(train) <= order) next
    mod_i <- lm(mean_angle ~ poly(x, order, raw = TRUE), data = train, weights = train$w)
    yhat  <- predict(mod_i, newdata = dat[i, , drop = FALSE])
    err2  <- (dat$mean_angle[i] - yhat)^2
    wi    <- dat$w[i]
    se_sum <- se_sum + wi * err2
    w_sum  <- w_sum + wi
    used   <- used + 1
  }
  if (used == 0 || w_sum == 0) return(NA_real_)
  sqrt(se_sum / w_sum)
}

# ---------- Fit & plot par groupe avec métriques robustes ----------
fit_group_poly <- function(data, group_name) {
  cat("\n==============================\n")
  cat(">>>", toupper(group_name), "\n")
  cat("==============================\n")
  
  results <- data.frame()
  color_fit <- ifelse(group_name == "in vivo", "#1b9e77", "#d95f02")
  
  for (order in 1:4) {
    # Fit polynomial pondéré (WLS)
    mod <- lm(mean_angle ~ poly(x, order, raw = TRUE), data = data, weights = data$w)
    
    # Metrics
    aicc  <- compute_aicc(mod, n_eff = nrow(data))
    cvrmse <- cv_rmse_weighted(data, order = order)
    
    cat(sprintf("\nOrder %d | AICc: %.2f | CV-RMSE_w: %.4f", order, aicc, cvrmse))
    
    results <- rbind(results, data.frame(
      Group = group_name,
      Order = order,
      AICc = aicc,
      CV_RMSE_w = cvrmse
    ))
    
    # Courbe prédite (uniquement sur le domaine observé du groupe)
    x_pred <- seq(min(data$x), max(data$x), length.out = 200)
    pred_df <- data.frame(x = x_pred)
    pred_df$y_pred <- predict(mod, newdata = pred_df)
    
    # Graphique: points gris + courbe modèle colorée
    p <- ggplot(data, aes(x = x, y = mean_angle)) +
      geom_point(color = "grey60", alpha = 0.6, size = 2) +
      geom_line(data = pred_df, aes(y = y_pred), color = color_fit, linewidth = 1.4) +
      labs(
        title = paste(group_name, "- Polynomial fit (order", order, ")"),
        subtitle = sprintf("AICc: %.2f | CV-RMSE_w: %.4f", aicc, cvrmse),
        x = "Humerothoracic angle (radians)",
        y = "Mean scapulothoracic angle (radians)"
      ) +
      theme_minimal(base_size = 14)
    print(p)
  }
  
  # Choix du meilleur ordre: min CV-RMSE_w, puis min AICc en cas d'égalité
  results <- results %>%
    arrange(CV_RMSE_w, AICc) %>%
    mutate(Best = ifelse(row_number() == 1, "★", ""))
  
  cat("\n\n-- Résumé métriques (triées) --\n"); print(results)
  invisible(results)
}

res_in_vivo <- fit_group_poly(filter(st_grouped, vivo_label == "in vivo"), "in vivo")
res_ex_vivo <- fit_group_poly(filter(st_grouped, vivo_label == "ex vivo"), "ex vivo")

# Tableau final
res_all <- bind_rows(res_in_vivo, res_ex_vivo)
cat("\n===== Récapitulatif global =====\n"); print(res_all)
