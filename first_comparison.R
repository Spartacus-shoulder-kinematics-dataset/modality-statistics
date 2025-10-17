# todo try: https://m-clark.github.io/clustered-data/latent-growth-curves.html
# Load necessary libraries
library(dplyr)
library(tidyr)
library(lavaan)
library(lme4)
library(here)


# --- Load and filter ---
my_df <- read.csv("corrected_confident_data.csv")

gh_frontal_plane_dof_1 <- my_df
  filter(joint == "glenohumeral",
         humeral_motion == "frontal plane elevation",
         degree_of_freedom == 1,
         unit == "rad",
         )

gh_frontal_plane_dof_1 <- gh_frontal_plane_dof_1 %>%
  mutate(combined_subject = paste(article, shoulder_id, sep = "_"))

# make y_trial_number
gh_frontal_plane_dof_1 <- gh_frontal_plane_dof_1 %>%
  mutate(y_trial_number = paste0("y", trial_number))


# --- Convert to wide format ---
gh_wide <- gh_frontal_plane_dof_1 %>%
  pivot_wider(id_cols = c(combined_subject, in_vivo),
              names_from = y_trial_number,
              values_from = value)

# Convert in_vivo to numeric (0 for ex_vivo, 1 for in_vivo)
gh_wide <- gh_wide %>%
  mutate(in_vivo = ifelse(in_vivo == "True", 1, 0))

# --- LGC Model Specification with Group Interaction ---

model_syntax <- '
    # Latent factors (intercept and slope)
    i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4

    s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4

    # Regress latent factors on in_vivo (group indicator)
    i ~ in_vivo
    s ~ in_vivo
    
    # Equal residual variances
    y1 ~~ resvar*y1
    y2 ~~ resvar*y2
    y3 ~~ resvar*y3
    y4 ~~ resvar*y4
    
'

# --- Model Estimation with lavaan ---
growth_model <- growth(model_syntax, data = gh_wide)
summary(growth_model, fit.measures=TRUE, standardized=TRUE)



# --- LGC Model Specification without Group Interaction ---

model_syntax_no_interaction <- '
    # Latent factors (intercept and slope)
    i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4

    s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4

    # Equal residual variances
    y1 ~~ resvar*y1
    y2 ~~ resvar*y2
    y3 ~~ resvar*y3
    y4 ~~ resvar*y4
    
'

# --- Model Estimation with lavaan ---
growth_model_no_interaction <- growth(model_syntax_no_interaction, data = gh_wide)
summary(growth_model_no_interaction, fit.measures=TRUE, standardized=TRUE)


# --- Compare Models ---
anova(growth_model_no_interaction, growth_model)