# Growth curve analysis
# a dummy example of (not latent) growth curve analysis
# based on:
# https://stats.stackexchange.com/questions/354732/latent-growth-curve-model-with-more-time-points-than-participants
library(lme4) # for mixed (or multilevel) model
library(lattice) # the plot library that goes hand-in-hand with lme4
library(here)

# need to install some libraries to make it work :)


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

