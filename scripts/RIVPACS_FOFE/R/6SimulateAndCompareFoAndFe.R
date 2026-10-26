# --------------------------------------------------------------------
# Test the null hypothesis that Fo = Fe.
# Hawkins & Yuan 2016 Freshwater Science
# --------------------------------------------------------------------


# --------------------------------------------------------------------
# Taxa probability threshold
pth <- 0.0001 # t?m? tarttisi tutkimusta !!
# Number of simulations (5000 in Hawkins & Yuan 2016)
nsimu <- 5000
# --------------------------------------------------------------------


# --------------------------------------------------------------------
# # Name the site set
# SiteSetName <- "validation"
# # RIVPACS-modelled taxa capture probabilities data frame
# capture_probs <- OERFvld$Capture.Probs
# # Observed taxa p/a data
# BIOvld.pa <- BIOvld
# BIOvld.pa[BIOvld.pa > 0] <- 1
# spe_data_pa <- BIOvld.pa
# env_data <- ENVvld
# # Run the run script to run all Fo vs Fe functions
# source("R/RunFoFeFunctions_ver1.2.r")

# --------------------------------------------------------------------

# SiteSetName <- "test_sites"
# capture_probs <- OERFtest$Capture.Probs
# BIOtest.pa <- BIOtest
# BIOtest.pa[BIOtest.pa > 0] <- 1
# spe_data_pa <- BIOtest.pa
# env_data <- ENVtest
# # Run the run script to run all Fo vs Fe functions
# source("R/RunFoFeFunctions_ver1.2.r")

# SiteSetName <- "Validation"
# capture_probs <- OERFvld$Capture.Probs
# BIOvld.pa <- BIOvld
# BIOvld.pa[BIOvld.pa > 0] <- 1
# spe_data_pa <- BIOvld.pa
# env_data <- ENVvld
# # Run the run script to run all Fo vs Fe functions
# source("R/RunFoFeFunctions_ver1.2.r")

# O/E-data
# OE_data <- OERFtest$OE.scores
# --------------------------------------------------------------------


# RUN SITE SETS ----------------------------------------------------------

BIOall.pa <- BIOall
BIOall.pa[BIOall.pa > 0] <- 1

BIOtest.pa <- BIOtest
BIOtest.pa[BIOtest.pa > 0] <- 1

plots <- list()

deviations <- list()

FoFeSimus <- list()

# Cross-validation data

SiteSetName <- "Cross-validation"
capture_probs <- capture_probs_cv
spe_data_pa <- spe_data_pa_cv
env_data <- env_data_cv

source("R/RunFoFeFunctions_ver1.2.r")

FoFeSimu_cv <- FoFeSimu

# =========================================================================
# Agriculture more than 10%
# =========================================================================
SiteSetName <- "Agriculture more than 10%"
# First, identify the sites that meet the new initial criteria
initial_selection <- all_sites[
  , "AGRICULT"
] > 10 &
  all_sites[, "FORLOSS0021"] < 15 &
  all_sites[, "ARTIFICI"] < 2

# Get the indices for the sites that meet the initial criteria
initial_indices <- which(initial_selection)

# Identify indices of sites from Sweden ("SE") and Finland ("FI")
# Identify indices of sites from Sweden ("SE") and Finland ("FI")
se_indices <- initial_indices[all_sites[initial_indices, "country"] == "SE"]

fi_indices <- initial_indices[all_sites[initial_indices, "country"] == "FI"]

se_fi_indices <- c(se_indices, fi_indices)

# Handle cases where there are fewer than 120 sites
num_se_fi_sites <- length(se_fi_indices)
num_to_sample <- min(120, num_se_fi_sites)

# Randomly select the sites
set.seed(42) # Optional: for reproducibility
selected_se_fi_indices <- sample(se_fi_indices, num_to_sample, replace = FALSE)

# Identify indices of sites from Norway ("NO")
no_indices <- initial_indices[all_sites[initial_indices, "country"] == "NO"]

# Combine the indices for the final selection
final_indices <- c(selected_se_fi_indices, no_indices)

# Now, use these combined indices to select your data
capture_probs <- OERFall$Capture.Probs[initial_indices, ]
spe_data_pa <- BIOall.pa[initial_indices, ]
env_data <- ENVall[initial_indices, ]
source("R/RunFoFeFunctions_ver1.2.r")

plots[[1]] <- p
deviations[[1]] <- c(increasing_percentage, decreasing_percentage)
FoFeSimus[[1]] <- FoFeSimu

# ====================================================================
# Forestry more than 20%
# ====================================================================

SiteSetName <- "Forestry more than 20%"
# First, identify the sites that meet the initial criteria
initial_selection <- all_sites[, "FORLOSS0021"] > 20 &
  all_sites[, "AGRICULT"] < 2 &
  all_sites[, "ARTIFICI"] < 1

# Get the indices for the sites that meet the initial criteria
initial_indices <- which(initial_selection)

# Identify indices of sites from Sweden ("SE") and Finland ("FI")
se_indices <- initial_indices[all_sites[initial_indices, "country"] == "SE"]

fi_indices <- initial_indices[all_sites[initial_indices, "country"] == "FI"]

se_fi_indices <- c(se_indices, fi_indices)

# Handle cases where there are fewer than 120 sites
num_se_fi_sites <- length(se_fi_indices)
num_to_sample <- min(120, num_se_fi_sites)

# Randomly select the sites
set.seed(42) # Optional: for reproducibility
selected_se_fi_indices <- sample(se_fi_indices, num_to_sample, replace = FALSE)

# Identify indices of sites from Norway ("NO")
no_indices <- initial_indices[all_sites[initial_indices, "country"] == "NO"]

# Combine the indices for the final selection
final_indices <- c(selected_se_fi_indices, no_indices)

# Now, use these combined indices to select your data
capture_probs <- OERFall$Capture.Probs[initial_indices, ]
spe_data_pa <- BIOall.pa[initial_indices, ]
env_data <- ENVall[initial_indices, ]
source("R/RunFoFeFunctions_ver1.2.r")

plots[[2]] <- p
deviations[[2]] <- c(increasing_percentage, decreasing_percentage)
FoFeSimus[[2]] <- FoFeSimu

# ====================================================================
# Urban more than 5%
# ====================================================================

SiteSetName <- "Urban more than 5%"
capture_probs <- OERFall$Capture.Probs[
  all_sites[, "ARTIFICI"] > 5 &
    all_sites[, "AGRICULT"] < 2 &
    all_sites[, "FORLOSS0021"] < 15,
]
spe_data_pa <- BIOall.pa[
  all_sites[, "ARTIFICI"] > 5 &
    all_sites[, "AGRICULT"] < 2 &
    all_sites[, "FORLOSS0021"] < 15,
]
env_data <- ENVall[
  all_sites[, "ARTIFICI"] > 5 &
    all_sites[, "AGRICULT"] < 2 &
    all_sites[, "FORLOSS0021"] < 15,
]
source("R/RunFoFeFunctions_ver1.2.r")

plots[[3]] <- p
deviations[[3]] <- c(increasing_percentage, decreasing_percentage)
FoFeSimus[[3]] <- FoFeSimu

# ====================================================================
# Multi-pressure sites
# ====================================================================

SiteSetName <- "Multiple pressures"

if (ecosystem == "lakes") {
  # Select sites based on pressure thresholds

  indices <- which(
    (all_sites[, "AGRICULT"] > 10 & all_sites[, "FORLOSS0021"] > 20) |
      (all_sites[, "AGRICULT"] > 10 & all_sites[, "ARTIFICI"] > 5) |
      (all_sites[, "FORLOSS0021"] > 20 & all_sites[, "ARTIFICI"] > 5)
  )
} else if (ecosystem == "rivers") {
  # Load previously selected river multiple-pressure sites

  load(here("data/multiple_pressure_site_ids.RData"))
  indices <- all_sites$site_ID %in% multiple_pressure_sites
} else {
  stop("Unknown ecosystem type")
}

# Apply selection

capture_probs <- OERFall$Capture.Probs[indices, ]
spe_data_pa <- BIOall.pa[indices, ]
env_data <- ENVall[indices, ]

source("R/RunFoFeFunctions_ver1.2.r")

plots[[4]] <- p
deviations[[4]] <- c(
  increasing_percentage,
  decreasing_percentage
)
FoFeSimus[[4]] <- FoFeSimu

# ==================================================================
# Generate zoomed in plots
# ==================================================================

global_ymax <- max(
  sapply(FoFeSimus, function(x) {
    max(x$Fo[x$Fe <= 0.1], na.rm = TRUE)
  })
)

zoom_plots <- lapply(
  FoFeSimus,
  PlotFoFe_zoom,
  ymax = global_ymax
)
