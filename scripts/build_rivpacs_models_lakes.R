# ============================================================================
# 04_build_rivpacs_models_lakes.R
#
# Purpose:
# Build RIVPACS models for lake reference sites using Random Forest and RIVPACS methodology.
#
# Inputs:
#   - Compiled reference site biological and environmental data (bio_data_ref_sites)
#   - All samples/sites biological and environmental data (bio_data_all_sites)
#   - RIVPACS model functions (model.predict.RanFor.4.2.r, rep.sam.sd.r)
#   - Clustered reference site data (ref_lakes_with_clusters.csv)
#
# Outputs:
#   - RData file for FoFe analysis
#   - Cross-validation results for RIVPACS model
#
# Manuscript:
# Differential land use effects on taxonomic and functional
# freshwater biodiversity in Nordic streams and lakes
# =============================================================================

# Download libraries

library(vegan)
library(cluster)
library(randomForest)
library(data.table)
library(ggplot2)
library(dplyr)
library(here)
library(tidyr)

input_dir <- here("data")
output_dir <- here("output")

# ============================================================================
# 1. Load data
# ============================================================================

bio_data_ref_sites <- fread(here("output/lbioREF_N207.csv"))

bio_data_all_sites <- fread(here("output/lbioAll_N817.csv"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lake_core_predictors <- c("altitude", "lake_area_km2", "catch_area_km2", "bio1_06", "bio12_06", "PEATLPoctop")

bio_data_ref_sites[!complete.cases(bio_data_ref_sites[, ..lake_core_predictors]), ]
missing_values <- bio_data_all_sites[!complete.cases(bio_data_all_sites[, ..lake_core_predictors]), ]

bio_data_all_sites <- bio_data_all_sites[complete.cases(bio_data_all_sites[, ..lake_core_predictors]), ]

# ============================================================================
# 2. Species filtering
# ============================================================================

taxa_cols <- 51:226

names(bio_data_ref_sites)

# Fine level: 51:253, coarse level: 51:226
species_data <- as.data.frame(bio_data_ref_sites[, ..taxa_cols])

# Remove taxa that are entirely zero

table(colSums(species_data) > 0)
species_data <- species_data[, colSums(species_data) > 0]
dim(species_data) # Check dimensions after filtering

# ============================================================================
# 3. Community clustering
# ============================================================================

# Apply Hellinger transformation to species data
species_data_hellinger <- decostand(species_data, method = "hellinger")
species_data_log <- log10(species_data + 1)

# Compute Bray-Curtis dissimilarity matrix
dissimilarity_matrix <- vegdist(species_data_hellinger, method = "bray")

# Perform hierarchical clustering using Ward's method
clustering_result <- agnes(dissimilarity_matrix, method = "ward")

# Hierarchical clustering using flexible beta linkage (Beta = -0.5, so Alpha = 0.75)
# clustering_result <- agnes(dissimilarity_matrix, method = "flexible", par.method = (1 - (-0.5)) / 2)

# Plot dendrogram to visualize clustering

pltree(clustering_result,
  main = "Dendrogram of Reference Sites Clustering",
  xlab = "Reference Sites", ylab = "Height", hang = -1
)

# Convert agnes to hclust object for cutting

hclust_result <- as.hclust(clustering_result)

# Plot the dendogram

plot(hclust_result,
  main = "Dendrogram of Reference Sites Clustering",
  xlab = "Reference Sites", ylab = "Height", hang = -1
)

# Draw cluster boxes (for 10 clusters)

rect.hclust(hclust_result, k = 6, border = "red")

# ---- Manual clustering by hand ----
#
# plot(clustering_result,
#   which.plots = 2, main = "Dendrogram of Reference Sites Clustering",
#   xlab = "Reference Sites", ylab = "Height", hang = -1
# )
#
# manual_cluster_points <- identify(as.hclust(clustering_result), N = 50, MAXCLUSTER = 50)
#
# group.id <- rep(NA, dim(species_data)[[1]])
# for (k in 1:length(manual_cluster_points)) {
#   group.id[manual_cluster_points[[k]]] <- k
# }
#
# table(group.id)
#
# cbind(row.names(species_data), group.id)
#
# bio_data_ref_sites$cluster <- group.id
#
# # Save clustering result for validation
#
# fwrite(
#   bio_data_ref_sites,
#   file = here("output", paste0("lbioREF_N", nrow(bio_data_ref_sites), "_clustering.csv")),
#   sep = ";",
#   bom = TRUE
# )

# Read the saved data with clusters

bio_data_ref_sites <- fread(file.path(input_dir, "lbioREF_N207_clustering.csv"))

# Draw map of reference sites with cluster coloring

# Get map data for Norway, Sweden, and Finland

nordic_countries <- map_data("world", region = c("Norway", "Sweden", "Finland"))

ggplot(data = bio_data_ref_sites, aes(y = latitude, x = longitude, col = as.factor(cluster))) +
  geom_polygon(data = nordic_countries, aes(x = long, y = lat, group = group), fill = "lightgrey", color = "black") +
  geom_point(size = 3) +
  coord_cartesian(xlim = c(5, 32), ylim = c(54, 72)) +
  labs(
    title = "Geographic distribution of lake REF sites with clusters",
    x = "Longitude",
    y = "Latitude",
    color = "Cluster"
  ) +
  theme_minimal()

# ---- Color dendrogram labels based on clusters ----

library(colorspace)
library(dendextend)

# 1. Get the dendrogram
dend <- as.dendrogram(clustering_result)

# 2. Get labels and map clusters
dend_labels <- labels(dend)
cluster_mapping <- setNames(bio_data_ref_sites$cluster, bio_data_ref_sites$site_ID)
ordered_clusters <- cluster_mapping[dend_labels] # Clusters in dendrogram order

# 3. Determine the number of clusters (k) and create colors
unique_clusters <- sort(unique(ordered_clusters))
k <- length(unique_clusters)
cluster_colors <- qualitative_hcl(k, palette = "Dark 3")
names(cluster_colors) <- unique_clusters

# --- Coloring ONLY the Labels ---

# 4. Color the labels based on your original cluster
label_colors <- cluster_colors[as.character(ordered_clusters)]
labels_colors(dend) <- label_colors

# Optional: Adjust label size or other label attributes
dend <- set(dend, "labels_cex", 0.8) # Adjust text size

# 5. Plot the dendrogram (branches will be default color)
plot(dend,
  main = "Dendrogram with Colored Labels",
  xlab = "Sites",
  ylab = "Distance"
) # Add axis labels for clarity


# 7. Add legend (reflecting label colors)
# legend("topright",
#   legend = names(cluster_colors),
#   fill = cluster_colors, # To represent the label colors
#   title = "Cluster",
#   bty = "n",
#   cex = 0.8
# )


# ============================================================================
# 4. Random Forest model
# ============================================================================

predictor_candidates_all <- c(
  "Log10_Lake_Area_km2", "Log10_catch_area_km2", "altitude", "bio1_06", "bio4_06", "bio12_06", "bio15_06", "WATER"
)

df_bio_data_ref_sites <- as.data.frame(bio_data_ref_sites)

rf_data <- data.frame(
  Log10_Lake_Area_km2 = log10(df_bio_data_ref_sites$lake_area_km2 + 1),
  Log10_catch_area_km2 = log10(df_bio_data_ref_sites$catch_area_km2 + 1),
  altitude = df_bio_data_ref_sites$altitude,
  bio1_06 = df_bio_data_ref_sites$bio1_06,
  bio4_06 = df_bio_data_ref_sites$bio4_06,
  bio12_06 = df_bio_data_ref_sites$bio12_06,
  bio15_06 = df_bio_data_ref_sites$bio15_06,
  PEATLPoctop = df_bio_data_ref_sites$PEATLPoctop,
  WATER = df_bio_data_ref_sites$WATER,
  Groups = factor(df_bio_data_ref_sites$cluster)
)

dim(rf_data) # Check dimensions of the RF data
str(rf_data) # Check structure of the RF data
head(rf_data) # View the first few rows of the RF data
table(rf_data$Groups) # Check the distribution of groups

# Train initial RF model

set.seed(123) # For reproducibility
rf_model <- randomForest(Groups ~ ., data = rf_data, ntree = 500, importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
print(rf_model)
varImpPlot(rf_model, type = 1)

# ---- Final RF model ----

# Select final predictors based on importance plot
# IMPORTANT: Choose the most influential variables for lakes.

final_predictors <- c(
  "Log10_Lake_Area_km2", "altitude", "bio4_06", "bio1_06"
)

final_rf_data <- data.frame(rf_data[, final_predictors], Groups = rf_data$Groups)

set.seed(123) # For reproducibility
rf_mod_final <- randomForest(Groups ~ ., data = final_rf_data, ntree = 500, importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
print(rf_mod_final)
varImpPlot(rf_mod_final, type = 1)

# ---- Environmental niches ----

env_data <- final_rf_data[, -which(names(final_rf_data) == "Groups")]
cluster_labels <- final_rf_data$Groups
adonis2(env_data ~ cluster_labels, method = "euclidean")

# Boxplot code using ggplot2 for top variables

top_vars <- final_predictors

plot_data <- env_data %>%
  mutate(cluster = as.factor(cluster_labels)) %>%
  dplyr::select(all_of(top_vars), cluster) %>%
  pivot_longer(cols = all_of(top_vars), names_to = "Variable", values_to = "Value")

# Create boxplots

ggplot(plot_data, aes(x = cluster, y = Value, fill = cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  facet_wrap(~Variable, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title = paste("Top Environmental Variables Differentiating Clusters"),
    x = "Cluster", y = "Value"
  )

# ============================================================================
# 5. RIVPACS model build
# ============================================================================

# Source RIVPACS model functions (ensure paths are correct)
source(here("scripts/RIVPACS model/model.predict.RanFor.4.2.r"))
source(here("scripts/RIVPACS model/rep.sam.sd.r")) # occprb.cal

# ============================================================================
# Prepare data for RIVPACS model
# ============================================================================

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

prob_threshold <- 0.5

group_labels <- final_rf_data$Groups

# ----------------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------------

to_pa <- function(x) {
  x <- as.data.frame(x)
  x[x > 0] <- 1
  x
}

prepare_predictor_data <- function(df, predictors) {
  df <- as.data.frame(df)

  df[
    complete.cases(df[, predictors]),
  ]
}

# ----------------------------------------------------------------------------
# Biological data
# ----------------------------------------------------------------------------

# Reference sites

bug_reference <- as.data.frame(
  bio_data_ref_sites[, ..taxa_cols]
)

bug_reference_pa <- to_pa(
  bug_reference
)

# All sites

bug_all <- as.data.frame(
  bio_data_all_sites[, ..taxa_cols]
)

fwrite(
  bug_all,
  file.path(here(output_dir, "observed_abundances_lakes.csv")),
  sep = ";",
  bom = TRUE
)

bug_all_pa <- to_pa(
  bug_all
)

# Impacted sites

bio_data_imp_sites <- bio_data_all_sites[
  STATUS == "IMP"
]

bio_data_imp_sites_df <- as.data.frame(
  bio_data_imp_sites
)

bug_test <- as.data.frame(
  bio_data_imp_sites[, ..taxa_cols]
)

bug_test_pa <- to_pa(
  bug_test
)

# Calibration data
# (same sites as the reference set)

bug_calibration <- bug_reference

bug_calibration_pa <- bug_reference_pa

# ----------------------------------------------------------------------------
# Predictor data
# ----------------------------------------------------------------------------

bio_data_ref_sites$Log10_Lake_Area_km2 <-
  log10(bio_data_ref_sites$lake_area_km2 + 1)

bio_data_all_sites$Log10_Lake_Area_km2 <-
  log10(bio_data_all_sites$lake_area_km2 + 1)

bio_data_imp_sites$Log10_Lake_Area_km2 <-
  log10(bio_data_imp_sites$lake_area_km2 + 1)

reference_data_df <- prepare_predictor_data(
  bio_data_ref_sites,
  final_predictors
)

bio_data_all_df <- prepare_predictor_data(
  bio_data_all_sites,
  final_predictors
)

bio_data_imp_sites_df <- prepare_predictor_data(
  bio_data_imp_sites,
  final_predictors
)

# ----------------------------------------------------------------------------
# Predictor matrices
# ----------------------------------------------------------------------------

predictors_reference <-
  reference_data_df[, final_predictors]

predictors_all <-
  bio_data_all_df[, final_predictors]

predictors_test <-
  bio_data_imp_sites_df[, final_predictors]

# ----------------------------------------------------------------------------
# Quick checks
# ----------------------------------------------------------------------------

cat(
  "Reference sites:", nrow(bug_reference),
  "\nAll sites:", nrow(bug_all),
  "\nImpacted sites:", nrow(bug_test),
  "\n"
)

# Function to evaluate OE scores and save results

evaluate_OE <- function(bug_cal_pa, group_labels, predictors, model, new_predictors, new_bugs, prob_threshold, cal_oob, output_name) {
  results <- model.predict.RanFor.4.2(
    bugcal.pa = bug_cal_pa,
    grps.final = group_labels,
    preds.final = predictors,
    ranfor.mod = model,
    prednew = new_predictors,
    bugnew = new_bugs,
    Pc = prob_threshold,
    Cal.OOB = cal_oob
  )

  # Calculate sample standard deviation for assessment
  rep_sample_sd <- rep.sam.sd(occprb.cal = results$Capture.Probs, Pc = prob_threshold)
  print(rep_sample_sd)

  # Save results
  output_file <- here("output", paste0(output_name, "_N", dim(results$OE.scores)[1], ".csv"))
  write.table(results$OE.scores, sep = ";", file = output_file)

  write.table(
    cbind(MASTERID_kk_vuosi = row.names(results$OE.scores), results$OE.scores),
    dec = ".", sep = ";", file = output_file, row.names = FALSE
  )

  # Plot observed vs. expected richness
  plot(
    results$OE.scores[, c("E", "O")],
    xlab = "Expected richness", ylab = "Observed richness",
    asp = 1, ylim = c(0, 30), xlim = c(0, 30),
    main = paste(output_name, "pt=", prob_threshold, "OOB=", cal_oob)
  )
  abline(0, 1, lty = 2)

  return(results)
}


# ============================================================================
# 6. Cross-validation
# ============================================================================

# 1. Setup CV Parameters
k <- 5
set.seed(123) # For reproducibility
clustered_data <- bio_data_ref_sites

# Create a 'fold' column: balanced across your biological clusters
clustered_data <- clustered_data %>%
  group_by(cluster) %>%
  mutate(fold = sample(rep(1:k, length.out = n()))) %>%
  ungroup()

# Save clustered_data and final_predictors to Rdata

temp_env <- new.env()
temp_env$clustered_data <- clustered_data
temp_env$final_predictors <- final_predictors

save(temp_env, file = here("output", "clustered_data_and_predictors_lakes.Rdata"))

# 1. Initialize lists to store objects from each fold
fold_results_list <- list()
bug_val_pa_list <- list()
env_val_list <- list()

# 2. Start the CV Loop
for (i in 1:k) {
  message(paste("--- Running Fold", i, "of", k, "---"))

  # Split data for this fold
  train_data <- clustered_data %>% filter(fold != i)
  valid_data <- clustered_data %>% filter(fold == i)

  # Prepare Biological Data
  bug_cal_pa <- to_pa(
    train_data[, taxa_cols]
  )

  bug_val_pa <- to_pa(
    valid_data[, taxa_cols]
  )


  # Prepare Predictors
  preds_cal <- as.data.frame(train_data[, final_predictors])
  preds_val <- as.data.frame(valid_data[, final_predictors])

  labels_cal <- train_data$cluster

  # --- TRAIN THE MODEL FOR THIS FOLD ---
  # Note: You must re-train the RF model inside the loop so it hasn't seen 'valid_data'
  # rf_fold <- randomForest(x = preds_cal, y = as.factor(labels_cal), ...)
  # Assuming your existing rf_mod_final training code goes here:
  rf_fold <- randomForest(x = preds_cal, y = as.factor(labels_cal), ntree = 500)

  # Run the assessment
  fold_results <- evaluate_OE(
    bug_cal_pa = bug_cal_pa,
    group_labels = labels_cal,
    predictors = final_predictors,
    model = rf_fold,
    new_predictors = preds_val,
    new_bugs = bug_val_pa,
    prob_threshold = prob_threshold,
    cal_oob = FALSE,
    output_name = paste0("CV_Fold_", i)
  )

  # 2. STORE RESULTS CORRECTLY
  # Use double brackets to store the entire list object
  fold_results_list[[i]] <- fold_results
  # Store the actual observed data for this fold to match the predictions
  bug_val_pa_list[[i]] <- bug_val_pa
  # STORE THE ENV DATA
  env_val_list[[i]] <- preds_val
}

# 1. Combine the O/E scores from all folds into one data frame
# We reach into each list element, grab the 'OE.scores' data frame, and stack them
all_oe_scores_combined <- do.call(rbind, lapply(fold_results_list, function(x) x$OE.scores))

# 2. Calculate the means and SDs from the combined columns
final_mean_oe <- mean(all_oe_scores_combined$OoverE, na.rm = TRUE)
final_sd_oe <- sd(all_oe_scores_combined$OoverE, na.rm = TRUE)

# 3. Repeat for the Null model scores
final_mean_null <- mean(all_oe_scores_combined$OoverE.null, na.rm = TRUE)
final_sd_null <- sd(all_oe_scores_combined$OoverE.null, na.rm = TRUE)

# 4. Print the results
message("--- FINAL CROSS-VALIDATION STATISTICS ---")
print(paste("Global Mean O/E:", round(final_mean_oe, 3)))
print(paste("Global SD O/E:", round(final_sd_oe, 3)))
print(paste("Global Mean O/E Null:", round(final_mean_null, 3)))
print(paste("Global SD O/E Null:", round(final_sd_null, 3)))

mean(abs(all_oe_scores_combined$OoverE - 1),
  na.rm = TRUE
)

mean(abs(all_oe_scores_combined$OoverE.null - 1),
  na.rm = TRUE
)

# --- Additional Analysis: Reference Taxa Frequency Validation ---

# 1. Prepare to collect capture probabilities from the CV loop
# You'll need to modify your loop slightly to store 'Capture.Probs'
# or use this logic if you already have them stored.

# Create a container for all predicted probabilities
all_cv_probs <- data.frame()
all_cv_obs <- data.frame()
all_cv_env <- data.frame()

# (Assuming you run this inside/after your loop)
# Note: Extracting from 'fold_results' in your previous code
for (i in 1:k) {
  # 'Capture.Probs' are the Pi values for each taxon at each site
  # 'bug_val_pa' are the actual observations for those same sites
  all_cv_probs <- rbind(all_cv_probs, fold_results_list[[i]]$Capture.Probs)
  all_cv_obs <- rbind(all_cv_obs, bug_val_pa_list[[i]])
  all_cv_env <- rbind(all_cv_env, env_val_list[[i]])
}

# 2. Calculate E (Expected) and O (Observed) per taxon
# E is the sum of probabilities across all sites
E_freq <- colSums(all_cv_probs, na.rm = TRUE)
# O is the sum of actual presences
O_freq <- colSums(all_cv_obs, na.rm = TRUE)

# 3. Calculate the 95% Confidence Interval for the Reference Group
# We use the binomial variance property: Var = sum(Pi * (1 - Pi))
variance_ref <- colSums(all_cv_probs * (1 - all_cv_probs), na.rm = TRUE)
sd_ref <- sqrt(variance_ref)

# 95% CI calculation (1.96 * SD)
Lower_CI_95 <- E_freq - (1.96 * sd_ref)
Upper_CI_95 <- E_freq + (1.96 * sd_ref)

# 4. Build the Reference Frequency Table
ref_freq_table <- data.frame(
  Taxon = names(E_freq),
  Observed_Freq = O_freq,
  Expected_Freq = E_freq,
  Lower_CI = Lower_CI_95,
  Upper_CI = Upper_CI_95
)

# 5. Flag taxa where the model fails to predict frequency in Reference
ref_freq_table$Model_Reliable <- ifelse(
  ref_freq_table$Observed_Freq >= ref_freq_table$Lower_CI &
    ref_freq_table$Observed_Freq <= ref_freq_table$Upper_CI,
  "Reliable", "Unreliable"
)

# Export to CSV
write.csv(ref_freq_table, file.path(output_dir, "Lakes_Reference_Taxa_Frequency_Validation.csv"), row.names = FALSE)

# ===========================================================================
# 7. Run RIVPACS model
# ===========================================================================

OE_assess_cal_OB <- evaluate_OE(
  bug_cal_pa = bug_calibration_pa,
  group_labels = group_labels,
  predictors = final_predictors,
  model = rf_mod_final,
  new_predictors = final_rf_data,
  new_bugs = bug_calibration_pa,
  prob_threshold = prob_threshold,
  cal_oob = TRUE,
  output_name = "OEassesscalOB"
)

# Evaluate model performance for all sites (Cal.OOB = FALSE)
OE_assess_all_OB <- evaluate_OE(
  bug_cal_pa = bug_calibration_pa,
  group_labels = group_labels,
  predictors = final_predictors,
  model = rf_mod_final,
  new_predictors = predictors_all,
  new_bugs = bug_all_pa,
  prob_threshold = prob_threshold,
  cal_oob = FALSE,
  output_name = "OEassessallOB"
)

# ---- Observed (O) vs Expected (E) Data Preparation ----

# Create a copy of reference site data
reference_data_2 <- bio_data_ref_sites

# Add O/E scores
reference_data_2[, `:=`(
  O = OE_assess_cal_OB$OE.scores$O,
  E = OE_assess_cal_OB$OE.scores$E,
  OoverE = OE_assess_cal_OB$OE.scores$OoverE
)]

# Save the updated dataset
fwrite(
  reference_data_2,
  file = file.path(output_dir, paste0("rbioREF_model_N", nrow(reference_data_2), ".csv")), sep = ";",
  bom = TRUE
)

bio_data_all_sites_2 <- bio_data_all_sites[!is.na(catch_area_km2) & !is.na(altitude), ]

length(OE_assess_all_OB$OE.scores$O)
bio_data_all_sites_2$country <- factor(bio_data_all_sites_2$country, levels = c("NO", "SE", "FI"))

bio_data_all_sites_2$O <- OE_assess_all_OB$OE.scores$O
bio_data_all_sites_2$E <- OE_assess_all_OB$OE.scores$E
bio_data_all_sites_2$OoverE <- OE_assess_all_OB$OE.scores$OoverE

save(bio_data_all_sites_2, file = here("output", "bio_data_all_sites_with_OE_lakes.Rdata"))

reference_data_2$OoverEclass <- cut(reference_data_2$OoverE, 5)
reference_data_2$log10ABU <- log10(reference_data_2$ABU + 1)
reference_data_2$log10ABUcat <- cut(reference_data_2$log10ABU, 4)

# ============================================================================
# 8. Export results
# ============================================================================

library(data.table)

bio_data_imp_sites <- as.data.table(bio_data_imp_sites_df)

BIO.cal.pa <- bug_calibration_pa
BIOcal <- bug_calibration
grps.finalAB2 <- group_labels
preds.finalAB2 <- final_predictors
rf.modAB <- rf_mod_final
ENVtest <- predictors_test
BIOtest <- bug_test
BIOtest.pa <- bug_test_pa
ENVall <- predictors_all
BIOall <- bug_all
BIOall.pa <- bug_all_pa
imp_data <- bio_data_imp_sites
all_sites <- bio_data_all_sites
capture_probs_cv <- all_cv_probs # From your CV loop
spe_data_pa_cv <- all_cv_obs # From your CV loop
env_data_cv <- all_cv_env # From your CV loop

# 1. Define the names of variables to be saved
vars_to_save <- c(
  "BIO.cal.pa",
  "BIOcal",
  "grps.finalAB2",
  "preds.finalAB2",
  "rf.modAB",
  "ENVtest",
  "BIOtest",
  "BIOtest.pa",
  "ENVall",
  "BIOall",
  "BIOall.pa",
  "imp_data",
  "all_sites",
  "capture_probs_cv",
  "spe_data_pa_cv",
  "env_data_cv"
)

# 2. Create a new, empty environment
temp_env <- new.env()

# 3. Copy only the desired variables from the global environment to the new environment
for (var_name in vars_to_save) {
  if (exists(var_name, envir = .GlobalEnv)) {
    assign(var_name, get(var_name, envir = .GlobalEnv), envir = temp_env)
  } else {
    warning(paste0("Variable '", var_name, "' not found in global environment. Skipping."))
  }
}

# Ensure the 'data' directory exists
if (!dir.exists("data")) {
  dir.create("data")
}

# 4. Save the new environment to the .RData file
# When saving an environment, all objects within it are saved.
save(list = ls(envir = temp_env), envir = temp_env, file = here("output", "Strict_Henna_RIVPACS_data_lakes.RData"))

# Verify by clearing environment and loading
rm(list = ls())
load(here("output", "Strict_Henna_RIVPACS_data_lakes.RData"))
print("Variables in environment after loading new RData:")
print(ls())
