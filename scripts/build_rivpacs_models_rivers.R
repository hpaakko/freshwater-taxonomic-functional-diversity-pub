# ============================================================================
# 04_build_rivpacs_models.R
#
# Purpose:
# Build RIVPACS models for river reference sites using Random Forest and RIVPACS methodology.
#
# Inputs:
#   - Compiled reference site biological and environmental data (bio_data_ref_sites)
#   - All samples/sites biological and environmental data (bio_data_all_sites)
#   - RIVPACS model functions (model.predict.RanFor.4.2.r, rep.sam.sd.r)
#   - Clustered reference site data (ref_rivers_with_clusters.csv)
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

input_dir <- here("data")
output_dir <- here("output")

# ============================================================================
# 1. Load data
# ============================================================================

# Read compiled reference site biological and environmental data
bio_data_ref_sites <- fread(file.path(output_dir, "rbioREF_N392.csv"))

# All samples/sites
bio_data_all_sites <- fread(file.path(output_dir, "rbioAll_N2659.csv"))

# Check if there are multiple samples per site

multiple_samples <- bio_data_all_sites %>%
  group_by(site_ID) %>%
  summarise(sample_count = n()) %>%
  filter(sample_count > 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Identify missing values (optional: to check before removing)
bio_data_ref_sites[!complete.cases(bio_data_ref_sites[, c("altitude", "catch_area_km2")]), ]

# Remove rows with missing altitude or catchment area
bio_data_ref_sites <- bio_data_ref_sites[complete.cases(bio_data_ref_sites[, c("altitude", "catch_area_km2")]), ]

bio_data_all_sites <- bio_data_all_sites[complete.cases(bio_data_all_sites[, c("altitude", "catch_area_km2")]), ]

# Count rows where MostRecentYear is < 2000
nrow(bio_data_ref_sites[bio_data_ref_sites$MostRecentYear < 2000, ])
nrow(bio_data_all_sites[bio_data_all_sites$MostRecentYear < 2000, ])

# Remove rows where MostRecentYear is < 2000

bio_data_all_sites <- bio_data_all_sites[bio_data_all_sites$MostRecentYear >= 2000, ]
bio_data_ref_sites <- bio_data_ref_sites[bio_data_ref_sites$MostRecentYear >= 2000, ]

# ============================================================================
# 2. Species filtering
# ============================================================================

taxa_cols <- 52:289

dim(bio_data_ref_sites)
table(bio_data_ref_sites$country, bio_data_ref_sites$STATUS, useNA = "always")

names(bio_data_ref_sites)

# Extract species data
species_data <- as.data.frame(bio_data_ref_sites[, ..taxa_cols])

# Print all NAs
sum(is.na(species_data))

# turn NA into 0

species_data[is.na(species_data)] <- 0

# Remove taxa that are entirely zero

# Turn columns into numeric (if not already)
species_data[] <- lapply(species_data, function(x) as.numeric(as.character(x)))

species_data <- species_data[, colSums(species_data) > 0]

# Check for empty samples (rows that sum to zero)
table(rowSums(species_data) == 0)

multiple_pressures <- bio_data_all_sites[which(
  (bio_data_all_sites[, "AGRICULT"] > 10 & bio_data_all_sites[, "FORLOSS0021"] > 20) |
    (bio_data_all_sites[, "AGRICULT"] > 10 & bio_data_all_sites[, "ARTIFICI"] > 5) |
    (bio_data_all_sites[, "FORLOSS0021"] > 20 & bio_data_all_sites[, "ARTIFICI"] > 5)
)]

# ============================================================================
# 3. Community clustering
# ============================================================================

# Apply Hellinger transformation and log transformation
species_data_hellinger <- decostand(species_data, method = "hellinger")
species_data_log <- log10(species_data + 1)

# Compute Bray-Curtis dissimilarity matrix
dissimilarity_matrix <- vegdist(species_data_hellinger, method = "bray") # Bray-Curtis distance

# Hierarchical clustering using flexible beta linkage (Beta = -0.5, so Alpha = 0.75)
clustering_result <- agnes(dissimilarity_matrix, method = "flexible", par.method = (1 - (-0.5)) / 2)

# Plot clustering tree
pltree(clustering_result, col = "forest green", pch = 1)

# Convert the 'agnes' object to 'hclust' for visualization
hclust_obj <- as.hclust(clustering_result)

# Plot the dendrogram
pltree(clustering_result, col = "forest green", pch = 1)

# Draw cluster boxes (for 6 clusters)
rect.hclust(hclust_obj, k = 12, border = 1:12)

# ---- Manual clustering by hand ----

plot(clustering_result, which.plots = 2, labels = row.names(species_data), cex = .4)

# Next, drag to spread out the plot horizontally
# Interactive ID of clusters on dendrogram;
manual_cluster_points <- identify(as.hclust(clustering_result), N = 50, MAXCLUSTER = 50)

# Cluster selection is complete
# Convert result to a group ID vector
group.id <- rep(NA, dim(species_data)[[1]])
for (k in 1:length(manual_cluster_points)) {
  group.id[manual_cluster_points[[k]]] <- k
}

# Count number of sites in each group
table(group.id)

# List calibration sites and their group assignments;
cbind(row.names(species_data), group.id)

# ---- Automatically assign clusters ----

bio_data_ref_sites$automatic_cluster <- group.id
# cutree(clustering_result, 12)

# ---- Save clustering results ----
# Save the clustering results to a CSV file
fwrite(
  bio_data_ref_sites,
  file = file.path(here("output", "ref_rivers_with_clusters.csv")),
  sep = ";",
  bom = TRUE
)

# Read the saved data with clusters

bio_data_ref_sites <- fread(
  file = file.path(here("output", "ref_rivers_with_clusters.csv")),
  sep = ";",
  header = TRUE
)

# ---- Draw map of reference sites with cluster coloring ----

# Get map data for Norway, Sweden, and Finland

nordic_countries <- map_data("world", region = c("Norway", "Sweden", "Finland"))

ggplot(data = bio_data_ref_sites, aes(y = latitude, x = longitude, col = as.factor(automatic_cluster))) +
  geom_polygon(data = nordic_countries, aes(x = long, y = lat, group = group), fill = "lightgrey", color = "black") +
  geom_point(size = 3) +
  coord_cartesian(xlim = c(5, 32), ylim = c(54, 72)) +
  labs(
    title = "Geographic distribution of river REF sites with clusters",
    x = "Longitude",
    y = "Latitude",
    color = "Cluster"
  ) +
  theme_minimal()

# ---- Color dendrogram labels based on clusters ----

library(colorspace)
library(dendextend)
library(data.table)

# 1. Get the dendrogram
dend <- as.dendrogram(clustering_result)

# 2. Get labels and map clusters
dend_labels <- labels(dend)
cluster_mapping <- setNames(bio_data_ref_sites$automatic_cluster, bio_data_ref_sites$site_ID)

ordered_clusters <- bio_data_ref_sites$automatic_cluster[
  as.numeric(dend_labels)
]

# 3. Determine the number of clusters (k) and create colors
unique_clusters <- sort(unique(ordered_clusters))
k <- length(unique_clusters)
cluster_colors <- qualitative_hcl(k, palette = "Dark 3")
names(cluster_colors) <- unique_clusters

# --- Coloring ONLY the Labels ---

# 4. Color the labels based on your original automatic_cluster
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
legend("topright",
  legend = names(cluster_colors),
  fill = cluster_colors, # To represent the label colors
  title = "Cluster",
  bty = "n",
  cex = 0.8
)

# ============================================================================
# 4. Random Forest model
# ============================================================================

predictor_candidates_all <- c(
  "Log10_catch_area_km2", "altitude", "PEATLPoctop", "bio4_06", "bio15_06", "bio1_06", "bio12_06"
)

df_bio_data_ref_sites <- as.data.frame(bio_data_ref_sites)

rf_data <- data.frame(
  altitude = df_bio_data_ref_sites$altitude,
  Log10_catch_area_km2 = log10(df_bio_data_ref_sites$catch_area_km2 + 1),
  WATER = df_bio_data_ref_sites$WATER,
  PEATLPoctop = df_bio_data_ref_sites$PEATLPoctop,
  FORESTS = df_bio_data_ref_sites$FORESTS,
  bio1_06 = df_bio_data_ref_sites$bio1_06,
  bio2_06 = df_bio_data_ref_sites$bio2_06,
  bio3_06 = df_bio_data_ref_sites$bio3_06,
  bio4_06 = df_bio_data_ref_sites$bio4_06,
  bio5_06 = df_bio_data_ref_sites$bio5_06,
  bio6_06 = df_bio_data_ref_sites$bio6_06,
  bio7_06 = df_bio_data_ref_sites$bio7_06,
  bio8_06 = df_bio_data_ref_sites$bio8_06,
  bio9_06 = df_bio_data_ref_sites$bio9_06,
  bio10_06 = df_bio_data_ref_sites$bio10_06,
  bio11_06 = df_bio_data_ref_sites$bio11_06,
  bio12_06 = df_bio_data_ref_sites$bio12_06,
  bio13_06 = df_bio_data_ref_sites$bio13_06,
  bio14_06 = df_bio_data_ref_sites$bio14_06,
  bio15_06 = df_bio_data_ref_sites$bio15_06,
  bio16_06 = df_bio_data_ref_sites$bio16_06,
  bio17_06 = df_bio_data_ref_sites$bio17_06,
  bio18_06 = df_bio_data_ref_sites$bio18_06,
  bio19_06 = df_bio_data_ref_sites$bio19_06,
  longitude = df_bio_data_ref_sites$longitude,
  latitude = df_bio_data_ref_sites$latitude,
  Groups = factor(df_bio_data_ref_sites$automatic_cluster)
)

# Train the Random Forest model

set.seed(123)
rf_model <- randomForest(Groups ~ ., data = rf_data, ntree = 500, importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
print(rf_model)

# Check variable importance

varImpPlot(rf_model, type = 1)

# ---- Final RF model ----

# Select final predictors based on importance, ecological relevance with no correlation between predictors
# Minimal human influence on predictors important

final_predictors <- c(
  "Log10_catch_area_km2", "altitude", "PEATLPoctop", "bio1_06", "bio12_06"
)

# Remove any sites with missing values in the final predictors

rf_data <- rf_data[complete.cases(rf_data[, final_predictors]), ]

final_rf_data <- data.frame(rf_data[, final_predictors], Groups = rf_data$Groups)

# Reproducibility:
# RIVPACS O/E scores depend on RF vote probabilities.
# Keep seed fixed so model outputs remain identical between runs.
set.seed(123)
rf_mod_final <- randomForest(Groups ~ ., data = final_rf_data, ntree = 500, importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
print(rf_mod_final)
varImpPlot(rf_mod_final, type = 1)

# ---- Environmental niches ----

# Separate your environmental variables and cluster labels
env_data <- final_rf_data[, -which(names(final_rf_data) == "Groups")] # or use "Groups" if that's the name
cluster_labels <- final_rf_data$Groups # or final_rf_data$Groups

# Run PERMANOVA
adonis_result <- adonis2(env_data ~ cluster_labels, method = "euclidean")

# Print result
print(adonis_result)

env_data %>%
  mutate(Cluster = as.factor(cluster_labels)) %>%
  group_by(Cluster) %>%
  summarise(across(everything(), list(mean = mean, sd = sd)))

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# Step 1: Calculate p-values for each variable using ANOVA
anova_results <- env_data %>%
  mutate(Cluster = as.factor(cluster_labels)) %>%
  summarise(across(
    .cols = -Cluster,
    .fns = list(p = ~ summary(aov(.x ~ Cluster))[[1]][["Pr(>F)"]][1]),
    .names = "{.col}_p"
  ))

# If you’re using the randomForest package:
var_imp <- importance(rf_mod_final, type = 1)
# type = 1 gives mean decrease in accuracy; type = 2 gives mean decrease in Gini.

# Turn into a tidy data frame
library(dplyr)
imp_df <- as.data.frame(var_imp) %>%
  tibble::rownames_to_column("Variable") %>%
  rename(Importance = MeanDecreaseAccuracy) %>%
  arrange(desc(Importance))

# See the variables
top_n <- 5
top_vars <- imp_df %>%
  slice(1:top_n) %>%
  pull(Variable)
print(imp_df[1:top_n, ])


# Step 2: Prepare data for plotting
plot_data <- env_data %>%
  mutate(Cluster = as.factor(cluster_labels)) %>%
  dplyr::select(all_of(top_vars), Cluster) %>%
  pivot_longer(cols = all_of(top_vars), names_to = "Variable", values_to = "Value")


# Step 3: Create boxplots
ggplot(plot_data, aes(x = Cluster, y = Value, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  facet_wrap(~Variable, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title = paste("Top", top_n, "Environmental Variables Differentiating Clusters"),
    x = "Cluster", y = "Value"
  )

# ============================================================================
# 5. RIVPACS model build
# ============================================================================

# Load required model functions
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
  file.path(here(output_dir, "observed_abundances_rivers.csv")),
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

bio_data_ref_sites$Log10_catch_area_km2 <-
  log10(bio_data_ref_sites$catch_area_km2 + 1)

bio_data_all_sites$Log10_catch_area_km2 <-
  log10(bio_data_all_sites$catch_area_km2 + 1)

bio_data_imp_sites$Log10_catch_area_km2 <-
  log10(bio_data_imp_sites$catch_area_km2 + 1)

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
  # Make parameters numeric

  bug_cal_pa[] <- lapply(bug_cal_pa, function(x) as.numeric(as.character(x)))
  new_bugs[] <- lapply(new_bugs, function(x) as.numeric(as.character(x)))
  new_predictors[] <- lapply(new_predictors, function(x) as.numeric(as.character(x)))
  new_bugs[] <- lapply(new_bugs, function(x) as.numeric(as.character(x)))

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
  output_file <- file.path(output_dir, paste0(output_name, "_N", dim(results$OE.scores)[1], ".csv"))
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
  group_by(automatic_cluster) %>%
  mutate(fold = sample(rep(1:k, length.out = n()))) %>%
  ungroup()

# Save clustered_data and final_predictors to Rdata

temp_env <- new.env()
temp_env$clustered_data <- clustered_data
temp_env$final_predictors <- final_predictors

save(temp_env, file = file.path(output_dir, "clustered_data_and_predictors_rivers.Rdata"))

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

  labels_cal <- train_data$automatic_cluster

  # --- TRAIN THE MODEL FOR THIS FOLD ---
  # Note: You must re-train the RF model inside the loop so it hasn't seen 'valid_data'
  rf_fold <- randomForest(x = preds_cal, y = as.factor(labels_cal), ntree = 500)

  randomForest(Groups ~ ., data = final_rf_data, ntree = 500, importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)

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
write.csv(ref_freq_table, file.path(output_dir, "Rivers_Reference_Taxa_Frequency_Validation.csv"), row.names = FALSE)

# ===========================================================================
# 7. Run RIVPACS model
# ===========================================================================

# Evaluate model performance for calibration sites (Cal.OOB = TRUE)
OE_assess_cal_OB <- evaluate_OE(
  bug_cal_pa = bug_reference_pa,
  group_labels = group_labels,
  predictors = final_predictors,
  model = rf_mod_final,
  new_predictors = predictors_reference,
  new_bugs = bug_reference_pa,
  prob_threshold = prob_threshold,
  cal_oob = TRUE,
  output_name = "OEassesscalOB"
)

# Evaluate model performance for all sites (Cal.OOB = FALSE)
OE_assess_all_OB <- evaluate_OE(
  bug_cal_pa = bug_reference_pa,
  group_labels = group_labels,
  predictors = final_predictors,
  model = rf_mod_final,
  new_predictors = predictors_all,
  new_bugs = bug_all_pa,
  prob_threshold = prob_threshold,
  cal_oob = FALSE,
  output_name = "OEassessallOB"
)

# PERFORMANCE ANALYSIS

# Create a status column for visualization
bio_data_all_sites$Group <- "Other/Intermediate"

bio_data_all_sites$Group[bio_data_all_sites$STATUS == "REF"] <- "Reference"
bio_data_all_sites$Group[bio_data_all_sites$STATUS == "IMP"] <- "Impacted"

ggplot(bio_data_all_sites, aes(x = bio4_06, y = Log10_catch_area_km2, color = Group)) +
  geom_point() +
  scale_color_manual(values = c("Impacted" = "red", "Reference" = "blue", "Other" = "grey"))

ggplot(bio_data_all_sites, aes(x = bio1_06, y = altitude, color = Group)) +
  geom_point() +
  scale_color_manual(values = c("Impacted" = "red", "Reference" = "blue", "Other" = "grey"))

ggplot(bio_data_all_sites, aes(x = bio12_06, y = PEATLPoctop, color = Group)) +
  geom_point() +
  scale_color_manual(values = c("Impacted" = "red", "Reference" = "blue", "Other" = "grey"))

# Look at the mean latitude of your groups
aggregate(latitude ~ Group, data = bio_data_all_sites, mean)

# Check the range of your most important predictor
aggregate(Log10_catch_area_km2 ~ Group, data = bio_data_all_sites, range)

aggregate(altitude ~ Group, data = bio_data_all_sites, range)
aggregate(bio4_06 ~ Group, data = bio_data_all_sites, range)
aggregate(bio1_06 ~ Group, data = bio_data_all_sites, range)
aggregate(bio12_06 ~ Group, data = bio_data_all_sites, range)

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
  file = file.path(output_dir, paste0("rbioREF_model_N", nrow(reference_data_2), ".csv")),
  sep = ";",
  bom = TRUE
)

bio_data_all_sites_2 <- bio_data_all_sites[!is.na(catch_area_km2) & !is.na(altitude), ]

bio_data_all_sites_2$country <- factor(bio_data_all_sites_2$country, levels = c("NO", "SE", "FI"))

bio_data_all_sites_2$O <- OE_assess_all_OB$OE.scores$O
bio_data_all_sites_2$E <- OE_assess_all_OB$OE.scores$E
bio_data_all_sites_2$OoverE <- OE_assess_all_OB$OE.scores$OoverE

save(bio_data_all_sites_2, file = file.path(output_dir, "bio_data_all_sites_with_OE_rivers.Rdata"))

# ============================================================================
# 8. Export results
# ============================================================================

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
save(list = ls(envir = temp_env), envir = temp_env, file = file.path(output_dir, "Strict_Henna_RIVPACS_data_rivers.RData"))

# Verify by clearing environment and loading
rm(list = ls())
load(file.path(here("output/Strict_Henna_RIVPACS_data_rivers.RData")))
print("Variables in environment after loading new RData:")
print(ls())
