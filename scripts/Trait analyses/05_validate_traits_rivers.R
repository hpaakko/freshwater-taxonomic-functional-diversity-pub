# ===========================================================================
# 05_validate_traits_rivers.R
#
# Purpose:
#   - Validate the trait abundance predictions using K-fold cross-validation on RIVPACS data.
#
# Input:
#   - clustered_data_and_predictors_rivers.Rdata
#   - traits_for_trait_abundance.RData
#
# Output:
#   - trait_validation_results_rivers.Rdata
#   - reliable_traits_rivers.Rdata
# ===========================================================================

# K-fold cross validation of functional trait prediction using RIVPACS data
# Premise: the abundance of functional trait is predicted by the abundance of taxa, and the abundance of taxa is predicted by the RIVPACS model
# The validation is done by comparing the predicted functional trait abundance with the observed functional trait abundance in the test set
# Each reference site will be used as validation once, as the reference set is divided into 5 folds, and each fold is used as validation once

library(data.table)
library(randomForest)
library(dplyr)
library(ggplot2)
library(vegan)
library(tidyr)
library(here)
library(stringr)

# ---- Cross-validation of Trait Abundance Model ----

source(here("scripts/Trait analyses/model.predict.RanFor.4.2_and_abuOE.r"))
source(here("functions/create_trait_abundance_matrix.r"))
source(here("functions/standardize_trait_abundances.r"))

load(here("output/clustered_data_and_predictors_rivers.Rdata"))

clustered_data <- temp_env$clustered_data
final_predictors <- temp_env$final_predictors

load(here("output/traits_for_trait_abundance.RData")) # Loads traits

species_trait_df <- temp_env$traits

k <- 5
set.seed(123)
# Assuming clustered_data and final_predictors are already defined from your RIVPACS setup

# Initialize storage for trait-specific results
trait_cv_results <- list()

for (i in 1:k) {
  message(paste("--- Running Fold", i, "of", k, "---"))

  # 1. Split Data
  train_data <- clustered_data %>% filter(fold != i)
  valid_data <- clustered_data %>% filter(fold == i)

  # 2. Prepare Inputs for the RIVPACS-style abundance prediction
  bugcal_pa <- as.data.frame(train_data[, 52:289])
  bugcal_pa[bugcal_pa > 0] <- 1

  bugcal_abu <- as.data.frame(train_data[, 52:289]) # Actual abundances for training

  # Valid data (Observed)
  bug_val_abu <- as.data.frame(valid_data[, 52:289])
  bug_val_pa <- bug_val_abu
  bug_val_pa[bug_val_pa > 0] <- 1

  # Predictors
  preds_cal <- as.data.frame(train_data[, final_predictors])
  preds_val <- as.data.frame(valid_data[, final_predictors])
  labels_cal <- train_data$automatic_cluster

  # 3. Train Fold Model
  rf_fold <- randomForest(x = preds_cal, y = as.factor(labels_cal), ntree = 500)

  # 4. Predict Taxon Abundances for Validation Sites
  # Using your sourced function
  pred_out <- model.predict.RanFor.4.2_abu(
    bugcal.abu = bugcal_abu,
    bugcal.pa = bugcal_pa,
    grps.final = labels_cal,
    preds.final = final_predictors,
    ranfor.mod = rf_fold,
    prednew = preds_val,
    bugnew = bug_val_abu,
    Pc = 0.5,
    Cal.OOB = FALSE
  )

  # 1. Predict Taxon Abundances (as before)
  pred_taxa_abu <- as.data.frame(pred_out$Exp.LogAbus)
  obs_taxa_abu <- as.data.frame(bug_val_abu)

  # Log-transform the observed abundances to match the scale
  obs_taxa_abu <- log10(as.data.frame(bug_val_abu) + 1)

  # 2. Convert to raw Trait Abundances
  pred_trait_raw <- create_trait_abundance_matrix(pred_taxa_abu, species_trait_df)
  obs_trait_raw <- create_trait_abundance_matrix(obs_taxa_abu, species_trait_df)

  # 3. STANDARDIZE (The New Step)
  # This ensures we are comparing predicted proportions vs observed proportions
  pred_trait_std <- standardize_trait_abundances(pred_trait_raw)
  obs_trait_std <- standardize_trait_abundances(obs_trait_raw)

  # 4. Store these for the final reliability table
  trait_cv_results[[i]] <- list(obs = obs_trait_std, pred = pred_trait_std)
}

all_trait_obs <- do.call(rbind, lapply(trait_cv_results, function(x) x$obs))
all_trait_pred <- do.call(rbind, lapply(trait_cv_results, function(x) x$pred))

# =========================================================
# TRAIT PREDICTABILITY: LM on O/E ratio if it equals 1
# =========================================================

# Calculate epsilon based on half of the smallest non-zero value in the observed trait data to avoid division by zero

trait_epsilons <- apply(all_trait_obs, 2, function(x) {
  min_val <- min(x[x > 0], na.rm = TRUE)
  # Use half of the smallest observed value as the epsilon
  return(min_val / 2)
})

# 1. Calculate ratios O/E for each trait and each site

OE_ratios <- all_trait_obs # Start with raw

# Calculate SD of the O/E ratios


for (trait in colnames(all_trait_obs)) {
  eps <- trait_epsilons[trait]
  OE_ratios[, trait] <- (all_trait_obs[, trait] + eps) / (all_trait_pred[, trait] + eps)
}

trait_sds <- apply(OE_ratios, 2, sd, na.rm = TRUE)

mean(trait_sds, na.rm = TRUE)
median(trait_sds, na.rm = TRUE)

summary(trait_sds)

# Calculate log(O/E) for better statistical properties (normality, homoscedasticity)

trait_names <- colnames(all_trait_obs)

trait_results <- data.frame(
  Trait = trait_names,
  Mean_Bias = NA,
  MAE = NA, # Mean Absolute Error (overall accuracy)
  t_statistic = NA,
  p_value = NA
)

for (j in 1:length(trait_names)) {
  ratios <- OE_ratios[, j]

  # Remove NAs for the check
  ratios_clean <- ratios[!is.na(ratios)]

  # Mean Bias
  trait_results$Mean_Bias[j] <- mean(ratios, na.rm = TRUE)

  # Corrected MAE: How far from 1.0 on average?
  trait_results$MAE[j] <- mean(abs(ratios - 1), na.rm = TRUE)

  # Check if we have enough data and if it's not constant
  if (length(ratios_clean) > 2 && var(ratios_clean) > 0) {
    test_result <- t.test(ratios, mu = 1)
    trait_results$t_statistic[j] <- test_result$statistic
    trait_results$p_value[j] <- test_result$p.value
  } else {
    # If constant, it's either perfectly 1.0 or systematically 0
    trait_results$t_statistic[j] <- NA
    trait_results$p_value[j] <- NA
  }
}

trait_results <- trait_results %>%
  mutate(
    Status = case_when(
      # 1. Handle the absent traits first
      is.na(p_value) ~ "Zero Variance (Absent)",

      # 2. Check P-value: If the T-test says 'No Difference', it's the gold standard
      p_value > 0.05 ~ "Unbiased",

      # # 3. Check Magnitude: If P < 0.05 but the error is tiny, it's 'Practically Unbiased'
      # abs(Mean_Bias - 1) < 0.1 ~ "Negligible Bias",

      # 4. Everything else is truly biased
      TRUE ~ "Biased"
    ),

    # Define Reliability: Both Category 2 and 3 are safe to use for your results
    Reliable = ifelse(Status %in% c("Unbiased"), TRUE, FALSE)
  )

print(trait_results)

# 1. Clean, Split, and Reorder the data
plot_data <- trait_results %>%
  filter(!is.na(Mean_Bias)) %>%
  # Split 'Feedinghabits_predator' into 'Major' and 'Minor'
  separate(Trait,
    into = c("Major_Group", "Minor_Group"), sep = "_",
    extra = "merge", remove = FALSE
  ) %>%
  mutate(
    # A. Clean the Major Group (Facet Labels)
    Major_Group = case_when(
      Major_Group == "Maxsize" ~ "Max body size",
      Major_Group == "Lifecycleduration" ~ "Life cycle duration",
      Major_Group == "Cycles" ~ "Cycles per year",
      Major_Group == "Aquaticstages" ~ "Aquatic stage",
      Major_Group == "Locosubstrate" ~ "Locomotion - substrate relation",
      Major_Group == "Feedinghabits" ~ "Feeding habit",
      Major_Group == "Currentvelocity" ~ "Current velocity preference",
      Major_Group == "Substrate" ~ "Substrate preference",
      Major_Group == "Temperature" ~ "Temperature preference",
      Major_Group == "pH" ~ "pH preference",
      TRUE ~ str_to_title(Major_Group)
    ),

    # B. Clean the Minor Group (X-axis Labels)
    Minor_Group_Clean = Minor_Group %>%
      str_replace_all("_", " ") %>%
      str_replace("less than", "<") %>%
      str_replace("more than", ">") %>%
      str_replace("equals", "=") %>%
      str_replace_all("\\bto\\b", "–") %>%
      str_replace("flags boulders cobbles pebbles", "flags/boulders/cobbles/pebbles") %>%
      str_replace("organic detritus litter", "Organic detritus/litter") %>%
      str_replace("twigs roots", "Twigs/roots") %>%
      str_to_sentence()
  )

plot_data <- plot_data %>%
  mutate(
    # 1. Create a hidden numeric sorting value for natural sequences
    Sort_Value = case_when(
      # pH sequence
      Minor_Group_Clean == "< 4" ~ 1,
      Minor_Group_Clean == "4 – 4.5" ~ 2,
      Minor_Group_Clean == "4.5 – 5" ~ 3,
      Minor_Group_Clean == "5 – 5.5" ~ 4,
      Minor_Group_Clean == "5.5 – 6" ~ 5,
      Minor_Group_Clean == "> 6" ~ 6,

      # Size sequence
      Minor_Group_Clean == "< 0.25" ~ 1,
      Minor_Group_Clean == "0.25 – 0.5" ~ 2,
      Minor_Group_Clean == "0.5 – 1" ~ 3,
      Minor_Group_Clean == "1 – 2" ~ 4,
      Minor_Group_Clean == "2 – 4" ~ 5,
      Minor_Group_Clean == "4 – 8" ~ 6,
      Minor_Group_Clean == "> 8" ~ 7,

      # Cycles per year
      Minor_Group_Clean == "< 1" ~ 1,
      Minor_Group_Clean == "= 1" ~ 2,
      Minor_Group_Clean == "> 1" ~ 3,

      # Life cycle duration
      Minor_Group_Clean == "< 1" & Major_Group == "Life cycle duration" ~ 1,
      Minor_Group_Clean == "> 1" & Major_Group == "Life cycle duration" ~ 2,

      # Default for categorical traits (alphabetical or by O/E mean)
      TRUE ~ 99
    ),

    # 2. Apply the ordering
    # We order by Sort_Value first, THEN by the cleaned name
    Trait_Ordered = reorder(Minor_Group_Clean, Sort_Value)
  )


# 2. Create the Faceted Plot
p_validation <- ggplot(plot_data, aes(x = Mean_Bias, y = Trait_Ordered)) +

  # # 1. SHADED ZONE: Shifted to center on 1.0 (representing +/- 10% error)
  # annotate("rect",
  #   xmin = 0.9, xmax = 1.1, ymin = -Inf, ymax = Inf,
  #   fill = "grey93", alpha = 0.8
  # ) +

  # 2. BASELINE: Vertical line at 1.0 (Observed = Expected)
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.6) +

  # Add points colored by Status
  geom_point(aes(color = Status), size = 2) +

  # SCALES: Standardizing your colors
  scale_color_manual(
    values = c(
      "Biased" = "#D55E00",
      "Unbiased" = "#009E73",
      "Zero Variance (Exclude)" = "grey70"
    )
  ) +

  # FACETING
  facet_grid(Major_Group ~ ., scales = "free_y", space = "free_y") +

  # 3. LABELS: Removed "Log10" and updated subtitle
  labs(
    title = "Trait Validation: Model Predictive Performance",
    x = expression(bold("Mean Bias (Observed / Expected Ratio)")),
    y = NULL,
    color = "Model Status"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold", hjust = 0),
    strip.background = element_rect(fill = "grey98"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.text.y = element_text(size = 10)
  )

print(p_validation)

ggsave(
  filename = here("results/strict_trait_validation_rivers.png"),
  width = 180,
  height = 240,
  units = "mm",
  dpi = 600,
  bg = "white",
  plot = p_validation
)


# Save a list of reliable traits into Rdata

reliable_traits <- trait_results %>%
  filter(Reliable == TRUE) %>%
  pull(Trait)

# Print summary of numbers of reliable and unrealiable traits

summary_results <- trait_results %>%
  group_by(Reliable) %>%
  summarise(Count = n())

# Add separation between unbiased and negligible bias for clarity

summary_results <- trait_results %>%
  group_by(Status) %>%
  summarise(Count = n())

save(reliable_traits, file = here("output/strict_reliable_traits_rivers.Rdata"))

save(trait_results, file = here("output/strict_trait_validation_results_rivers.Rdata"))
