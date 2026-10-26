# ============================================================================
# 02_normalize_trait_data.R
#
# Purpose:
#  - Normalize trait data for use in functional trait analyses.
#
# Input:
#  - Excel file containing trait data (combined_traits_Tachet_Hennalle.xlsx).
#
# Output:
# - CSV file containing the normalized trait data (normalized_trait_data.csv).
#
# Manuscript:
# Differential land use effects on taxonomic and functional
# freshwater biodiversity in Nordic streams and lakes
# ===========================================================================

# Script to create a trait data frame and normalize it

library(readxl)
library(tidyr)
library(dplyr)
library(here)

# Step 1: Import the trait data

trait_data <- read_excel(here("data/combined_traits_Tachet_Hennalle.xlsx"), sheet = "clean_data")

# ---- FIX: Add this line to remove incomplete rows ----
trait_tibble <- trait_data %>% drop_na(fine_level)

# Step 2: Reshape the data to a long format.
long_format_data <- trait_tibble %>%
  pivot_longer(
    cols = -fine_level,
    names_to = "trait_category",
    values_to = "affinity_score"
  )

# Step 3: Extract the trait names from the category names.
normalized_data <- long_format_data %>%
  separate(trait_category, into = c("trait_name", "category"), sep = "_", extra = "merge") %>%
  distinct()
# Step 4: Group by `fine_level` and `trait_name`, then normalize.
normalized_data <- normalized_data %>%
  group_by(fine_level, trait_name) %>%
  mutate(normalized_score = if (sum(affinity_score, na.rm = TRUE) == 0) {
    0
  } else {
    affinity_score / sum(affinity_score, na.rm = TRUE)
  }) %>%
  ungroup()

# Step 5: Reshape the data back to its original "wide" format.
final_normalized_tibble <- normalized_data %>%
  select(-affinity_score) %>%
  pivot_wider(
    names_from = c(trait_name, category),
    values_from = normalized_score,
    names_sep = "_"
  )

# Step 6: Convert all columns (except the first) to numeric using a for loop.
for (col in names(final_normalized_tibble)[-1]) {
  final_normalized_tibble[[col]] <- as.numeric(unlist(final_normalized_tibble[[col]]))
}

print("Final Normalized Tibble with numeric columns:")
print(final_normalized_tibble)

# Step 7: Verification
print("\nVerifying row sums for 'Maxsize' trait:")
final_normalized_tibble %>%
  select(-fine_level) %>%
  select(starts_with("Maxsize")) %>%
  rowSums(na.rm = TRUE)

# You can also verify the class of the columns
print("\nClass of the final columns:")
sapply(final_normalized_tibble, class)

# Step 8: Save the final normalized tibble to a CSV file.

write.csv(final_normalized_tibble, here("output/normalized_trait_data.csv"), row.names = FALSE)
