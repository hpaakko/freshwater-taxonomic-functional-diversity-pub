# This R script provides a function to calculate a trait abundance matrix.
# It takes a species abundance matrix and a species trait matrix as input,
# where the species names in the trait data are in a column named "fine_level".

#' @title Create Trait Abundance Matrix
#' @description This function calculates a trait abundance matrix by performing
#' matrix multiplication on a species abundance matrix and a species trait matrix.
#' It correctly aligns the data by using the "fine_level" column in the trait data.
#'
#' @param species_abundance_df A data frame with sites as rows and species as columns.
#' @param species_trait_df A data frame with species as rows (with names in "fine_level" column)
#' and trait categories as other columns.
#'
#' @return A data frame representing the trait abundance matrix.
#' The rows are sites and the columns are trait categories.
#' @examples
#' # Example species abundance data
#' species_abundance <- data.frame(
#'   laj1 = c(13, 0),
#'   laj2 = c(32, 56),
#'   laj3 = c(7, 1)
#' )
#' rownames(species_abundance) <- c("site1", "site2")
#'
#' # Example species trait data with names in a column
#' species_trait <- data.frame(
#'   fine_level = c("laj1", "laj2", "laj3"),
#'   categ1 = c(1, 0, 1),
#'   categ2 = c(3, 1, 0),
#'   categ3 = c(2, 3, 1)
#' )
#'
#' # Calculate the trait abundance matrix
#' trait_abundance_matrix <- create_trait_abundance_matrix(
#'   species_abundance_df = species_abundance,
#'   species_trait_df = species_trait
#' )
#'
#' # Print the resulting matrix
#' print(trait_abundance_matrix)

create_trait_abundance_matrix <- function(species_abundance_df, species_trait_df) {
  # Step 1: Find shared species between abundance data and trait data
  shared_species <- intersect(colnames(species_abundance_df), species_trait_df$fine_level)

  if (length(shared_species) == 0) {
    stop("No common species found between the abundance and trait data.")
  }

  # Step 2: Filter and reorder both data frames to contain only the shared species.
  # This makes the alignment robust.
  species_abundance_aligned <- species_abundance_df[, shared_species, drop = FALSE]
  trait_df_aligned <- species_trait_df[match(shared_species, species_trait_df$fine_level), , drop = FALSE]

  # Step 3: Set the "fine_level" column as the row names for the trait data frame
  # and then remove the column itself.
  rownames(trait_df_aligned) <- trait_df_aligned$fine_level
  trait_df_aligned$fine_level <- NULL

  # Step 4: Identify and select only the numeric columns
  # Force to data.frame to avoid data.table indexing issues
  species_abundance_aligned <- as.data.frame(species_abundance_aligned)
  trait_df_aligned <- as.data.frame(trait_df_aligned)

  is_numeric_abundance <- sapply(species_abundance_aligned, is.numeric)
  species_abundance_aligned <- species_abundance_aligned[, is_numeric_abundance, drop = FALSE]

  is_numeric_trait <- sapply(trait_df_aligned, is.numeric)
  trait_df_aligned <- trait_df_aligned[, is_numeric_trait, drop = FALSE]

  # Step 5: Ensure all numeric columns are of type numeric (double).
  species_abundance_aligned[] <- lapply(species_abundance_aligned, as.numeric)
  trait_df_aligned[] <- lapply(trait_df_aligned, as.numeric)

  # Replace NAs with 0 so they don't break the matrix math
  species_abundance_aligned[is.na(species_abundance_aligned)] <- 0
  trait_df_aligned[is.na(trait_df_aligned)] <- 0

  # Step 6: Convert the data frames to matrices for efficient matrix multiplication.
  abundance_matrix <- as.matrix(species_abundance_aligned)
  trait_matrix <- as.matrix(trait_df_aligned)

  # Step 7: Perform the matrix multiplication using the `%*%` operator.
  trait_abundance_matrix <- abundance_matrix %*% trait_matrix

  # Step 8: Assign meaningful row and column names to the output matrix.
  rownames(trait_abundance_matrix) <- rownames(species_abundance_df)
  colnames(trait_abundance_matrix) <- colnames(trait_df_aligned)

  # Step 9: Return the final result as a data frame.
  return(as.data.frame(trait_abundance_matrix))
}

# --- Simple Test Case using data from the provided screenshot ---

# 1. Input Data
test_species_abundance <- data.frame(
  laj1 = c(13, 0),
  laj2 = c(32, 56),
  laj3 = c(7, 1),
  row.names = c("site1", "site2")
)

test_species_trait <- data.frame(
  fine_level = c("laj1", "laj2", "laj3"),
  categ1 = c(0, 0.16667, 1),
  categ2 = c(0.75, 0.5, 0),
  categ3 = c(0.25, 0.33333, 0)
)

# 2. Manual Calculation Logic (Log10(x+1))
# Site 1 Logs:
# laj1: log10(13+1) = 1.146128
# laj2: log10(32+1) = 1.518514
# laj3: log10(7+1)  = 0.903090

# Site 2 Logs:
# laj1: log10(0+1)  = 0
# laj2: log10(56+1) = 1.755875
# laj3: log10(1+1)  = 0.301030

# Site 1 Matrix Mult:
# categ1: (1.1461 * 0) + (1.5185 * 0.16667) + (0.9031 * 1) = 1.156173
# categ2: (1.1461 * 0.75) + (1.5185 * 0.5) + (0.9031 * 0) = 1.618853
# categ3: (1.1461 * 0.25) + (1.5185 * 0.33333) + (0.9031 * 0) = 0.792701

# Site 2 Matrix Mult:
# categ1: (0 * 0) + (1.7559 * 0.16667) + (0.3010 * 1) = 0.593683
# categ2: (0 * 0.75) + (1.7559 * 0.5) + (0.3010 * 0) = 0.877937
# categ3: (0 * 0.25) + (1.7559 * 0.33333) + (0.3010 * 0) = 0.585292

expected_result <- data.frame(
  categ1 = c(1.156173, 0.593683),
  categ2 = c(1.618853, 0.877937),
  categ3 = c(0.792701, 0.585292),
  row.names = c("site1", "site2")
)

# 3. Execution
# Ensure your 'create_trait_abundance_matrix' function now includes:
# abundance_matrix <- log10(as.matrix(species_abundance_aligned) + 1)

actual_result <- create_trait_abundance_matrix(
  species_abundance_df = test_species_abundance,
  species_trait_df = test_species_trait
)

# 4. Verification
diff_report <- all.equal(actual_result, expected_result, tolerance = 1e-6)

cat("Test Result: ", diff_report, "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Load necessary libraries

library(data.table)
library(here)

# Load data

traits <- fread(here("output/normalized_trait_data.csv"))

# Create a list of unique trait categories by using the first part of colname before "_"

unique_traits <- unique(sapply(strsplit(colnames(traits), "_"), `[`, 1))

# Select columns 1:13, 26:29, 35:39, 57:64, 83:95, 101:103, 109:114

traits <- traits[, c(1:17, 26:29, 35:47, 57:64, 83:95, 101:103, 109:114), with = FALSE]

library(tidyverse)
library(gt)

# Your vector of raw names
raw_traits <- colnames(traits)[-1] # Exclude the first column (fine_level)

clean_trait_table <- data.frame(RawName = raw_traits) %>%
  separate(RawName, into = c("Trait", "Modality"), sep = "_", extra = "merge", remove = FALSE) %>%
  mutate(
    # 1. Clean the Trait Group Names
    Trait = case_when(
      Trait == "Maxsize" ~ "Max Body Size",
      Trait == "Lifecycleduration" ~ "Life cycle duration",
      Trait == "Cycles" ~ "Cycles per year",
      Trait == "Aquaticstages" ~ "Aquatic stage",
      Trait == "Locosubstrate" ~ "Locomotion-substrate relation",
      Trait == "Feedinghabits" ~ "Feeding habit",
      Trait == "Currentvelocity" ~ "Current velocity preference",
      Trait == "Substrate" ~ "Substrate preference",
      Trait == "Temperature" ~ "Temperature preference",
      Trait == "pH" ~ "pH preference",
      TRUE ~ str_to_title(Trait)
    ),

    # 2. Clean the Modalities (The "Shorthand")
    Modality = Modality %>%
      str_replace_all("_", " ") %>%
      str_replace("less than", "<") %>%
      str_replace("more than", ">") %>%
      str_replace("equals", "=") %>%
      str_replace("to", "–") %>% # Uses the professional en-dash
      str_replace("flags boulders cobbles pebbles", "Lithal (boulders–pebbles)") %>%
      str_replace("organic detritus litter", "Organic detritus") %>%
      str_to_sentence()
  )

# View the result
print(clean_trait_table)

# 1. Prepare and collapse the data
trait_summary <- data.frame(RawName = colnames(traits)[-1]
) %>%
  separate(RawName, into = c("Trait", "Modality"), sep = "_", extra = "merge", remove = FALSE) %>%
  dplyr::mutate(
    # 1. Clean the Trait Group Names
    Trait = case_when(
      Trait == "Maxsize" ~ "Max body size",
      Trait == "Lifecycleduration" ~ "Life cycle duration",
      Trait == "Cycles" ~ "Cycles per year",
      Trait == "Aquaticstages" ~ "Aquatic stage",
      Trait == "Locosubstrate" ~ "Locomotion-substrate relation",
      Trait == "Feedinghabits" ~ "Feeding habit",
      Trait == "Currentvelocity" ~ "Current velocity preference",
      Trait == "Substrate" ~ "Substrate preference",
      Trait == "Temperature" ~ "Temperature preference",
      Trait == "pH" ~ "pH preference",
      TRUE ~ str_to_title(Trait)
    ),

    # 2. Clean the Modalities (The "Shorthand")
    Modality = Modality %>%
      str_replace_all("_", " ") %>%
      str_replace("less than", "<") %>%
      str_replace("more than", ">") %>%
      str_replace("equals", "=") %>%
      str_replace_all("\\bto\\b", "–") %>% # Uses the professional en-dash
      str_replace("flags boulders cobbles pebbles", "flags/boulders/cobbles/pebbles") %>%
      str_replace("organic detritus litter", "Organic detritus/litter") %>%
      str_replace("twigs roots", "Twigs/roots") %>%
      str_to_sentence()
  ) %>%
  dplyr::group_by(Trait) %>%
  # Combine all modalities for a single trait into one row
  dplyr::summarise(Modalities = paste(Modality, collapse = ", "), .groups = "drop")

# 2. Create the professional gt table
trait_table <- trait_summary %>%
  gt() %>%
  cols_label(
    Trait = "Functional Trait",
    Modalities = "Modalities (Categories)"
  ) %>%
  # Apply professional academic styling
  tab_options(
    table.font.names = "Times New Roman",
    column_labels.font.weight = "bold",
    table.border.top.color = "black",
    table.border.top.width = px(2),
    table.border.bottom.color = "black",
    table.border.bottom.width = px(2),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(1),
    data_row.padding = px(4)
  ) %>%
  cols_width(
    Trait ~ px(180),
    Modalities ~ px(450)
  )

trait_table

temp_env <- new.env()
temp_env$traits <- traits
save(temp_env, file = file.path(here("output", "traits_for_trait_abundance.RData")))

# Rivers

expected_abundances_river <- fread(here("output/expected_abundances_rivers.csv"))
observed_abundances_river <- fread(here("output/observed_abundances_rivers.csv"))

# Convert to a data frame (ensure site IDs are preserved if they are rownames)
obs_df_rivers <- as.data.frame(observed_abundances_river)

# CRITICAL STEP: Apply log10(x + 1) to match the RIVPACS Exp.LogAbus scale
observed_abundances_river <- log10(obs_df_rivers + 1)

# Lakes

expected_abundances_lake <- fread(here("output/expected_abundances_lakes.csv"))
observed_abundances_lake <- fread(here("output/observed_abundances_lakes.csv"))

# Convert to a data frame (ensure site IDs are preserved if they are rownames)
obs_df_lakes <- as.data.frame(observed_abundances_lake)

# CRITICAL STEP: Apply log10(x + 1) to match the RIVPACS Exp.LogAbus scale
observed_abundances_lake <- log10(obs_df_lakes + 1)

# Remove first column in abundance data frames

expected_abundances_river <- expected_abundances_river[, -1, with = FALSE]
expected_abundances_lake <- expected_abundances_lake[, -1, with = FALSE]

# Edit species names (fine_level column) in species trait data so that space is replaces with .

traits$fine_level <- gsub(" ", ".", traits$fine_level)

# Ensure that the species names in the abundance data match the species names in the trait data
# Print species names which are not found in both abundance and trait data
# Species names are in columns of abundances and fine_level column in traits

species_abundance_names <- colnames(expected_abundances_river)
species_trait_names <- traits$fine_level
missing_in_trait <- setdiff(species_abundance_names, species_trait_names)

# Edit column "Dixiidae" in abundance data to "Dixidae"

expected_abundances_river$Dixidae <- expected_abundances_river$Dixiidae
observed_abundances_river$Dixidae <- observed_abundances_river$Dixiidae

expected_abundances_river$Dixiidae <- NULL
observed_abundances_river$Dixiidae <- NULL

expected_abundances_river$Nepidae <- NULL
observed_abundances_river$Nepidae <- NULL

# Edit the data to be in dataframe format

expected_abundances_lake <- as.data.frame(expected_abundances_lake)
expected_abundances_river <- as.data.frame(expected_abundances_river)
observed_abundances_lake <- as.data.frame(observed_abundances_lake)
observed_abundances_river <- as.data.frame(observed_abundances_river)
traits <- as.data.frame(traits)

# Calculate trait abundance matrices for rivers and lakes

expected_trait_abundances_rivers <- create_trait_abundance_matrix(
  species_abundance_df = expected_abundances_river,
  species_trait_df = traits
)

observed_trait_abundances_rivers <- create_trait_abundance_matrix(
  species_abundance_df = observed_abundances_river,
  species_trait_df = traits
)

expected_trait_abundances_lakes <- create_trait_abundance_matrix(
  species_abundance_df = expected_abundances_lake,
  species_trait_df = traits
)

observed_trait_abundances_lakes <- create_trait_abundance_matrix(
  species_abundance_df = observed_abundances_lake,
  species_trait_df = traits
)

# Save the trait abundance matrices to CSV files

fwrite(expected_trait_abundances_rivers, file.path(here("output", "trait_expected_abundance_rivers.csv")))
fwrite(expected_trait_abundances_lakes, file.path(here("output", "trait_expected_abundance_lakes.csv")))

fwrite(observed_trait_abundances_rivers, file.path(here("output", "trait_observed_abundance_rivers.csv")))
fwrite(observed_trait_abundances_lakes, file.path(here("output", "trait_observed_abundance_lakes.csv")))

