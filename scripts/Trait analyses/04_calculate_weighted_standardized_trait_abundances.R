# This script provides a function to standardize trait abundances by normalizing
# each trait category within its trait group.

#' @title Standardize Trait Abundances
#' @description This function takes a trait abundance data frame and standardizes
#' it so that the sum of all trait categories within a given trait group equals 1
#' for each site.
#'
#' @param trait_abundance_df A data frame with sites as rows and trait categories
#' as columns. Column names should be in the format "traitgroup_trait_category".
#'
#' @return A data frame of standardized trait abundances, where each cell
#' represents the proportion of that trait category within its group.
#'
#' @examples
#' # Example trait abundance data
#' # This data would be the output from the previous create_trait_abundance_matrix function.
#' test_trait_abundance <- data.frame(
#'   Maxsize_less_than_0.25 = c(10, 15),
#'   Maxsize_0.25_to_0.5 = c(5, 5),
#'   Maxsize_0.5_to_1 = c(2, 10),
#'   Food_detritus = c(20, 30),
#'   Food_plants = c(10, 10),
#'   row.names = c("site1", "site2")
#' )
#'
#' # Standardize the trait abundances
#' standardized_traits <- standardize_trait_abundances(test_trait_abundance)
#'
#' # Print the standardized data frame
#' print(standardized_traits)
#'
#' # Check the row sums for a trait group (e.g., Maxsize) to verify it sums to 1.
#' cat("\nRow sums for 'Maxsize' trait group:\n")
#' maxsize_cols <- startsWith(names(standardized_traits), "Maxsize")
#' print(rowSums(standardized_traits[, maxsize_cols]))
#'
#' # Expected output for site1: 0.5882353 + 0.2941176 + 0.1176471 = 1
#' # Expected output for site2: 0.4285714 + 0.1428571 + 0.4285714 = 1

standardize_trait_abundances <- function(trait_abundance_df) {
  # Get all column names from the data frame.
  trait_columns <- colnames(trait_abundance_df)

  # Extract the trait group names from the column names.
  # The group name is the part before the first underscore.
  trait_groups <- sapply(strsplit(trait_columns, "_"), `[`, 1)

  # Create a list of column names for each unique trait group.
  # This groups the columns that need to be normalized together.
  grouped_columns <- split(trait_columns, trait_groups)

  # Create an empty data frame to store the standardized results.
  standardized_df <- as.data.frame(matrix(NA, nrow = nrow(trait_abundance_df), ncol = 0))
  rownames(standardized_df) <- rownames(trait_abundance_df)

  # Loop through each trait group.
  for (group_name in names(grouped_columns)) {
    cols_to_normalize <- grouped_columns[[group_name]]

    # Select the columns for the current trait group.
    group_data <- trait_abundance_df[, cols_to_normalize, drop = FALSE]

    # Calculate the sum of each row for the current group.
    # We use rowSums() to get the total abundance for each site within the group.
    # Turn into numeric columns (if not already)

    group_data[] <- lapply(group_data, as.numeric)

    group_sums <- rowSums(group_data)

    # Avoid division by zero by setting sums of 0 to 1 for normalization.
    # This ensures that if a site has no abundance for any category in a group,
    # the proportions will remain 0 rather than becoming NaN.
    group_sums[group_sums == 0] <- 1

    # Normalize each column by dividing by the row sums.
    normalized_data <- group_data / group_sums

    # Add the normalized columns to the final data frame.
    standardized_df <- cbind(standardized_df, normalized_data)
  }

  # Return the final standardized data frame.
  return(standardized_df)
}

# --- Simple Test Case ---

# Create a simplified trait abundance data frame based on the example.
test_trait_abundance <- data.frame(
  trait1_categ1 = c(12.3333, 10.3333),
  trait1_categ2 = c(25.75, 28),
  trait1_categ3 = c(13.9167, 18.6667),
  trait2_categ1 = c(10, 15),
  trait2_categ2 = c(20, 30),
  row.names = c("site1", "site2")
)

# Manually calculate the expected standardized result.
expected_result <- data.frame(
  trait1_categ1 = c(0.23715, 0.18129),
  trait1_categ2 = c(0.49519, 0.49123),
  trait1_categ3 = c(0.26763, 0.32749),
  trait2_categ1 = c(10/30, 15/45),
  trait2_categ2 = c(20/30, 30/45),
  row.names = c("site1", "site2")
)

# Run the function with the test data.
actual_result <- standardize_trait_abundances(test_trait_abundance)

# Compare the actual result with the expected result.
diff_report <- all.equal(actual_result, expected_result)

# Print the results of the test.
cat("Expected Result (from simple matrix multiplication):\n")
print(expected_result)
cat("\nActual Result (from your function):\n")
print(actual_result)
cat("\nTest Result: ", diff_report, "\n")

# Read in the data from the provided CSV file.

library(data.table)
library(here)

trait_expected_abundance_rivers <- fread(here("output/trait_expected_abundance_rivers.csv"))
trait_observed_abundance_rivers <- fread(here("output/trait_observed_abundance_rivers.csv"))

trait_expected_abundance_lakes <- fread(here("output/trait_expected_abundance_lakes.csv"))
trait_observed_abundance_lakes <- fread(here("output/trait_observed_abundance_lakes.csv"))

# Edit the data tables to be data frames

trait_expected_abundance_rivers <- as.data.frame(trait_expected_abundance_rivers)
trait_expected_abundance_lakes <- as.data.frame(trait_expected_abundance_lakes)

trait_observed_abundance_rivers <- as.data.frame(trait_observed_abundance_rivers)
trait_observed_abundance_lakes <- as.data.frame(trait_observed_abundance_lakes)

# Add site_ids

temp_env <- new.env()

load(
  here(
    "scripts",
    "RIVPACS_FOFE",
    "data",
    "Strict_Henna_RIVPACS_data_rivers.RData"
  ),
  envir = temp_env
)

all_sites_rivers <- temp_env$all_sites

trait_expected_abundance_rivers$site_ID <- all_sites_rivers$site_ID
trait_observed_abundance_rivers$site_ID <- all_sites_rivers$site_ID

# Standardize the trait abundances for rivers and lakes.

expected_standardized_traits_rivers <- standardize_trait_abundances(trait_expected_abundance_rivers)
expected_standardized_traits_lakes <- standardize_trait_abundances(trait_expected_abundance_lakes)

observed_standardized_traits_rivers <- standardize_trait_abundances(trait_observed_abundance_rivers)
observed_standardized_traits_lakes <- standardize_trait_abundances(trait_observed_abundance_lakes)

# Check that all trait groups sum to 1 for each site.
check_sums_rivers <- sapply(split(colnames(expected_standardized_traits_rivers),
                                  sub("_.*", "", colnames(expected_standardized_traits_rivers))),
                            function(cols) rowSums(expected_standardized_traits_rivers[, cols, drop = FALSE]))
check_sums_lakes <- sapply(split(colnames(expected_standardized_traits_lakes),
                                 sub("_.*", "", colnames(expected_standardized_traits_lakes))),
                           function(cols) rowSums(expected_standardized_traits_lakes[, cols, drop = FALSE]))


# Save the standardized trait abundances to CSV files.

fwrite(expected_standardized_traits_rivers, here("output/expected_standardized_trait_abundance_rivers.csv"))
fwrite(expected_standardized_traits_lakes, here("output/expected_standardized_trait_abundance_lakes.csv"))

fwrite(observed_standardized_traits_rivers, here("output/observed_standardized_trait_abundance_rivers.csv"))
fwrite(observed_standardized_traits_lakes, here("output/observed_standardized_trait_abundance_lakes.csv"))
