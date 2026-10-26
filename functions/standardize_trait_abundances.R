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
