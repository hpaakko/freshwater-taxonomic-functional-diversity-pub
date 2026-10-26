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
