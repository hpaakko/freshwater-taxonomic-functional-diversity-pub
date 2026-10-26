# FoFe analysis results reporting

library(data.table)
library(tidyverse)
library(gt)
library(here)
library(taxize)

# Read in the data

# Rivers

agriculture_rivers <- fread(here::here("results/Rivers plots/FoFeSimu_Agriculture more than 10%.csv"))
forestry_rivers <- fread(here::here("results/Rivers plots/FoFeSimu_Forestry more than 20%.csv"))
urban_rivers <- fread(here::here("results/Rivers plots/FoFeSimu_Urban more than 5%.csv"))
multiple_pressures_rivers <- fread(here::here("results/Rivers plots/FoFeSimu_Multiple pressures.csv"))

# Lakes

agriculture_lakes <- fread(here::here("results/Lakes plots new/FoFeSimu_Agriculture more than 10%.csv"))
forestry_lakes <- fread(here::here("results/Lakes plots new/FoFeSimu_Forestry more than 20%.csv"))
urban_lakes <- fread(here::here("results/Lakes plots new/FoFeSimu_Urban more than 5%.csv"))
multiple_pressures_lakes <- fread(here::here("results/Lakes plots new/FoFeSimu_Multiple pressures.csv"))

# Add pressure_type column to each data table

agriculture_rivers[, pressure_type := "Agriculture"]
forestry_rivers[, pressure_type := "Forestry"]
urban_rivers[, pressure_type := "Urban"]
multiple_pressures_rivers[, pressure_type := "Multiple pressures"]

agriculture_lakes[, pressure_type := "Agriculture"]
forestry_lakes[, pressure_type := "Forestry"]
urban_lakes[, pressure_type := "Urban"]
multiple_pressures_lakes[, pressure_type := "Multiple pressures"]

# Create a table of the significant increasers Df > 0

increasers_rivers <- rbindlist(list(
  agriculture_rivers[nulltest == "fail" & Fe > 0.1 & Df > 0],
  forestry_rivers[nulltest == "fail" & Fe > 0.1 & Df > 0],
  urban_rivers[nulltest == "fail" & Fe > 0.1 & Df > 0],
  multiple_pressures_rivers[nulltest == "fail" & Fe > 0.1 & Df > 0]
))

increasers_lakes <- rbindlist(list(
  agriculture_lakes[nulltest == "fail" & Fe > 0.1 & Df > 0],
  forestry_lakes[nulltest == "fail" & Fe > 0.1 & Df > 0],
  urban_lakes[nulltest == "fail" & Fe > 0.1 & Df > 0],
  multiple_pressures_lakes[nulltest == "fail" & Fe > 0.1 & Df > 0]
))

table_increasers_rivers <- increasers_rivers[, .(Row.names, pressure_type, Df, Fe)]
table_increasers_lakes <- increasers_lakes[, .(Row.names, pressure_type, Df, Fe)]

# Edit the Row.names column to be more readable (replace . with space)

table_increasers_rivers[, Row.names := gsub("\\.", " ", Row.names)]
table_increasers_lakes[, Row.names := gsub("\\.", " ", Row.names)]

# Create a table of significant decreasers

decreasers_rivers <- rbindlist(list(
  agriculture_rivers[nulltest == "fail" & Fe > 0.1 & Df < 0],
  forestry_rivers[nulltest == "fail" & Fe > 0.1 & Df < 0],
  urban_rivers[nulltest == "fail" & Fe > 0.1 & Df < 0],
  multiple_pressures_rivers[nulltest == "fail" & Fe > 0.1 & Df < 0]
))

decreasers_lakes <- rbindlist(list(
  agriculture_lakes[nulltest == "fail" & Fe > 0.1 & Df < 0],
  forestry_lakes[nulltest == "fail" & Fe > 0.1 & Df < 0],
  urban_lakes[nulltest == "fail" & Fe > 0.1 & Df < 0],
  multiple_pressures_lakes[nulltest == "fail" & Fe > 0.1 & Df < 0]
))

table_decreasers_rivers <- decreasers_rivers[, .(Row.names, pressure_type, Df, Fe)]
table_decreasers_lakes <- decreasers_lakes[, .(Row.names, pressure_type, Df, Fe)]

table_decreasers_rivers[, Row.names := gsub("\\.", " ", Row.names)]
table_decreasers_lakes[, Row.names := gsub("\\.", " ", Row.names)]

increasers_rivers <- table_increasers_rivers
increasers_lakes <- table_increasers_lakes
decreasers_rivers <- table_decreasers_rivers
decreasers_lakes <- table_decreasers_lakes

# Round Df and Fe to 3 decimal places

table_decreasers_lakes[, Df := round(Df, 3)]
table_decreasers_lakes[, Fe := round(Fe, 3)]

table_decreasers_rivers[, Df := round(Df, 3)]
table_decreasers_rivers[, Fe := round(Fe, 3)]

table_increasers_lakes[, Df := round(Df, 3)]
table_increasers_lakes[, Fe := round(Fe, 3)]

table_increasers_rivers[, Df := round(Df, 3)]
table_increasers_rivers[, Fe := round(Fe, 3)]

# Save as csv

fwrite(table_increasers_rivers, here::here("results/Summary of taxa responses - Rivers increasers.csv"))
fwrite(table_increasers_lakes, here::here("results/Summary of taxa responses - Lakes increasers.csv"))
fwrite(table_decreasers_rivers, here::here("results/Summary of taxa responses - Rivers decreasers.csv"))
fwrite(table_decreasers_lakes, here::here("results/Summary of taxa responses - Lakes decreasers.csv"))

all_taxa_long <- bind_rows(
  increasers_rivers %>% mutate(Status = "Inc", Ecosystem = "River"),
  decreasers_rivers %>% mutate(Status = "Dec", Ecosystem = "River"),
  increasers_lakes %>% mutate(Status = "Inc", Ecosystem = "Lake"),
  decreasers_lakes %>% mutate(Status = "Dec", Ecosystem = "Lake")
)

all_taxa_long <- all_taxa_long %>%
  dplyr::mutate(Row.names = dplyr::case_when(
    Row.names == "Elmis aena" ~ "Elmis aenea",
    Row.names == "Baetis subalpinus_vernus_fuscatus" ~ "Baetis subalpinus",
    Row.names == "Baetis niger_digitatus" ~ "Baetis niger",
    Row.names == "Radix_Lymnaea" ~ "Lymnaeidae",
    TRUE ~ as.character(Row.names) # Ensures types match exactly
  ))

add_taxonomy <- function(df, taxon_col = "Row.names", db = "ncbi") {
  # 1. Clean the dots automatically
  df[[taxon_col]] <- gsub("\\.", " ", df[[taxon_col]])
  df[[taxon_col]] <- gsub("_", " ", df[[taxon_col]])

  # 2. Extract unique taxa
  unique_taxa <- unique(df[[taxon_col]])
  message(paste("Querying", length(unique_taxa), "unique taxa from", db, "..."))

  # 3. Get classification hierarchy
  hier <- classification(unique_taxa, db = db, rows = 1)

  not_found <- names(hier)[sapply(hier, function(x) all(is.na(x)))]

  # Search from itis if not found in ncbi

  if (length(not_found) > 0) {
    message(paste("Querying", length(not_found), "taxa from ITIS as fallback..."))
    hier_itis <- classification(not_found, db = "itis", rows = 1)

    # Update the original hierarchy list with ITIS results
    for (taxon in not_found) {
      if (!is.na(hier_itis[[taxon]])) {
        hier[[taxon]] <- hier_itis[[taxon]]
      }
    }
  }

  # 4. Convert the list to a table using bind_rows
  # We filter for only the successful matches first
  hierarchy_table <- hier[!is.na(hier)] %>%
    dplyr::bind_rows(.id = "query") %>%
    as_tibble() %>%
    dplyr::filter(rank %in% c("class", "order", "family", "genus", "subclass", "class")) %>%
    # Pivot so each rank is a column
    pivot_wider(id_cols = query, names_from = rank, values_from = name) %>%
    rename(!!taxon_col := query)

  # 5. Join back to your original data
  final_df <- df %>%
    left_join(hierarchy_table, by = taxon_col)

  return(final_df)
}

# --- 2. Run your taxonomy function ---
# This will clean the names and add columns for Class, Order, Family, etc.
message("Fetching taxonomy data...")
all_taxa_with_info <- add_taxonomy(all_taxa_long, taxon_col = "Row.names", db = "ncbi")

# --- 3. Pivot into the Summary Matrix ---
# Note: We include the new taxonomic columns (order, family, class, etc.)
# in the selection so they stay as columns on the left side of the matrix.

test_df <- all_taxa_with_info %>%
  mutate(
    Pressure_Group = paste(Ecosystem, pressure_type, sep = "_"),
    Response_Value = paste0(Status, " (", round(Df, 3), ")")
  ) %>%
  dplyr::select(
    any_of(c("class", "order", "family", "Row.names")),
    Pressure_Group, Response_Value
  )

sapply(test_df, class)

summary_matrix <- all_taxa_with_info %>%
  # 1. Create the unique identifier for columns
  mutate(Pressure_Group = paste(Ecosystem, pressure_type, sep = "_")) %>%
  # 2. Create a combined label: "Inc (0.45)" or "Dec (-0.32)"
  # We round Df to 2 or 3 decimal places for cleanliness
  mutate(Response_Value = paste0(Status, " (", round(Df, 3), ")")) %>%
  # 3. Select Taxonomy, Taxon Name, the Pressure Group, and our new Response Value
  dplyr::select(any_of(c("class", "order", "family", "Row.names")), Pressure_Group, Response_Value) %>%
  distinct() %>%
  # 4. Pivot wider
  pivot_wider(
    names_from = Pressure_Group,
    values_from = Response_Value,
    values_fn = first,
    values_fill = list(Response_Value = "-")
  ) %>%
  # 5. Sort by Taxonomy
  arrange(across(any_of(c("order", "family", "Row.names"))))

# --- 4. Save the Summary Matrix ---

write_csv(summary_matrix, here::here("results/Summary of taxa responses.csv"))
