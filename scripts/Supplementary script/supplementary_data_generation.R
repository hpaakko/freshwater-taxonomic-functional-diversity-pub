library(here)
library(data.table)
library(dplyr)

# Load data

# Rivers

load(here("scripts/RIVPACS_FOFE/data/Strict_Henna_RIVPACS_data_rivers.RData"))

river_sites <- all_sites

# Lakes

load(here("scripts/RIVPACS_FOFE/data/Strict_Henna_RIVPACS_data_lakes.RData"))

lake_sites <- all_sites

# Load land use categories

# Rivers

land_use_rivers <- fread(here("output/land_use_categories_rivers.csv"))

# Lakes

land_use_lakes <- fread(here("output/land_use_categories_lakes.csv"))

# Get predictors for rivers and lakes

load(here("output/clustered_data_and_predictors_rivers.RData"))

predictors_rivers <- temp_env$final_predictors
reference_rivers <- temp_env$clustered_data$site_ID

clustered_rivers <- temp_env$clustered_data

load(here("output/clustered_data_and_predictors_lakes.RData"))

predictors_lakes <- temp_env$final_predictors
reference_lakes <- temp_env$clustered_data$site_ID

clustered_lakes <- temp_env$clustered_data

# Add "Reference" value to land use data based on reference site IDs

land_use_rivers[, LandUse_Category := ifelse(site_ID %in% reference_rivers, "Reference", LandUse_Category)]
land_use_lakes[, LandUse_Category := ifelse(site_ID %in% reference_lakes, "Reference", LandUse_Category)]

# Remove sites with "" in LandUse_Category column

land_use_rivers <- land_use_rivers[LandUse_Category != ""]
land_use_lakes <- land_use_lakes[LandUse_Category != ""]


# Get predictor columns from site data data into land use data

river_sites_predictors <- river_sites[, c("site_ID", "catch_area_km2", "AGRICULT", "ARTIFICI", "FORLOSS0021", predictors_rivers), with = FALSE]
lake_sites_predictors <- lake_sites[, c("site_ID", "lake_area_km2", "AGRICULT", "ARTIFICI", "FORLOSS0021", predictors_lakes), with = FALSE]

# Generate supplementary data for rivers and lakes

# Get median and range of MostRecentYear

print(river_sites[, .(median_year = median(MostRecentYear), min_year = min(MostRecentYear), max_year = max(MostRecentYear))])
print(lake_sites[, .(median_year = median(MostRecentYear), min_year = min(MostRecentYear), max_year = max(MostRecentYear))])

# Generate table of median and range of predictors for rivers and lakes by land use category

river_sites_predictors$Log10_catch_area_km2 <- NULL
lake_sites_predictors$Log10_Lake_Area_km2 <- NULL

river_sites_predictors <- merge(river_sites_predictors, land_use_rivers[, .(site_ID, LandUse_Category)], by = "site_ID", all.x = TRUE)
lake_sites_predictors <- merge(lake_sites_predictors, land_use_lakes[, .(site_ID, LandUse_Category)], by = "site_ID", all.x = TRUE)

# Remove rows with NA in LandUse_Category

river_sites_predictors <- river_sites_predictors[!is.na(LandUse_Category)]
lake_sites_predictors <- lake_sites_predictors[!is.na(LandUse_Category)]

table(river_sites_predictors$LandUse_Category)
table(lake_sites_predictors$LandUse_Category)

# Get median and range of AGRICULT, ARTIFICI, FORLOSS0021 per LandUse_Category

landuse_vars <- c("AGRICULT", "ARTIFICI", "FORLOSS0021")

landuse_summary <- river_sites_predictors[
  ,
  lapply(.SD, function(x) {
    sprintf(
      "%.1f (%.1f–%.1f)",
      median(x, na.rm = TRUE),
      min(x, na.rm = TRUE),
      max(x, na.rm = TRUE)
    )
  }),
  by = LandUse_Category,
  .SDcols = landuse_vars
]

landuse_summary

landuse_summary_lakes <- lake_sites_predictors[
  ,
  lapply(.SD, function(x) {
    sprintf(
      "%.1f (%.1f–%.1f)",
      median(x, na.rm = TRUE),
      min(x, na.rm = TRUE),
      max(x, na.rm = TRUE)
    )
  }),
  by = LandUse_Category,
  .SDcols = landuse_vars
]

landuse_summary_lakes

write.table(landuse_summary,
  here("output/landuse_summary_rivers.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(landuse_summary_lakes,
  here("output/landuse_summary_lakes.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Rivers

# Variables to summarize
vars <- setdiff(
  names(river_sites_predictors),
  c("site_ID", "LandUse_Category")
)

summary_table <- river_sites_predictors[
  ,
  lapply(.SD, function(x) {
    sprintf(
      "%.1f (%.1f–%.1f)",
      median(x, na.rm = TRUE),
      min(x, na.rm = TRUE),
      max(x, na.rm = TRUE)
    )
  }),
  by = LandUse_Category,
  .SDcols = vars
]

summary_table

write.table(summary_table,
  here("output/river_predictor_summary.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Lakes

vars <- setdiff(
  names(lake_sites_predictors),
  c("site_ID", "LandUse_Category")
)

summary_table_lakes <- lake_sites_predictors[
  ,
  lapply(.SD, function(x) {
    sprintf(
      "%.1f (%.1f–%.1f)",
      median(x, na.rm = TRUE),
      min(x, na.rm = TRUE),
      max(x, na.rm = TRUE)
    )
  }),
  by = LandUse_Category,
  .SDcols = vars
]

summary_table_lakes

write.table(summary_table_lakes,
  here("output/lake_predictor_summary.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
