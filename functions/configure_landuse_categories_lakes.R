library(here)

temp_env <- new.env()

load(
  here(
    "scripts",
    "RIVPACS_FOFE",
    "data",
    "Strict_Henna_RIVPACS_data_lakes.RData"
  ),
  envir = temp_env
)

all_sites_lakes <- temp_env$all_sites

# Generate ids for different land use criteria

agriculture_ids_lakes <- which(all_sites_lakes$AGRICULT > 10 & all_sites_lakes$FORLOSS0021 < 15 & all_sites_lakes$ARTIFICI < 2)
forestry_ids_lakes <- which(all_sites_lakes$FORLOSS0021 > 20 & all_sites_lakes$AGRICULT < 2 & all_sites_lakes$ARTIFICI < 1)
urban_ids_lakes <- which(all_sites_lakes$ARTIFICI > 5 & all_sites_lakes$AGRICULT < 2 & all_sites_lakes$FORLOSS0021 < 15)

# This code finds sites with (Ag & For) OR (Ag & Art) OR (For & Art).
# A site with all three pressures will satisfy all three of these OR conditions.
multiple_pressures_ids_lakes <- which(
  (all_sites_lakes$AGRICULT > 10 & all_sites_lakes$FORLOSS0021 > 20) |
    (all_sites_lakes$AGRICULT > 10 & all_sites_lakes$ARTIFICI > 5) |
    (all_sites_lakes$FORLOSS0021 > 20 & all_sites_lakes$ARTIFICI > 5)
)

# Initialize a new column for land use categories

all_sites_lakes$LandUse_Category <- NA_character_

# Assign land use categories based on your predefined criteria
all_sites_lakes$LandUse_Category[agriculture_ids_lakes] <- "Agriculture"
all_sites_lakes$LandUse_Category[forestry_ids_lakes] <- "Forestry"
all_sites_lakes$LandUse_Category[urban_ids_lakes] <- "Urban"
all_sites_lakes$LandUse_Category[multiple_pressures_ids_lakes] <- "Multiple pressures"

land_use_data <- all_sites_lakes[, c("site_ID", "LandUse_Category")]

# Save the land use data to a CSV file

fwrite(land_use_data, here("output/land_use_categories_lakes.csv"))
