library(here)

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

# Generate ids for different land use criteria

# Load multiple pressure ids for rivers

# temp_mp <- new.env()
#
# load(
#   here("data", "multiple_pressure_site_ids.RData"),
#   envir = temp_mp
# )
#
# multiple_pressure_ids <- temp_mp$multiple_pressure_sites
#
# multiple_pressures_ids_rivers <- which(
#   all_sites_rivers$site_ID %in% multiple_pressure_ids
# )

agriculture_ids_rivers <- which(all_sites_rivers$AGRICULT > 10 & all_sites_rivers$FORLOSS0021 < 15 & all_sites_rivers$ARTIFICI < 2)
forestry_ids_rivers <- which(all_sites_rivers$FORLOSS0021 > 20 & all_sites_rivers$AGRICULT < 2 & all_sites_rivers$ARTIFICI < 1)
urban_ids_rivers <- which(all_sites_rivers$ARTIFICI > 5 & all_sites_rivers$AGRICULT < 2 & all_sites_rivers$FORLOSS0021 < 15)

multiple_pressures_ids_rivers <- which(
  (all_sites_rivers$AGRICULT > 10 & all_sites_rivers$FORLOSS0021 > 20) |
    (all_sites_rivers$AGRICULT > 10 & all_sites_rivers$ARTIFICI > 5) |
    (all_sites_rivers$FORLOSS0021 > 20 & all_sites_rivers$ARTIFICI > 5)
)

# Initialize a new column for land use categories
all_sites_rivers$LandUse_Category <- NA_character_

# Assign land use categories based on your predefined criteria
all_sites_rivers$LandUse_Category[agriculture_ids_rivers] <- "Agriculture"
all_sites_rivers$LandUse_Category[forestry_ids_rivers] <- "Forestry"
all_sites_rivers$LandUse_Category[urban_ids_rivers] <- "Urban"
all_sites_rivers$LandUse_Category[multiple_pressures_ids_rivers] <- "Multiple pressures"

land_use_data <- all_sites_rivers[, c("site_ID", "LandUse_Category")]

# Save the land use data to a CSV file

fwrite(land_use_data, here("output/land_use_categories_rivers.csv"))
