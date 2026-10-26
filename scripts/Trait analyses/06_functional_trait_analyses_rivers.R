library(tidyverse)
library(vegan) # For Multivariate Analysis (RDA, PCA)
library(rstatix) # For easier statistical testing (ANOVA, post-hoc)
library(data.table)
library(Rmisc)
library(lme4)
library(lmerTest)
library(ggplot2)
library(patchwork)
library(here)

temp_env <- new.env()

load(
  here::here(
    "scripts",
    "RIVPACS_FOFE",
    "data",
    "Strict_Henna_RIVPACS_data_rivers.RData"
  ),
  envir = temp_env
)

all_sites_rivers <- temp_env$all_sites


table(all_sites_rivers$country)

ref_sites <- fread(here::here("data/ref_rivers_with_clusters.csv"))

ref_site_ids <- ref_sites$site_ID

# --- 2. Data Loading and Preparation ------------------------------------------

# Observed Trait Abundance Matrix (O_T)
# Values are observed abundances (O_Sp).
O_T_matrix <- fread(here::here("output/observed_standardized_trait_abundance_rivers.csv"))

# Expected Trait Abundance Matrix (E_T)
# Values are the RIVPACS model's predicted Expected Abundances/Probabilities (E_Sp).
E_T_matrix <- fread(here::here("output/expected_standardized_trait_abundance_rivers.csv"))

# Assign the Site ID column (as per your instructions)
O_T_matrix$site_ID <- all_sites_rivers$site_ID
E_T_matrix$site_ID <- all_sites_rivers$site_ID

# Convert to data.frame and set site_ID as row names for matrix operations
O_T_matrix <- O_T_matrix %>%
  column_to_rownames(var = "site_ID")

E_T_matrix <- E_T_matrix %>%
  column_to_rownames(var = "site_ID")

save(O_T_matrix, file = here::here("output/observed_trait_matrix_rivers.Rdata"))
save(E_T_matrix, file = here::here("output/expected_trait_matrix_rivers.Rdata"))

# Environmental/Land-Use Data (Site metadata)

# Must contain site_ID and the categorical land-use group (e.g., "Reference", "Agriculture > 20%")
Land_Use_Data <- fread(here::here("output/land_use_categories_rivers.csv")) %>%
  # Ensure site_ID is a character column
  mutate(site_ID = as.character(site_ID))

GROUPING_COLUMN_NAME <- "LandUse_Category"

# Add new LandUse_Category Reference to sites that are in ref_site_ids

Land_Use_Data <- Land_Use_Data %>%
  mutate(LandUse_Category = ifelse(site_ID %in% ref_site_ids, "Reference", LandUse_Category))

# Create O/E matrix

trait_epsilons <- apply(O_T_matrix, 2, function(x) {
  min_val <- min(x[x > 0], na.rm = TRUE)
  # Use half of the smallest observed value as the epsilon
  return(min_val / 2)
})

OE_T_matrix <- O_T_matrix # Start with raw
for (trait in colnames(O_T_matrix)) {
  eps <- trait_epsilons[trait]
  OE_T_matrix[, trait] <- ifelse(E_T_matrix[, trait] < 0.01, NA, (O_T_matrix[, trait] + eps) / (E_T_matrix[, trait] + eps))
}

save(OE_T_matrix, file = here::here("output/OE_trait_matrix_rivers_raw.Rdata"))

# OE_T_matrix <- log10(OE_T_matrix)
#
# save(OE_T_matrix, file = here::here("output/OE_trait_matrix_rivers.Rdata")

# Create the grouping data frame, ensuring site_ID is the row key.
grouping_df <- Land_Use_Data %>%
  filter(site_ID %in% rownames(OE_T_matrix)) %>%
  data.frame() %>%
  column_to_rownames(var = "site_ID")

grouping_factor_final <- as.factor(grouping_df[[GROUPING_COLUMN_NAME]])

# 1. Create a "Shell" dataframe from your matrix Site IDs
matrix_sites <- data.frame(Site = rownames(OE_T_matrix))

metadata <- matrix_sites %>%
  left_join(Land_Use_Data %>% select(site_ID, Land_Use = LandUse_Category), by = c("Site" = "site_ID"))

# Combine metadata Land_Use column based on Site columnd to o_T_matrix based on rowname

O_T_matrix_landuse <- as.data.frame(O_T_matrix) %>%
  dplyr::mutate(Site = rownames(O_T_matrix)) %>%
  left_join(metadata, by = "Site")

E_T_matrix_landuse <- as.data.frame(E_T_matrix) %>%
  mutate(Site = rownames(E_T_matrix)) %>%
  left_join(metadata, by = "Site")

# Put Site and Land_Use columns first in order

O_T_matrix_landuse <- O_T_matrix_landuse %>%
  select(Site, Land_Use, everything())

E_T_matrix_landuse <- E_T_matrix_landuse %>%
  select(Site, Land_Use, everything())

fwrite(O_T_matrix_landuse, here::here("output/river_trait_O_with_landuse.csv"), row.names = FALSE)

fwrite(E_T_matrix_landuse, here::here("output/river_trait_E_with_landuse.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------'
# T-TEST FOR EACH TRAIT TO ASSESS SIGNIFICANT DEVIATIONS FROM EXPECTED (O - E) PER LAND USE Land_Use, mu=0
# ------------------------------------------------------------------------------

Land_Use_groups <- c("Forestry", "Agriculture", "Urban", "Multiple pressures")


metadata <- metadata %>%
  mutate(Land_Use = as.character(Land_Use)) %>%
  filter(Land_Use %in% Land_Use_groups)

# Remove Land_Use NAs from metadata and corresponding rows from OE_T_matrix

metadata <- metadata %>%
  filter(!is.na(Land_Use))

OE_T_matrix <- OE_T_matrix[rownames(OE_T_matrix) %in% metadata$Site, ]

# Import reliable trait list from the trait reliability analysis

load(here::here("output/strict_reliable_traits_rivers.Rdata"))

load(here::here("output/strict_trait_validation_results_rivers.Rdata"))

test_results <- data.frame(
  Trait = character(),
  Land_Use = character(),
  Mean_OE = numeric(),
  P_Value_Adj = numeric(),
  stringsAsFactors = FALSE
)

# 1. Prepare the Long Format Data
long_ratios <- as.data.frame(OE_T_matrix[, reliable_traits]) %>%
  mutate(Site = rownames(.)) %>%
  left_join(metadata %>% select(Site, Land_Use), by = "Site") %>%
  pivot_longer(
    cols = all_of(reliable_traits),
    names_to = "Trait",
    values_to = "Ratio"
  )

long_ratios_all <- as.data.frame(OE_T_matrix) %>%
  mutate(Site = rownames(.)) %>%
  left_join(metadata %>% select(Site, Land_Use), by = "Site") %>%
  pivot_longer(
    cols = -c(Site, Land_Use),
    names_to = "Trait",
    values_to = "Ratio"
  )

test_results <- long_ratios %>%
  group_by(Trait, Land_Use) %>%
  dplyr::summarise(
    n_total = n(),
    # Count how many non-NA Ratios we actually have
    n_valid = sum(!is.na(Ratio)),
    Mean_OE_Ratio = mean(Ratio, na.rm = TRUE),
    SD_Ratio = sd(Ratio, na.rm = TRUE),

    # NEW LOGIC:
    # 1. Must have at least 3 valid observations
    # 2. Standard deviation must be greater than 0
    P_Value = if (n_valid > 2 && !is.na(SD_Ratio) && SD_Ratio > 0) {
      t.test(Ratio, mu = 1)$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(P_Value_Adj = p.adjust(P_Value, method = "BH"))

# 3. Apply the "Ecological Significance" Threshold (0.1 Log Units)
test_results <- test_results %>%
  mutate(Status = case_when(
    # Increaser: Statistically significant AND > 25% increase (0.1 log units)
    P_Value_Adj < 0.05 & Mean_OE_Ratio > 1 ~ "Increaser",

    # Decreaser: Statistically significant AND > 20% decrease (-0.1 log units)
    P_Value_Adj < 0.05 & Mean_OE_Ratio < 1 ~ "Decreaser",
    TRUE ~ "Not Significant"
  ))

test_results_all <- long_ratios_all %>%
  group_by(Trait, Land_Use) %>%
  dplyr::summarise(
    n_total = n(),
    # Count how many non-NA Ratios we actually have
    n_valid = sum(!is.na(Ratio)),
    Mean_OE_Ratio = mean(Ratio, na.rm = TRUE),
    SD_Ratio = sd(Ratio, na.rm = TRUE),

    # NEW LOGIC:
    # 1. Must have at least 3 valid observations
    # 2. Standard deviation must be greater than 0
    P_Value = if (n_valid > 2 && !is.na(SD_Ratio) && SD_Ratio > 0) {
      t.test(Ratio, mu = 1)$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(P_Value_Adj = p.adjust(P_Value, method = "BH"))

test_results_all <- test_results_all %>%
  mutate(Status = case_when(
    # Increaser: Statistically significant AND > 25% increase (0.1 log units)
    P_Value_Adj < 0.05 & Mean_OE_Ratio > 0 ~ "Increaser",

    # Decreaser: Statistically significant AND > 20% decrease (-0.1 log units)
    P_Value_Adj < 0.05 & Mean_OE_Ratio < 0 ~ "Decreaser",

    # Significant but small: Statistically different but below our 0.1 threshold
    P_Value_Adj < 0.05 ~ "Negligible Shift",
    TRUE ~ "Not Significant"
  ))

# --- 1. DATA PREPARATION ---
# Ensure we have all necessary columns for the plot
plot_data <- test_results %>%
  # Join validation info
  left_join(trait_results %>% select(Trait, Validation_Status = Status), by = "Trait") %>%
  # Split the raw Trait name (e.g., Maxsize_less_than_0.25)
  separate(Trait,
    into = c("Major_Group", "Minor_Group"), sep = "_",
    extra = "merge", remove = FALSE
  ) %>%
  mutate(
    # A. Clean the Major Group (Facet Labels)
    Major_Group = case_when(
      Major_Group == "Maxsize" ~ "Max body size",
      Major_Group == "Lifecycleduration" ~ "LCD",
      Major_Group == "Cycles" ~ "CPY",
      Major_Group == "Aquaticstages" ~ "Aquatic stage",
      Major_Group == "Locosubstrate" ~ "Locomotion",
      Major_Group == "Feedinghabits" ~ "Feeding habit",
      Major_Group == "Currentvelocity" ~ "Current",
      Major_Group == "Substrate" ~ "Substrate",
      Major_Group == "Temperature" ~ "Temperature",
      Major_Group == "pH" ~ "pH",
      TRUE ~ str_to_title(Major_Group)
    ),

    # B. Clean the Minor Group (X-axis Labels)
    Minor_Group_Clean = Minor_Group %>%
      str_replace_all("_", " ") %>%
      str_replace("less than", "<") %>%
      str_replace("more than", ">") %>%
      str_replace("equals", "=") %>%
      str_replace_all("\\bto\\b", "–") %>%
      str_replace("flags boulders cobbles pebbles", "Coarse") %>%
      str_replace("organic detritus litter", "Organic detritus/litter") %>%
      str_replace("twigs roots", "Twigs/roots") %>%
      str_to_sentence(),

    # C. Order Land Use
    Land_Use = factor(Land_Use, levels = c(
      "Agriculture", "Forestry", "Urban", "Multiple pressures"
    )),

    # D. Stats and Status
    se = SD_Ratio / sqrt(as.numeric(n_valid)),
    Lower_CI = Mean_OE_Ratio - (1.96 * se),
    Upper_CI = Mean_OE_Ratio + (1.96 * se), ,
    Impact_Status = case_when(
      P_Value_Adj < 0.05 & Mean_OE_Ratio > 1 ~ "Increaser",
      P_Value_Adj < 0.05 & Mean_OE_Ratio < 1 ~ "Decreaser",
      TRUE ~ "Not Significant"
    )
  ) %>%
  # Ordering for the X-axis within facets
  mutate(Trait_Ordered = reorder(Minor_Group_Clean, Mean_OE_Ratio))

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
      Minor_Group_Clean == "< 1" & Major_Group == "LCD" ~ 1,
      Minor_Group_Clean == "> 1" & Major_Group == "LCD" ~ 2,

      # Default for categorical traits (alphabetical or by O/E mean)
      TRUE ~ 99
    ),

    # 2. Apply the ordering
    # We order by Sort_Value first, THEN by the cleaned name
    Trait_Ordered = reorder(Minor_Group_Clean, Sort_Value)
  )

# Remove rows with NA in Mean_OE_Ratio to avoid plotting issues

plot_data <- plot_data %>%
  filter(!is.na(Mean_OE_Ratio))


# Select only the traits that had a significant shift in at least one land use for the final plot

plot_data <- plot_data %>%
  group_by(Trait) %>%
  filter(any(Impact_Status %in% c("Increaser", "Decreaser"))) %>%
  ungroup()

# Print plot data for cleaned results

plot_data %>%
  filter(Impact_Status %in% c("Increaser", "Decreaser")) %>%
  select(Trait, Major_Group, Minor_Group_Clean, Land_Use, Mean_OE_Ratio, Lower_CI, Upper_CI, Impact_Status) %>%
  print(n = Inf)

plot_data <- plot_data %>%
  mutate(
    Impact_Status = factor(
      Impact_Status,
      levels = c("Increaser", "Decreaser", "Not Significant")
    )
  )

color_winner <- "#0072B2" # Deep Blue
color_loser <- "#D55E00" # Vermillion
color_ns <- "grey80" # Neutral Grey

fwrite(plot_data, here::here("output/traits_rivers_plot_data.csv"), row.names = FALSE)

# Facet grid

final_plot <- ggplot(plot_data, aes(x = Trait_Ordered, y = Mean_OE_Ratio)) +
  # Reference Zone (Shaded)
  # annotate("rect",
  #   xmin = -Inf, xmax = Inf, ymin = 0.8, ymax = 1.2,
  #   fill = "grey95", alpha = 0.8
  # ) +

  # # Reference Lines
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.3) +
  # geom_hline(yintercept = c(0.8, 1.2), linetype = "dotted", colour = "grey70", linewidth = 0.2) +

  # Data Layers
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI),
    colour = "grey30", width = 0.3, linewidth = 0.4
  ) +
  geom_point(aes(fill = Impact_Status, shape = Impact_Status),
    size = 2.5, stroke = 0.3, colour = "black"
  ) +

  # # Outlier Label (Positioned slightly lower to stay within 3.0 limit)
  # geom_label(
  #   data = subset(plot_data, Mean_OE_Ratio > 3),
  #   aes(y = 2.8, label = round(Mean_OE_Ratio, 1)),
  #   vjust = 0, color = color_winner, size = 2.8, fontface = "bold",
  #   label.padding = unit(0.1, "lines"), fill = "white", alpha = 0.8
  # ) +

  # Faceting & Scaling
  # Note: space = "free_x" is removed to keep facet widths equal so titles aren't cut
  facet_grid(Land_Use ~ Major_Group,
    scales = "free_x",
    space = "free_x",
    labeller = labeller(Major_Group = label_wrap_gen(width = 10))
  ) +

  # Scales
  scale_fill_manual(values = c(
    "Increaser" = color_winner,
    "Decreaser" = color_loser,
    "Not Significant" = color_ns
  )) +
  scale_shape_manual(values = c(
    "Increaser" = 21,
    "Decreaser" = 21,
    "Not Significant" = 21
  )) +
  scale_y_continuous(breaks = seq(0, 2, by = 0.5), expand = expansion(mult = c(0, 0.05))) +

  # Theme Customization
  theme_classic(base_size = 11) +
  theme(
    # Panel & Facet borders
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "grey92", colour = "grey50", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 9, lineheight = 0.9, angle = 90),
    strip.placement = "outside",
    strip.text.x = element_text(size = 9, margin = ggplot2::margin(t = 5, b = 5)),
    strip.text.y = element_text(size = 9, margin = ggplot2::margin(r = 5, l = 5), angle = 90),

    # Axis Alignment (90 degree fix)
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = "black", size = 10),
    axis.text.y = element_text(colour = "black", size = 10),
    axis.title.y = element_text(face = "bold", margin = ggplot2::margin(r = 10), size = 10),

    # Legend - Safe at the bottom
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),

    # Spacing - panel.spacing.y provides room for rows
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.5, "lines"),
    plot.margin = ggplot2::margin(5, 5, 5, 5)
  ) +
  labs(x = NULL, y = expression(bold("Observed / Expected Ratio"))) +
  coord_cartesian(ylim = c(0, 2))

print(final_plot)

ggsave(here::here("results/strict_traits_rivers.tiff"),
  plot = final_plot,
  width = 180,
  height = 240,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

# Plot all traits for supplement

plot_data_all <- test_results_all %>%
  # Join validation info
  left_join(trait_results %>% dplyr::select(Trait, Validation_Status = Status), by = "Trait") %>%
  # Split the raw Trait name (e.g., Maxsize_less_than_0.25)
  separate(Trait,
    into = c("Major_Group", "Minor_Group"), sep = "_",
    extra = "merge", remove = FALSE
  ) %>%
  mutate(
    # A. Clean the Major Group (Facet Labels)
    Major_Group = case_when(
      Major_Group == "Maxsize" ~ "Max body size",
      Major_Group == "Lifecycleduration" ~ "LCD",
      Major_Group == "Cycles" ~ "CPY",
      Major_Group == "Aquaticstages" ~ "Aquatic stage",
      Major_Group == "Locosubstrate" ~ "Locomotion",
      Major_Group == "Feedinghabits" ~ "Feeding habit",
      Major_Group == "Currentvelocity" ~ "Current",
      Major_Group == "Substrate" ~ "Substrate",
      Major_Group == "Temperature" ~ "Temperature",
      Major_Group == "pH" ~ "pH",
      TRUE ~ str_to_title(Major_Group)
    ),

    # B. Clean the Minor Group (X-axis Labels)
    Minor_Group_Clean = Minor_Group %>%
      str_replace_all("_", " ") %>%
      str_replace("less than", "<") %>%
      str_replace("more than", ">") %>%
      str_replace("equals", "=") %>%
      str_replace_all("\\bto\\b", "–") %>%
      str_replace("flags boulders cobbles pebbles", "Coarse") %>%
      str_replace("organic detritus litter", "Organic detritus/litter") %>%
      str_replace("twigs roots", "Twigs/roots") %>%
      str_to_sentence(),

    # C. Order Land Use
    Land_Use = factor(Land_Use, levels = c(
      "Agriculture", "Forestry", "Urban", "Multiple pressures"
    )),

    # D. Stats and Status
    se = SD_Ratio / sqrt(as.numeric(n_valid)),
    Lower_CI = Mean_OE_Ratio - (1.96 * se),
    Upper_CI = Mean_OE_Ratio + (1.96 * se), ,
    Impact_Status = case_when(
      P_Value_Adj < 0.05 & Mean_OE_Ratio > 1 ~ "Increaser",
      P_Value_Adj < 0.05 & Mean_OE_Ratio < 1 ~ "Decreaser",
      TRUE ~ "Not Significant"
    )
  ) %>%
  # Ordering for the X-axis within facets
  mutate(Trait_Ordered = reorder(Minor_Group_Clean, Mean_OE_Ratio))

plot_data_all <- plot_data_all %>%
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
      Minor_Group_Clean == "< 1" & Major_Group == "LCD" ~ 1,
      Minor_Group_Clean == "> 1" & Major_Group == "LCD" ~ 2,

      # Default for categorical traits (alphabetical or by O/E mean)
      TRUE ~ 99
    ),

    # 2. Apply the ordering
    # We order by Sort_Value first, THEN by the cleaned name
    Trait_Ordered = reorder(Minor_Group_Clean, Sort_Value)
  )

# Remove rows with NA in Mean_OE_Ratio to avoid plotting issues

plot_data_all <- plot_data_all %>%
  filter(!is.na(Mean_OE_Ratio))

plot_data_all <- plot_data_all %>%
  mutate(
    Impact_Status = factor(
      Impact_Status,
      levels = c("Increaser", "Decreaser", "Not Significant")
    )
  )

final_plot_all <- ggplot(plot_data_all, aes(x = Trait_Ordered, y = Mean_OE_Ratio)) +
  # Reference Zone (Shaded)
  # annotate("rect",
  #   xmin = -Inf, xmax = Inf, ymin = 0.8, ymax = 1.2,
  #   fill = "grey95", alpha = 0.8
  # ) +

  # Reference Lines
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.3) +
  # geom_hline(yintercept = c(0.8, 1.2), linetype = "dotted", colour = "grey70", linewidth = 0.2) +

  # Data Layers
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI),
    colour = "grey30", width = 0.3, linewidth = 0.4
  ) +
  geom_point(aes(fill = Impact_Status, shape = Impact_Status),
    size = 2.5, stroke = 0.3, colour = "black"
  ) +

  # # Outlier Label (Positioned slightly lower to stay within 3.0 limit)
  # geom_label(
  #   data = subset(plot_data, Mean_OE_Ratio > 3),
  #   aes(y = 2.8, label = round(Mean_OE_Ratio, 1)),
  #   vjust = 0, color = color_winner, size = 2.8, fontface = "bold",
  #   label.padding = unit(0.1, "lines"), fill = "white", alpha = 0.8
  # ) +

  # Faceting & Scaling
  # Note: space = "free_x" is removed to keep facet widths equal so titles aren't cut
  facet_grid(Land_Use ~ Major_Group,
    scales = "free_x",
    space = "free_x",
    labeller = labeller(Major_Group = label_wrap_gen(width = 10))
  ) +

  # Scales
  scale_fill_manual(values = c(
    "Increaser" = color_winner,
    "Decreaser" = color_loser,
    "Not Significant" = color_ns
  )) +
  scale_shape_manual(values = c(
    "Increaser" = 21,
    "Decreaser" = 21,
    "Not Significant" = 21
  )) +
  scale_y_continuous(breaks = seq(0, 2, by = 0.5), expand = expansion(mult = c(0, 0.05))) +

  # Theme Customization
  theme_classic(base_size = 11) +
  theme(
    # Panel & Facet borders
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "grey92", colour = "grey50", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 9, lineheight = 0.9, angle = 90),
    strip.placement = "outside",
    strip.text.x = element_text(size = 9, margin = ggplot2::margin(t = 5, b = 5)),
    strip.text.y = element_text(size = 9, margin = ggplot2::margin(r = 5, l = 5), angle = 90),

    # Axis Alignment (90 degree fix)
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = "black", size = 8),
    axis.text.y = element_text(colour = "black", size = 10),
    axis.title.y = element_text(face = "bold", margin = ggplot2::margin(r = 10), size = 10),

    # Legend - Safe at the bottom
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),

    # Spacing - panel.spacing.y provides room for rows
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.5, "lines"),
    plot.margin = ggplot2::margin(5, 5, 5, 5)
  ) +
  labs(x = NULL, y = expression(bold("Observed / Expected Ratio"))) +
  coord_cartesian(ylim = c(0, 2))

print(final_plot_all)

ggsave(here::here("results/traits_rivers_all.tiff"),
  plot = final_plot_all,
  width = 180,
  height = 240,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

# --- 1. CREATE THE CATEGORIZED SUMMARY ---
results_summary <- test_results %>%
  # Bring in validation to keep things organized
  dplyr::left_join(trait_results %>% dplyr::select(Trait, Validation_Status = Status), by = "Trait") %>%
  mutate(
    # Define the 5 key ecological categories
    Category = case_when(
      P_Value_Adj < 0.05 & Mean_OE_Ratio > 0.1 ~ "Major Increaser",
      P_Value_Adj < 0.05 & Mean_OE_Ratio < -0.1 ~ "Major Decreaser",
      P_Value_Adj < 0.05 & abs(Mean_OE_Ratio) <= 0.1 ~ "Minor Significant Shift",
      TRUE ~ "Not Significant"
    )
  )

# --- 2. CREATE A FREQUENCY TALLY (For your Results Section) ---
# This shows exactly how many traits fall into each bucket per Land Use
category_tally <- results_summary %>%
  group_by(Land_Use, Category) %>%
  tally() %>%
  pivot_wider(names_from = Land_Use, values_from = n, values_fill = 0)

# --- 3. CREATE A LIST OF TRAITS PER CATEGORY (For your Appendix) ---
# This makes it easy to see "Which traits are the Major Losers in Agriculture?"
trait_lookup_table <- results_summary %>%
  select(Land_Use, Category, Trait, Mean_OE_Ratio, P_Value_Adj) %>%
  arrange(Land_Use, Category, desc(abs(Mean_OE_Ratio))) %>%
  filter(Category %in% c("Major Increaser", "Major Decreaser", "Minor Significant Shift"))

# Print the tally to your console
print("Summary Tally of Trait Responses:")
print(category_tally)
