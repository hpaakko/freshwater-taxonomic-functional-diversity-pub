# Site specific O/E ratio comparison taxa vs trait

library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(vegan)
library(here)

# Import O/E ratios

# Streams

load(here::here("output/bio_data_all_sites_with_OE_rivers.Rdata"))

taxa_data_rivers <- bio_data_all_sites_2 %>%
  dplyr::select(site_ID, OoverE)

load(here::here("output/OE_trait_matrix_rivers_raw.Rdata"))

trait_data_rivers_raw <- OE_T_matrix

trait_data_rivers_raw <- trait_data_rivers_raw %>%
  rownames_to_column("site_ID")

# Combine taxa and trait data for rivers

combined_rivers <- taxa_data_rivers %>%
  left_join(trait_data_rivers_raw, by = "site_ID")

Land_Use_Data <- fread(here::here("output/land_use_categories_rivers.csv")) %>%
  # Ensure site_ID is a character column
  mutate(site_ID = as.character(site_ID))

combined_rivers <- combined_rivers %>%
  left_join(Land_Use_Data, by = "site_ID")

# Remove sites with LandUse_Category empty string

combined_rivers <- combined_rivers %>%
  filter(LandUse_Category != "")

# Put site_ID and LandUse Category at the front

combined_rivers <- combined_rivers %>%
  dplyr::select(site_ID, LandUse_Category, everything())

# Rename OverE column to Taxa_OE for clarity

colnames(combined_rivers)[which(colnames(combined_rivers) == "OoverE")] <- "Taxa_OE"

# Lakes

load(here::here("output/bio_data_all_sites_with_OE_lakes.Rdata"))

taxa_data_lakes <- bio_data_all_sites_2 %>%
  dplyr::select(site_ID, OoverE)

load(here::here("output/OE_trait_matrix_lakes_raw.Rdata"))

trait_data_lakes_raw <- OE_T_matrix

trait_data_lakes_raw <- trait_data_lakes_raw %>%
  rownames_to_column("site_ID")

combined_lakes <- taxa_data_lakes %>%
  left_join(trait_data_lakes_raw, by = "site_ID")

Land_Use_Data <- fread(here::here("output/land_use_categories_lakes.csv")) %>%
  # Ensure site_ID is a character column
  mutate(site_ID = as.character(site_ID))

combined_lakes <- combined_lakes %>%
  left_join(Land_Use_Data, by = "site_ID")

# Remove sites with LandUse_Category empty string

combined_lakes <- combined_lakes %>%
  filter(LandUse_Category != "")

# Put site_ID and LandUse Category at the front

combined_lakes <- combined_lakes %>%
  dplyr::select(site_ID, LandUse_Category, everything())

# Rename OverE column to Taxa_OE for clarity

colnames(combined_lakes)[which(colnames(combined_lakes) == "OoverE")] <- "Taxa_OE"

fwrite(combined_rivers, here::here("output/combined_rivers_taxa_trait_OE.csv"))

fwrite(combined_lakes, here::here("output/combined_lakes_taxa_trait_OE.csv"))

# Correlation function

analyze_oe_correlation <- function(trait_data, taxa_data,
                                   method = "spearman",
                                   group_name = "All Sites") {
  library(dplyr)

  # --- 1. Remove site_ID and keep matrix ---
  trait_matrix <- trait_data %>%
    dplyr::select(-site_ID) %>%
    dplyr::select(-LandUse_Category)

  # --- 2. Calculate absolute deviation ---
  trait_dev <- abs(trait_matrix - 1)

  # --- 3. Mean deviation per site ---
  functional_dev <- trait_dev %>%
    rowMeans(na.rm = TRUE)

  # --- 4. Combine with site_ID ---
  functional_dev_df <- data.frame(
    site_ID = trait_data$site_ID,
    functional_dev = functional_dev
  )

  # --- 5. Merge with taxonomic data ---
  combined <- functional_dev_df %>%
    left_join(taxa_data, by = "site_ID") %>%
    filter(!is.na(OoverE))

  # --- 6. Optional: taxonomic deviation
  combined <- combined %>%
    mutate(taxa_dev = 1 - OoverE)

  # --- 7. Correlation test ---
  cor_result <- cor.test(
    combined$OoverE,
    combined$functional_dev,
    method = method,
    use = "complete.obs"
  )

  # --- 8. Plot ---


  format_p <- function(p) {
    if (is.na(p)) {
      return("p = NA")
    } else if (p < 0.001) {
      return("p < 0.001")
    } else {
      return(paste0("p = ", signif(p, 2)))
    }
  }


  # Combine
  plot_data <- combined %>%
    mutate(system = group_name)

  base_size <- 20


  p <- ggplot(plot_data, aes(x = OoverE, y = functional_dev)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    facet_wrap(~system, nrow = 1) +
    labs(
      x = "Taxonomic condition (O/E taxa)",
      y = "Functional deviation (mean |O/E − 1|)"
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1.5)) +
    theme_bw(base_size = base_size) +
    theme(
      strip.text = element_text(size = base_size, face = "bold"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size * 0.9),
      panel.grid = element_blank()
    )


  rho <- round(cor_result$estimate, 2)

  annot_df <- data.frame(
    x = c(1.1),
    y = c(0.8)
  )

  p <- p + geom_text(
    data = annot_df,
    aes(x = x, y = y),
    label = paste0(
      "rho = ", round(as.numeric(rho), 2),
      "\n", format_p(cor_result$p.value)
    ),
    size = base_size * 0.3,
  )


  print(p)
  print(cor_result$p.value)


  # --- 9. Output ---
  return(list(
    data = combined,
    correlation = cor_result,
    plot = p
  ))
}


# Run analyses for both rivers and lakes

# RIVERS

# Divide data into land use groups

Land_Use_Data <- fread(here::here("output/land_use_categories_rivers.csv")) %>%
  # Ensure site_ID is a character column
  mutate(site_ID = as.character(site_ID))

ref_sites_rivers <- fread(here::here("data/ref_rivers_with_clusters.csv")) %>%
  # Ensure site_ID is a character column
  mutate(site_ID = as.character(site_ID))

ref_sites_rivers_ids <- ref_sites_rivers$site_ID

Land_Use_Data <- Land_Use_Data %>%
  mutate(LandUse_Category = ifelse(site_ID %in% ref_sites_rivers_ids, "Reference", LandUse_Category))

# Merge land use data with taxa and trait summaries

taxa_data_rivers <- taxa_data_rivers %>%
  left_join(Land_Use_Data, by = "site_ID")

trait_data_rivers_raw <- trait_data_rivers_raw %>%
  left_join(Land_Use_Data, by = "site_ID")

# Remove sites with LandUse_Category NA

taxa_data_rivers <- taxa_data_rivers %>%
  filter(LandUse_Category != "")

trait_data_rivers <- trait_data_rivers_raw %>%
  filter(LandUse_Category != "")

all_sites_results_rivers <- analyze_oe_correlation(
  taxa_data = taxa_data_rivers,
  trait_data = trait_data_rivers,
  group_name = "All Streams"
)

print(all_sites_results_rivers$correlation)

# Analysis for Agricultural sites
agri_results_rivers <- analyze_oe_correlation(
  taxa_data = taxa_data_rivers %>% filter(LandUse_Category == "Agriculture"),
  trait_data = trait_data_rivers %>% filter(LandUse_Category == "Agriculture"),
  group_name = "Agricultural Streams"
)

# View the plot
print(agri_results_rivers$plot)

# View the correlation coefficients
print(agri_results_rivers$correlation)

# Analysis for Forestry sites

for_results_rivers <- analyze_oe_correlation(
  taxa_data = taxa_data_rivers %>% filter(LandUse_Category == "Forestry"),
  trait_data = trait_data_rivers %>% filter(LandUse_Category == "Forestry"),
  group_name = "Forestry Streams"
)

print(for_results_rivers$plot)

# Analysis for Urban sites

urb_results_rivers <- analyze_oe_correlation(
  taxa_data = taxa_data_rivers %>% filter(LandUse_Category == "Urban"),
  trait_data = trait_data_rivers %>% filter(LandUse_Category == "Urban"),
  group_name = "Urban Streams"
)

print(urb_results_rivers$plot)

# Analysis for Multiple pressures sites

multi_results_rivers <- analyze_oe_correlation(
  taxa_data = taxa_data_rivers %>% filter(LandUse_Category == "Multiple pressures"),
  trait_data = trait_data_rivers %>% filter(LandUse_Category == "Multiple pressures"),
  group_name = "Multiple Pressure Streams"
)

print(multi_results_rivers$plot)

# Combine correlation results into one file

combined_correlation_results_rivers <- tibble(
  Group = c(
    "All Streams",
    "Agriculture",
    "Forestry",
    "Urban",
    "Multiple pressures"
  ),
  Rho = c(
    as.numeric(all_sites_results_rivers$correlation$estimate),
    as.numeric(agri_results_rivers$correlation$estimate),
    as.numeric(for_results_rivers$correlation$estimate),
    as.numeric(urb_results_rivers$correlation$estimate),
    as.numeric(multi_results_rivers$correlation$estimate)
  ),
  P_Value = c(
    all_sites_results_rivers$correlation$p.value,
    agri_results_rivers$correlation$p.value,
    for_results_rivers$correlation$p.value,
    urb_results_rivers$correlation$p.value,
    multi_results_rivers$correlation$p.value
  )
)

fwrite(combined_correlation_results_rivers, here::here("output/Correlation_results_rivers.csv"))

# Combine all plots into one figure for rivers

agri_results_rivers$plot <- agri_results_rivers$plot + theme(axis.title.x = element_blank())
for_results_rivers$plot <- for_results_rivers$plot + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank()
)
multi_results_rivers$plot <- multi_results_rivers$plot + theme(axis.title.y = element_blank())

combined_plot_rivers <- cowplot::plot_grid(
  agri_results_rivers$plot,
  for_results_rivers$plot,
  urb_results_rivers$plot,
  multi_results_rivers$plot,
  ncol = 2
)

print(combined_plot_rivers)

# Save plot

ggsave(
  filename = here::here("results/Fig_taxa_function_coupling_rivers.png"),
  plot = combined_plot_rivers,
  width = 17,
  height = 12,
  dpi = 300
)

# LAKES

# Divide data into land use groups

Land_Use_Data <- fread(here::here("output/land_use_categories_lakes.csv")) %>%
  # Ensure site_ID is a character column
  mutate(site_ID = as.character(site_ID))

ref_sites_lakes <- fread(here::here("data/lbioREF_N207_clustering.csv")) %>%
  # Ensure site_ID is a character column
  mutate(site_ID = as.character(site_ID))

ref_sites_lakes_ids <- ref_sites_lakes$site_ID

Land_Use_Data <- Land_Use_Data %>%
  mutate(LandUse_Category = ifelse(site_ID %in% ref_sites_lakes_ids, "Reference", LandUse_Category))

# Merge land use data with taxa and trait summaries

taxa_data_lakes <- taxa_data_lakes %>%
  left_join(Land_Use_Data, by = "site_ID")

trait_data_lakes_raw <- trait_data_lakes_raw %>%
  left_join(Land_Use_Data, by = "site_ID")

# Remove sites with LandUse_Category NA

taxa_data_lakes <- taxa_data_lakes %>%
  filter(LandUse_Category != "")

trait_data_lakes <- trait_data_lakes_raw %>%
  filter(LandUse_Category != "")

all_sites_results_lakes <- analyze_oe_correlation(
  taxa_data = taxa_data_lakes,
  trait_data = trait_data_lakes,
  group_name = "All Lakes"
)

print(all_sites_results_lakes$correlation)
print(all_sites_results_lakes$correlation$p.value)

# Analysis for Agricultural sites

agri_results_lakes <- analyze_oe_correlation(
  taxa_data = taxa_data_lakes %>% filter(LandUse_Category == "Agriculture"),
  trait_data = trait_data_lakes %>% filter(LandUse_Category == "Agriculture"),
  group_name = "Agricultural Lakes"
)

print(agri_results_lakes$plot)


# Analysis for Forestry sites

for_results_lakes <- analyze_oe_correlation(
  taxa_data = taxa_data_lakes %>% filter(LandUse_Category == "Forestry"),
  trait_data = trait_data_lakes %>% filter(LandUse_Category == "Forestry"),
  group_name = "Forestry Lakes"
)

print(for_results_lakes$plot)

# Analysis for Urban sites

urb_results_lakes <- analyze_oe_correlation(
  taxa_data = taxa_data_lakes %>% filter(LandUse_Category == "Urban"),
  trait_data = trait_data_lakes %>% filter(LandUse_Category == "Urban"),
  group_name = "Urban Lakes"
)

print(urb_results_lakes$plot)

# Analysis for Multiple pressures sites

multi_results_lakes <- analyze_oe_correlation(
  taxa_data = taxa_data_lakes %>% filter(LandUse_Category == "Multiple pressures"),
  trait_data = trait_data_lakes %>% filter(LandUse_Category == "Multiple pressures"),
  group_name = "Multiple Pressure Lakes"
)

print(multi_results_lakes$plot)

# Combine correlation results into one file

combined_correlation_results_lakes <- tibble(
  Group = c(
    "All Lakes",
    "Agriculture",
    "Forestry",
    "Urban",
    "Multiple pressures"
  ),
  Rho = c(
    as.numeric(all_sites_results_lakes$correlation$estimate),
    as.numeric(agri_results_lakes$correlation$estimate),
    as.numeric(for_results_lakes$correlation$estimate),
    as.numeric(urb_results_lakes$correlation$estimate),
    as.numeric(multi_results_lakes$correlation$estimate)
  ),
  P_Value = c(
    all_sites_results_lakes$correlation$p.value,
    agri_results_lakes$correlation$p.value,
    for_results_lakes$correlation$p.value,
    urb_results_lakes$correlation$p.value,
    multi_results_lakes$correlation$p.value
  )
)

fwrite(combined_correlation_results_lakes, here::here("output/Correlation_results_lakes.csv"))

agri_results_lakes$plot <- agri_results_lakes$plot + theme(axis.title.x = element_blank())
for_results_lakes$plot <- for_results_lakes$plot + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
multi_results_lakes$plot <- multi_results_lakes$plot + theme(axis.title.y = element_blank())

# Combine all plots into one figure for lakes

combined_plot_lakes <- cowplot::plot_grid(
  agri_results_lakes$plot,
  for_results_lakes$plot,
  urb_results_lakes$plot,
  multi_results_lakes$plot,
  ncol = 2
)

print(combined_plot_lakes)

# Save plot

ggsave(
  filename = here::here("results/Fig_taxa_function_coupling_lakes.png"),
  plot = combined_plot_lakes,
  width = 17,
  height = 12,
  dpi = 300
)

all_river_data <- bind_rows(
  agri_results_rivers$data %>% mutate(Group = "Agriculture"),
  for_results_rivers$data %>% mutate(Group = "Forestry"),
  urb_results_rivers$data %>% mutate(Group = "Urban"),
  multi_results_rivers$data %>% mutate(Group = "Multiple")
)

streams_plot <- all_sites_results_rivers$data %>%
  mutate(system = "Streams")

lakes_plot <- all_sites_results_lakes$data %>%
  mutate(system = "Lakes")

# Combine
plot_data <- bind_rows(streams_plot, lakes_plot) %>%
  mutate(
    system = factor(system, levels = c("Streams", "Lakes"))
  )


p <- ggplot(plot_data, aes(x = OoverE, y = functional_dev)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  facet_wrap(
    ~ factor(system, levels = c("Streams", "Lakes")),
    nrow = 1
  ) +
  labs(
    x = "Taxonomic condition (O/E taxa)",
    y = "Functional deviation (mean |O/E − 1|)"
  ) +
  scale_y_continuous(limits = c(0, 1.3)) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    panel.grid = element_blank()
  )


rho_streams <- round(all_sites_results_rivers$correlation$estimate, 2)
rho_lakes <- round(all_sites_results_lakes$correlation$estimate, 2)

p_value_streams <- all_sites_results_rivers$correlation$p.value
p_value_lakes <- all_sites_results_lakes$correlation$p.value

annot_df <- data.frame(
  system = c("Streams", "Lakes"),
  x = c(1.1, 1.1),
  y = c(
    max(1.1, na.rm = TRUE),
    max(1.1, na.rm = TRUE)
  ),
  label = c(
    paste0("rho == ", round(as.numeric(rho_streams), 2)),
    paste0("rho == ", round(as.numeric(rho_lakes), 2))
  )
)

annot_df$system <- factor(
  annot_df$system,
  levels = c("Streams", "Lakes")
)

p <- p +
  geom_text(
    data = annot_df,
    aes(x = x, y = y, label = label),
    parse = TRUE
  )


print(p)

# Save as png

ggsave(
  filename = here::here("results/Fig_taxa_function_coupling.png"),
  plot = p,
  width = 170, # mm
  height = 85, # mm (nice 2:1 ratio)
  units = "mm",
  dpi = 300
)
