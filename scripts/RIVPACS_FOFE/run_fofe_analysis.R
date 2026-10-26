# ============================================================================
# Analysis settings
# ============================================================================

# 1. Setup
.libPaths()

rm(list = ls())

# ecosystem <- "rivers"

ecosystem <- "lakes"

config <- list(
  ecosystem = ecosystem
)

library(here)

# Repository root
repo_root <- here()

# FoFe workflow root
fofe_dir <- here("scripts", "RIVPACS_FOFE")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd(fofe_dir)

source("R/load_libraries.R")

# 2. Describe the current session and make folders
sesdes <- c("Lake update plots")
sesname <- c("Lake update plots")
source("R/1SetVersionControl.R")

list.files("data")
list.files("R")

# 3. Run scripts
source("R/3ReadModelData.R")
source("R/4RivpacsRunModels.R")
source("R/6SimulateAndCompareFoAndFe.R")

# =======================================================================
# Cross-validation metrics
# =======================================================================

# Only include taxa with Fe > 0.1

FoFeSimu_cv_filtered <- FoFeSimu_cv %>%
  filter(Fe > 0.1)

round(mean(abs(FoFeSimu_cv_filtered$Df)), 3)
round(sd(FoFeSimu_cv_filtered$Df), 3)
round(mean(FoFeSimu_cv_filtered$Fo / FoFeSimu_cv_filtered$Fe), 3)

# Print table of deviations

deviations_table <- data.frame(
  SiteSet = c("Agriculture", "Forestry", "Urban", "Multiple pressures"),
  Increasers = c(deviations[[1]][[1]], deviations[[2]][[1]], deviations[[3]][[1]], deviations[[4]][[1]]),
  Decreasers = c(deviations[[1]][[2]], deviations[[2]][[2]], deviations[[3]][[2]], deviations[[4]][[2]]),
  Total = c(
    deviations[[1]][[1]] + deviations[[1]][[2]],
    deviations[[2]][[1]] + deviations[[2]][[2]],
    deviations[[3]][[1]] + deviations[[3]][[2]],
    deviations[[4]][[1]] + deviations[[4]][[2]]
  )
)

deviations_table


library(patchwork)
library(ggtext)


plots[[1]] <- plots[[1]] +
  ggtitle("Agriculture") +
  theme(axis.title.x = element_blank())

plots[[2]] <- plots[[2]] +
  ggtitle("Forestry") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

plots[[3]] <- plots[[3]] +
  ggtitle("Urban")

plots[[4]] <- plots[[4]] +
  ggtitle("Multiple pressures") +
  theme(axis.title.y = element_blank())

combined_plot <- (plots[[1]] | plots[[2]]) /
  (plots[[3]] | plots[[4]]) +
  plot_layout(guides = "collect") &
  theme(
    plot.title = element_textbox_simple(
      fill = "grey92",
      box.colour = "grey50",
      colour = "black",
      size = 20,
      face = "bold",
      halign = 0.5,
      padding = margin(4, 4, 4, 4),
      margin = margin(0, 0, 5, 0)
    ),
    legend.position = "bottom",
    text = element_text(size = 18),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )


combined_plot

titles <- c(
  "Agriculture",
  "Forestry",
  "Urban",
  "Multiple pressures"
)

for (i in seq_along(zoom_plots)) {
  zoom_plots[[i]] <- zoom_plots[[i]] +
    ggtitle(titles[i])
}

zoom_plots[[1]] <- zoom_plots[[1]] +
  theme(axis.title.x = element_blank())

zoom_plots[[2]] <- zoom_plots[[2]] +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

zoom_plots[[4]] <- zoom_plots[[4]] +
  theme(axis.title.y = element_blank())

zoom_combined <- (zoom_plots[[1]] | zoom_plots[[2]]) /
  (zoom_plots[[3]] | zoom_plots[[4]]) +
  plot_layout(guides = "collect") &
  theme(
    plot.title = element_textbox_simple(
      fill = "grey92",
      box.colour = "grey50",
      colour = "black",
      size = 20,
      face = "bold",
      halign = 0.5,
      padding = margin(4, 4, 4, 4),
      margin = margin(0, 0, 5, 0)
    ),
    legend.position = "bottom",
    text = element_text(size = 18),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )

zoom_combined

# # Calculate mean Fo/Fe ratio with SD and bias
#
# FoFeSimu$ratio <- FoFeSimu$Fo / FoFeSimu$Fe
#
# mean_ratio <- mean(FoFeSimu$ratio, na.rm = TRUE)
# sd_ratio   <- sd(FoFeSimu$ratio, na.rm = TRUE)
#
# bias <- mean(FoFeSimu$Fo - FoFeSimu$Fe, na.rm = TRUE)
#
# print(paste("Mean Fo/Fe ratio:", round(mean_ratio, 3)))
# print(paste("SD Fo/Fe ratio:", round(sd_ratio, 3)))
# print(paste("Bias Fo/Fe ratio:", round(bias, 3)))

if (config$ecosystem == "lakes") {
  ggsave(here("results/FoFe_combined_lakes.png"), combined_plot,
    width = 16, height = 10, dpi = 300
  )
  ggsave(here("results/FoFe_combined_lakes_zoom.png"), zoom_combined,
    width = 16, height = 10, dpi = 300
  )
} else if (config$ecosystem == "rivers") {
  ggsave(here("results/FoFe_combined_rivers.png"), combined_plot,
    width = 16, height = 10, dpi = 300
  )
  ggsave(here("results/FoFe_combined_rivers_zoom.png"), zoom_combined,
    width = 16, height = 10, dpi = 300
  )
}

#
source("R/9WriteSessionDescription.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
