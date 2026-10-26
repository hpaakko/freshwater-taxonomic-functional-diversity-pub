# Make sure sdirfigs is defined in your global environment or passed to the function
# Example (replace with your actual path):
# sdirfigs <- "C:/Your/Path/To/Figures/"
# SiteSetName <- "YourSiteName" # Example

# Ensure the directory exists (good practice!)
# if (!dir.exists(sdirfigs)) {
#   dir.create(sdirfigs, recursive = TRUE)
# }


PlotFoFe_single <- function(FoFeSimu, SiteSetName, captureprobs) {
  color_winner <- "#0072B2"
  color_loser <- "#D55E00"
  color_stable <- "grey50"

  taxon_colors <- ifelse(FoFeSimu$Fo > FoFeSimu$p975, color_winner,
    ifelse(FoFeSimu$Fo < FoFeSimu$p025, color_loser,
      color_stable
    )
  )

  df <- FoFeSimu

  df$trend <- ifelse(df$Fo > df$p975, "Increaser",
    ifelse(df$Fo < df$p025, "Decreaser", "Not Significant")
  )

  df$trend <- factor(df$trend, levels = c("Increaser", "Decreaser", "Not Significant"))


  base_size <- 20

  p <- ggplot(df, aes(x = Fe, y = Fo)) +
    geom_point(aes(color = trend)) +
    scale_color_manual(
      values = c(
        "Increaser" = "#0072B2",
        "Decreaser" = "#D55E00",
        "Not Significant" = "grey50"
      ),
      breaks = c("Increaser", "Decreaser", "Not Significant"),
      name = "Taxon response"
    ) +
    geom_line(aes(y = p025), linetype = "dashed", colour = "grey30") +
    geom_line(aes(y = p975), linetype = "dashed", colour = "grey30") +
    # Add red line at Fe = 0.1
    geom_vline(xintercept = 0.1, linetype = "dashed", colour = "red") +
    geom_abline(slope = 1, intercept = 0, colour = "grey80", linetype = "dotted") +
    labs(
      x = "Expected frequency (Fe)",
      y = "Observed frequency (Fo)"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(size = base_size * 0.9),
      axis.title = element_text(size = base_size),
    )

  NrOf02taxa <- dim(FoFeSimu[FoFeSimu$Fe > 0.1, ])[1] # Corrected subsetting
  NrOf02decr <- dim(FoFeSimu[FoFeSimu$nulltest == "fail" & FoFeSimu$Fe > 0.1 & FoFeSimu$Df < 0, ])[1] # Corrected subsetting
  NrOf02incr <- dim(FoFeSimu[FoFeSimu$nulltest == "fail" & FoFeSimu$Fe > 0.1 & FoFeSimu$Df > 0, ])[1] # Corrected subsetting
  Perc02decr <- round(NrOf02decr / NrOf02taxa * 100, 0)
  Perc02incr <- round(NrOf02incr / NrOf02taxa * 100, 0)

  # Return plot and deviation percentages

  return(list(p, Perc02incr, Perc02decr))
}

PlotFoFe_zoom <- function(FoFeSimu, ymax) {
  base_size <- 20

  df <- FoFeSimu

  df$trend <- ifelse(df$Fo > df$p975, "Increaser",
    ifelse(df$Fo < df$p025, "Decreaser", "Not Significant")
  )

  df$trend <- factor(
    df$trend,
    levels = c("Increaser", "Decreaser", "Not Significant")
  )

  p <- ggplot(df, aes(x = Fe, y = Fo)) +
    geom_point(aes(color = trend)) +
    scale_color_manual(
      values = c(
        "Increaser" = "#0072B2",
        "Decreaser" = "#D55E00",
        "Not Significant" = "grey50"
      ),
      breaks = c("Increaser", "Decreaser", "Not Significant"),
      name = "Taxon response"
    ) +
    geom_line(aes(y = p025),
      linetype = "dashed",
      colour = "grey30"
    ) +
    geom_line(aes(y = p975),
      linetype = "dashed",
      colour = "grey30"
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      colour = "grey80",
      linetype = "dotted"
    ) +
    labs(
      x = "Expected frequency (Fe)",
      y = "Observed frequency (Fo)"
    ) +
    coord_cartesian(
      xlim = c(0, 0.1),
      ylim = c(0, ymax * 1.05)
    ) +
    theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(size = base_size * 0.9),
      axis.title = element_text(size = base_size)
    )

  return(p)
}


# Function to plot taxa Fo and Fe results and save as pdf, emf, and png
PlotFoFe <- function(FoFeSimu, SiteSetName, captureprobs) { # Added sdirfigs as argument
  # or ensure it's globally available

  print("")
  print("*4. Plot FoFe function")
  print("")

  # Original SiteSetName (the source of the error):
  # SiteSetName <- "Agriculture more than 10%"

  # Step 1: Replace "%" with "perc" or "percent"
  SafeSiteSetName <- gsub("%", "perc", SiteSetName, fixed = TRUE)

  # Step 2: Replace spaces with underscores or hyphens
  SafeSiteSetName <- gsub(" ", "_", SafeSiteSetName, fixed = TRUE)

  # --- Plot 1: Fo_Fe without taxon labels ---

  # 1. Save as PDF
  pdf(paste(sdirfigs, "Fo_Fe_", SafeSiteSetName, ".pdf", sep = ""), onefile = TRUE)
  plot(FoFeSimu$Fe, FoFeSimu$Fo,
    main = paste(SiteSetName, " (N=", dim(capture_probs)[1], ")", sep = ""),
    ylab = "Observed frequence", xlab = "Expected frequence", ylim = c(0, 1), xlim = c(0, 1),
    las = 1, box(bty = "l"), pch = 19, cex = 1.3, cex.lab = 1.3, cex.axis = 1.3,
    frame.plot = FALSE, col = 1
  )
  lines(FoFeSimu$Fe, FoFeSimu$p025, lty = 2)
  lines(FoFeSimu$Fe, FoFeSimu$p975, lty = 2)
  NrOf02taxa <- dim(FoFeSimu[FoFeSimu$Fe > 0.1, ])[1] # Corrected subsetting
  NrOf02decr <- dim(FoFeSimu[FoFeSimu$nulltest == "fail" & FoFeSimu$Fe > 0.1 & FoFeSimu$Df < 0, ])[1] # Corrected subsetting
  NrOf02incr <- dim(FoFeSimu[FoFeSimu$nulltest == "fail" & FoFeSimu$Fe > 0.1 & FoFeSimu$Df > 0, ])[1] # Corrected subsetting
  Perc02decr <- round(NrOf02decr / NrOf02taxa * 100, 0)
  Perc02incr <- round(NrOf02incr / NrOf02taxa * 100, 0)
  dev.off() # Close the PDF device

  # 2. Save as EMF (Enhanced Metafile - good for Word/PowerPoint on Windows)
  # Use win.metafile() for EMF on Windows
  win.metafile(paste(sdirfigs, "Fo_Fe_", SafeSiteSetName, ".emf", sep = ""))
  # Re-plot the exact same plot content for EMF
  plot(FoFeSimu$Fe, FoFeSimu$Fo,
    main = paste(SiteSetName, " (N=", dim(capture_probs)[1], ")", sep = ""),
    ylab = "Observed frequence", xlab = "Expected frequence", ylim = c(0, 1), xlim = c(0, 1),
    las = 1, box(bty = "l"), pch = 19, cex = 1.3, cex.lab = 1.3, cex.axis = 1.3,
    frame.plot = FALSE, col = 1
  )
  lines(FoFeSimu$Fe, FoFeSimu$p025, lty = 2)
  lines(FoFeSimu$Fe, FoFeSimu$p975, lty = 2)
  dev.off() # Close the EMF device

  # 1. Define the colors to match your trait plot
  color_winner <- "#0072B2"
  color_loser <- "#D55E00"
  color_stable <- "grey50"

  # 2. Create a color vector based on the O/E relationship
  # Logic: If Observed (Fo) > Upper Bound (p975) = Red
  #        If Observed (Fo) < Lower Bound (p025) = Blue
  taxon_colors <- ifelse(FoFeSimu$Fo > FoFeSimu$p975, color_winner,
    ifelse(FoFeSimu$Fo < FoFeSimu$p025, color_loser,
      color_stable
    )
  )

  # 3. Save as PNG
  png(paste(sdirfigs, "Fo_Fe_", SafeSiteSetName, ".png", sep = ""), width = 800, height = 600, res = 100)

  plot(FoFeSimu$Fe, FoFeSimu$Fo,
    main = paste(SiteSetName, " (N=", dim(capture_probs)[1], ")", sep = ""),
    ylab = "Observed frequency",
    xlab = "Expected frequency",
    ylim = c(0, 1), xlim = c(0, 1),
    las = 1,
    pch = 19,
    cex = 1.3,
    cex.lab = 1.3,
    cex.axis = 1.3,
    frame.plot = FALSE,
    col = taxon_colors, # Applied the color vector here
    axes = FALSE # Handled manually for 'l' box style
  )

  # Custom axes for the 'l' box look
  axis(1)
  axis(2, las = 1)

  # Add the 1:1 reference line (Optional but helpful for O/E plots)
  abline(0, 1, col = "grey80", lty = 3)

  # Add the 95% Confidence Interval lines
  lines(FoFeSimu$Fe, FoFeSimu$p025, lty = 2, col = "grey30")
  lines(FoFeSimu$Fe, FoFeSimu$p975, lty = 2, col = "grey30")

  dev.off()

  # --- Plot 2: Fo_Fe with short taxon names ---

  # make short taxon names using vegan and add them to plots
  # (assuming vegan package is installed and loaded)
  # library(vegan) # Uncomment if not already loaded
  shortnames <- make.cepnames(FoFeSimu$Row.names)

  # 1. Save as PDF
  pdf(paste(sdirfigs, "Fo_Fe_taxa_", SafeSiteSetName, ".pdf", sep = ""), onefile = TRUE)
  plot(FoFeSimu$Fe, FoFeSimu$Fo,
    main = paste(SiteSetName, " (N=", dim(capture_probs)[1], ")", sep = ""),
    ylab = "Observed frequence", xlab = "Expected frequence", ylim = c(0, 1), xlim = c(0, 1),
    las = 1, box(bty = "l"), pch = 19, cex = 1.3, cex.lab = 1.3, cex.axis = 1.3,
    frame.plot = FALSE, col = rgb(0, 0.34509803921568627, 0.32941176470588235)
  )
  lines(FoFeSimu$Fe, FoFeSimu$p025, lty = 2)
  lines(FoFeSimu$Fe, FoFeSimu$p975, lty = 2)
  text(FoFeSimu$Fe, FoFeSimu$Fo, labels = shortnames, cex = .5, pos = 3, col = "darkgray")
  dev.off() # Close the PDF device

  # 2. Save as EMF
  win.metafile(paste(sdirfigs, "Fo_Fe_taxa_", SafeSiteSetName, ".emf", sep = ""))
  # Re-plot for EMF
  plot(FoFeSimu$Fe, FoFeSimu$Fo,
    main = paste(SiteSetName, " (N=", dim(capture_probs)[1], ")", sep = ""),
    ylab = "Observed frequence", xlab = "Expected frequence", ylim = c(0, 1), xlim = c(0, 1),
    las = 1, box(bty = "l"), pch = 19, cex = 1.3, cex.lab = 1.3, cex.axis = 1.3,
    frame.plot = FALSE, col = rgb(0, 0.34509803921568627, 0.32941176470588235)
  )
  lines(FoFeSimu$Fe, FoFeSimu$p025, lty = 2)
  lines(FoFeSimu$Fe, FoFeSimu$p975, lty = 2)
  text(FoFeSimu$Fe, FoFeSimu$Fo, labels = shortnames, cex = .5, pos = 3, col = "darkgray")
  dev.off() # Close the EMF device

  # 3. Save as PNG
  png(paste(sdirfigs, "Fo_Fe_taxa_", SafeSiteSetName, ".png", sep = ""), width = 800, height = 600, res = 100)
  # Re-plot for PNG
  plot(FoFeSimu$Fe, FoFeSimu$Fo,
    main = paste(SiteSetName, " (N=", dim(capture_probs)[1], ")", sep = ""),
    ylab = "Observed frequence", xlab = "Expected frequence", ylim = c(0, 1), xlim = c(0, 1),
    las = 1, box(bty = "l"), pch = 19, cex = 1.3, cex.lab = 1.3, cex.axis = 1.3,
    frame.plot = FALSE, col = rgb(0, 0.34509803921568627, 0.32941176470588235)
  )
  lines(FoFeSimu$Fe, FoFeSimu$p025, lty = 2)
  lines(FoFeSimu$Fe, FoFeSimu$p975, lty = 2)
  text(FoFeSimu$Fe, FoFeSimu$Fo, labels = shortnames, cex = .5, pos = 3, col = "darkgray")
  dev.off() # Close the PNG device

  print("Plots saved successfully!")
}
