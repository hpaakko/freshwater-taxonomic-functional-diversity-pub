 # Function to plot taxa Fo and Fe results and save as pdf

PlotFoFe <- function(FoFeSimu, SiteSetName, captureprobs) {

print("")
print("*4. Plot FoFe function")
print("")

pdf(paste(sdirfigs, "Fo_Fe_", SiteSetName, ".pdf", sep=""), onefile=T)
plot(FoFeSimu$Fe, FoFeSimu$Fo, main=paste("O vs E taxon freq. with simul. CI in ", SiteSetName, " (N=", dim(capture_probs)[1],")", sep=""))

# add confidence intervals for Fe from simulation
lines(FoFeSimu$Fe, FoFeSimu$p025, lty=2)
lines(FoFeSimu$Fe, FoFeSimu$p975, lty=2)
# The ellipse represents the 95% confidence intervals 
# associated with estimates of Fe.

# add taxon names
text(FoFeSimu$Fe, FoFeSimu$Fo, labels=FoFeSimu$Row.names, cex=.4, pos=3, col="darkgray")

#print(plot1)

dev.off()

}

# end of function





 
