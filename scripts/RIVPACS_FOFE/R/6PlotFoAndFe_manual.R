#plot(FoFeSimu$Fe, FoFeSimu$Fo, main=paste("O vs E taxon freq. with simul. CI in ", SiteSetName, " (N=", dim(capture_probs)[1],")", sep=""))

### Kaikki, jätä < .2 nimet pois
plot(FoFeSimu$Fe, FoFeSimu$Fo, xlab="Mallinnettu esiintymisfrekvenssi", ylab="Havaittu esiintymisfrekvenssi", sep="", las=1, bty="l", xlim=c(0,1), ylim=c(0,1))

# add confidence intervals for Fe from simulation
lines(FoFeSimu$Fe, FoFeSimu$p025, lty=2)
lines(FoFeSimu$Fe, FoFeSimu$p975, lty=2)
# The ellipse represents the 95% confidence intervals 
# associated with estimates of Fe.

# add taxon names
text(FoFeSimu[Fe>.1]$Fe, FoFeSimu[Fe>.1]$Fo, labels=FoFeSimu[Fe>.1]$Row.names, cex=.8, pos=3, col="blue")


### Vain 0-0.2
plot(FoFeSimu$Fe, FoFeSimu$Fo, xlab="Mallinnettu esiintymisfrekvenssi", ylab="Havaittu esiintymisfrekvenssi", sep="", las=1, bty="l", xlim=c(0,.3), ylim=c(0,.3))

# add confidence intervals for Fe from simulation
lines(FoFeSimu$Fe, FoFeSimu$p025, lty=2)
lines(FoFeSimu$Fe, FoFeSimu$p975, lty=2)
# The ellipse represents the 95% confidence intervals 
# associated with estimates of Fe.

# add taxon names
text(FoFeSimu[Fe>.02]$Fe, FoFeSimu[Fe>.02]$Fo, labels=FoFeSimu[Fe>.02]$Row.names, cex=.8, pos=3, col="blue")

