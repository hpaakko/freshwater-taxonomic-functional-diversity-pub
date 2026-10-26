# Ver 1.0
# 4.9.2017 15:35:46
# Jukka Aroviita, with help of Kaisa Huttunen
# Hawkins & Yuan 2016 Freshwater Science
# Taxa simulation based on predictive model p:s 
# to test the null hypothesis that Fo = Fe 

# example run:
# nsimu <- 5000
# pth <- 0.5   (this may need some research!)
# captureprobs <- loocvabu$Capture.Probs
# Sosimu <- simufunk(captureprobs=OEassesscalOB$Capture.Probs, nsimu=nsimu, pth=pth)

# Fo simulation function

simufunk <- function(captureprobs, nsimu, pth, SiteSetName){

print("")
print("*1. Simulation function")
print("")

# Include only capture probabilities of those taxa predicted with p>p-threshold (pth)
# I.e., set zeros to all p < pth:
captureprobs[captureprobs<pth] <- 0  

# Set up size of data:
ntaxa <- ncol(captureprobs)
nsites <- nrow(captureprobs)
ncells <- length(captureprobs)

Sosimu <- matrix(NA, nrow=nsimu, ncol=ntaxa)
colnames(Sosimu) <- colnames(captureprobs)

# Make matrix from capture p:s:
captureprobsm <- as.matrix(captureprobs, nrow=nsites, ncol=ntaxa)

for(i in 1:nsimu){

# Random number (probability pi) between 0 and 1 for each taxon and each site:
simu_p_all <- matrix(runif(ncells), nrow=nsites, ncol=ntaxa)    

# p/a matrix if the random number pi is < predicted p (1=smaller, 0=larger):
# = Matrix of p/a of those taxa that are simulated to be observed:
simu_pa <- ifelse(simu_p_all<captureprobs, yes=1, no=0)
#str(simu_pa);sum(simu_pa)

# Calculation of So:
# Hawkinsin & Yuanin 2016 FWS jutussa: "We then computed a value of So (i.e., the sum of pi for a given taxon over all sampled sites) based on these simulated data." Misleading!! Kuvaus vähintään harhaanjohtava.
# Yuan mail 1/9/2017: that phrase in the parentheses is confusing.  The value of So is the number of sites at which the taxon is designated as being present in procedure described in the preceding sentence: “That is, for each taxon and each site, we drew a random number between 0 and 1, and if this random number was <p, we designated that taxon as being present at the site.”
# So = the number of sites at which the taxon is designated as
#  being present by the simulation
Sosimu[i,1:ntaxa] <- colSums(simu_pa)

}

# Calculation of frequency Fo=So/N:
Fosimu <- Sosimu/nsites

# Calculation of confidence limites (2.5%, 97,5% from So/N):
#   "The 2.5th and 97.5th percentiles of So/N 
#   then provided the 95% confidence limits for 
#   Fo in the case in which Fo = Fe."
Fosimu_conflimits <- apply(Fosimu,2,quantile,probs = c(0.025, 0.5, 0.975))
str(Fosimu_conflimits)
nimet <- names(Fosimu_conflimits[,2])
Fosimu_conflimits <- data.frame(t(Fosimu_conflimits))
names(Fosimu_conflimits)<-nimet
Fosimu_conflimits <- Fosimu_conflimits[order(Fosimu_conflimits[,2], decreasing=T),]
names(Fosimu_conflimits) <- c("p025", "p500", "p975")

print(paste("Simulation function finished for", SiteSetName, "with", nsimu, "simulations"))
print(paste(ntaxa, "taxa")); print(paste(nsites, "sites"))

list(Sosimu=Sosimu, Fosimu=Fosimu, Fosimu_conflimits=Fosimu_conflimits)

} 
# end of function


# Yuan mail 1.9.2017: Dear Jukka, I can see how that phrase in the parentheses is confusing.  The value of So is the number of sites at which the taxon is designated as being present in procedure described in the preceding sentence: 
#
#“That is, for each taxon and each site, we drew a random number between 0 and 1, and if this random number was <p, we designated that taxon as being present at the site.”
#
#Se is the sum of the pi over all sites for a given taxon.  Over many iterations, the algorithm yields a distribution of So for each Se.  
#
#Sorry about the confusing text! Please let me know if you have further questions. Lester

