#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate Fo and Fe function
# Hawkins & Yuan 2016 Freshwater Science
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#captureprobs <- capture_probs
#spedatapa <- spe_data_pa


fofefunc <- function(captureprobs, spedatapa, pth, SiteSetName){

print("")
print("*2. Calculate Fo and Fe function")
print("")

# Select only capture probabilities of those taxa predicted with p>p-threshold (pth)
capture_probs_pth <- captureprobs
capture_probs_pth[capture_probs_pth<pth] <- 0       # set zeros to all P:s < pth

#plot(rowSums(capture_probs_pth)~OE_data$E, main=SiteSetName)  # check that is same as E

# remove those taxa not predicted by the model (i.e. not present in model calibration sites), i.e. select only modelled taxa
spedatapa <- spedatapa[,colnames(capture_probs_pth)]

# capture probabilities of those taxa predicted with p>pth AND observed 
capture_probs_pth_obs <- capture_probs_pth * spedatapa   

# p/a matrix of those taxa predicted (with p>pth) and observed 
pa_pth_obs <- capture_probs_pth_obs
pa_pth_obs[capture_probs_pth_obs>=pth] <- 1

#plot(rowSums(pa_pth_obs)~OE_data$O, main=SiteSetName)   # check that is same as O       

# Calculate observed frequency of taxa in the site set
Fo <- colSums(pa_pth_obs) / dim(pa_pth_obs)[1]

# Calculate expected frequency of taxa in the site set
# Fe is derived by aggregating site-specific estimates of pi
# as Fe = Si/N, where Si is the sum of pi for a given taxon 
# over all sampled sites, and N is the number of sites sampled (Hawkins & Yuan 2016) 
Fe <- colSums(capture_probs_pth) / dim(capture_probs_pth)[1]

print(paste("Fo and Fe calculated for", SiteSetName))

data.frame(Fo,Fe)
} 
# end of function

# Assuming perfect detection, the difference between observed (Fo) and expected (Fe)
# frequencies, Df , measures the degree to which a species’
# actual frequency of collection in samples collected within a
# region agrees with that expected under reference conditions (Hawkins & Yuan 2016).




