# A run script to run all Fo vs Fe functions


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print("Look at bug data: pa- data; predicted capture probs")
print(dim(spe_data_pa))
print(dim(capture_probs))
row.names(spe_data_pa) == row.names(capture_probs)

# Just look also at site set env data; abu data
# dim(pred_pohjis_REF); dim(bug_pohjis_REF)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Run taxa simulation based on RIVPACS-modelled p

# Read simulation function
# source("D:/Aroviita/3Rproj/1Functions/Taxon_frequencies/FoFe_simufunk_ver1.r")
# source("M:/GVESIBIO/3Rproj/1Functions/Taxon_frequencies/FoFe_simufunk_ver1.r")
source("R/FoFe_simufunk_ver1.r")

set.seed(123)
# Run the simulation function using the RIVPACS-modelled data
FeConf <- simufunk(captureprobs = capture_probs, nsimu = nsimu, pth = pth, SiteSetName = SiteSetName)
str(FeConf)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Calculate Fo

# Read Fo and Fe calculation function
# source("D:/Aroviita/3Rproj/1Functions/Taxon_frequencies/CalculateFoAndFe_ver1.r")
# source("M:/GVESIBIO/3Rproj/1Functions/Taxon_frequencies/CalculateFoAndFe_ver1.r")
source("R/CalculateFoAndFe_ver1.r")

# Run the function
FoFe <- fofefunc(captureprobs = capture_probs, spedatapa = spe_data_pa, pth = pth, SiteSetName = SiteSetName)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Merge simulated confidence limits and Fo and Fe

# Read Fo and Fe merge function
# source("D:/Aroviita/3Rproj/1Functions/Taxon_frequencies/MergeSimulAndFoFe_ver1.r")
# source("M:/GVESIBIO/3Rproj/1Functions/Taxon_frequencies/MergeSimulAndFoFe_ver1.r")
source("R/MergeSimulAndFoFe_ver1.r")

# Run the function
FoFeSimu <- MergeFoFe(FeConf = FeConf, FoFe = FoFe, SiteSetName = SiteSetName)
FoFeSimu

print(SiteSetName)
head(FoFeSimu)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Plot Fo to Fe
# Read the function
# source("D:/Aroviita/3Rproj/1Functions/Taxon_frequencies/PlotFoFe_ver1.r")
# source("M:/GVESIBIO/3Rproj/1Functions/Taxon_frequencies/PlotFoFe_ver1.r")
# source("R/PlotFoFe_ver1.r")
source("R/PlotFoFe_ver2.r") # uudet kuvamuodot 24.1.2024

# Run the function
# PlotFoFe(FoFeSimu=FoFeSimu, SiteSetName=SiteSetName, captureprobs=capture_probs)
result <- PlotFoFe_single(FoFeSimu = FoFeSimu, SiteSetName = SiteSetName, captureprobs = capture_probs)
p <- result[[1]]
increasing_percentage <- result[[2]]
decreasing_percentage <- result[[3]]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. Tallenna FoFeSimu
write.table(FoFeSimu, sep = ";", col.names = NA, file = paste(sdir, "FoFeSimu_", SiteSetName, ".csv", sep = ""))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Tallenna viel? k?ytetyt datat
write.table(spe_data_pa, sep = ";", col.names = NA, file = paste(sdir, "spe_data_pa_", SiteSetName, ".csv", sep = ""))
write.table(capture_probs, sep = ";", col.names = NA, file = paste(sdir, "capture_probs_", SiteSetName, ".csv", sep = ""))
write.table(env_data, sep = ";", col.names = NA, file = paste(sdir, "env_data_", SiteSetName, ".csv", sep = ""))
