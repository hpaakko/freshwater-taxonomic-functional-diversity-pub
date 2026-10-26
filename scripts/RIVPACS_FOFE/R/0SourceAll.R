#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FRESHABIT A5
# Source all scripts 
# 1-6/2019 Aroviita/SYKE
# R/0SourceAll.R
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 0 Load libraries
source("D:/Aroviita/3Rproj/1Functions/Packages/load_libraries.R")

# 1 Set version control
# source("R/1SetVersionControl.R")

# 2 Read all data
source("R/2ReadData.R")


# 3 Prepare and read data
# Basics of bryophyte data
source("R/3BryophyteDataBasics.R")

# Read Marian malli data
source("R/3ReadMariaModelData.R")

# Prepare all sites' data for RIVPACS
source("R/3PrepareModelDataAllSites.R")

# 4 RIVPACS
# Prepare data for DF model RIVPACS
source("R/4RivpacsPrepareDataForDFmodel.R") 

# Run RIVPACS models
source("R/4RivpacsRunModels.R")

# Run RIVPACS model LOOCV
source("R/4RivpacsRunLOOCVmodel.R")

# 5 Status classification
source("R/5ClassifyOEResults.R")

# 6 Taxa capture probabilities and frequencies
source("R/6SimulateAndCompareFoAndFe.R")



# Plot data


# 7 Plot maps
source("R/7PlotMaps.R")

# Stressor plots

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 9 List of files
# At end write description and list all files
source("R/9WriteSessionDescription.R")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




# Analyses and figures
# source("R/MyScript.R")

# Marian pohjis RIVPACS LOOCV
# source("R/maria_mallin_ajo_LOOCV_loop.R")

# Marian pohjis abundance RIVPACS
# source("R/Marian_mallin_abundanssilasku.R")



