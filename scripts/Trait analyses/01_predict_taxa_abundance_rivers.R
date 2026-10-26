# This R script uses a custom function to predict taxon abundances for new sites.

# -----------------------------------------------------------------------------
# STEP 1: Load the necessary data and model components.
#
# For a full RIVPACS model, you would typically have:
#   - bugcal.pa: a matrix of presence/absence data for your calibration sites.
#   - bugcal.abu: a matrix of abundance data for your calibration sites.
#   - grps.final: a vector of group memberships for your calibration sites.
#   - preds.final: a vector of the names of your predictor variables.
#   - ranfor.mod: the trained Random Forest model object.
#
# For your new test sites, you will need:
#   - prednew: a data frame of environmental predictors for the new sites.
#   - bugnew: a matrix of observed abundances for the new sites.
# -----------------------------------------------------------------------------

library(vegan)
library(cluster)
library(randomForest)
library(data.table)
library(ggplot2)
library(here)

rm(list = ls())

trait_dir <- here("scripts", "Trait analyses")

setwd(trait_dir)

load(here("scripts/RIVPACS_FOFE/data/Strict_Henna_RIVPACS_data_rivers.RData"))

bugcal.pa <- BIO.cal.pa
bugcal.abu <- BIOcal
grps.final <- grps.finalAB2
preds.final <- preds.finalAB2
ranfor.mod <- rf.modAB

prednew <- ENVall
bugnew <- BIOall

# -----------------------------------------------------------------------------
# STEP 2: Source the R function file.
#
# You MUST change the path below to the location where you saved the function
# from the previous conversation.
# -----------------------------------------------------------------------------
source("model.predict.RanFor.4.2_and_abuOE.r")

# -----------------------------------------------------------------------------
# STEP 3: Run the prediction function on your test data.
#
# Pc = 0.5 is a common threshold for capture probability.
# Cal.OOB = FALSE is used for new, non-calibration sites.
# -----------------------------------------------------------------------------
prediction_results <- model.predict.RanFor.4.2_abu(
  bugcal.abu = bugcal.abu,
  bugcal.pa = bugcal.pa,
  grps.final = grps.final,
  preds.final = preds.final,
  ranfor.mod = ranfor.mod,
  prednew = prednew,
  bugnew = bugnew,
  Pc = 0.5,
  Cal.OOB = FALSE
)

# -----------------------------------------------------------------------------
# STEP 5: Inspect and save the results.
#
# The `prediction_results` object is a list. The expected abundances are
# stored in the `Exp.Abus` element of this list.
# ----------------------------------------------------dim(-------------------------

# Print a summary of the OE scores
cat("\nSummary of OE Scores for new sites:\n")
print(summary(prediction_results$OE.scores))

# Print the matrix of expected abundances (Exp.Abus)
cat("\nExpected Abundances for each Taxon at each Test Site:\n")
print(prediction_results$Exp.LogAbus)

# To save the expected abundances to a CSV file:
write.csv(prediction_results$Exp.LogAbus, file.path(here("output/expected_abundances_rivers.csv")), row.names = TRUE)

write.csv(prediction_results$Exp.Abus, file.path(here("output/expected_abundances_rivers_raw.csv")), row.names = TRUE)
