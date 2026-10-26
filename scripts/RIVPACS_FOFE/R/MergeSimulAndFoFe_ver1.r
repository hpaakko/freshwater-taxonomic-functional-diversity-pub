# Function to merge simulated confidence limits and Fo and Fe

MergeFoFe <- function(FeConf, FoFe, SiteSetName){


print("")
print("*3. Merge FoFe function")
print("")

FoFeSimu <- merge(FeConf$Fosimu_conflimits, FoFe, by="row.names")
FoFeSimu <- as.data.table(FoFeSimu)
# Sort in decreasing order of Fe
FoFeSimu <- FoFeSimu[order(-Fe)]
FoFeSimu$Df <- FoFeSimu$Fo-FoFeSimu$Fe
# Assuming perfect detection, the difference between observed (Fo) and expected (Fe)
# frequencies, Df , measures the degree to which a species’
# actual frequency of collection in samples collected within a
# region agrees with that expected under reference conditions.

FoFeSimu$nulltest <- ifelse(FoFeSimu$Fo >= FoFeSimu$p025 & FoFeSimu$Fo <= FoFeSimu$p975, "between", "fail")
# Values of Fo that were outside of these confidence limits failed the null hypothesis test.

# str(FoFeSimu)

print(SiteSetName)
print("Mean and SD of difference between observed (Fo) and expected (Fe) frequencies = Df")

print(mean(FoFeSimu$Df))
print(sd(FoFeSimu$Df))

#hist(FoFeSimu$Df)

print("Taxa that fail null test i.e. whose frequency differ from expected:")
print(FoFeSimu[nulltest == "fail"])
print("The output is a data.table sorted in decreasing order of Fe")
FoFeSimu
} 
# end of function
