# Set version control
adate <- Sys.Date()
atime <- Sys.time()
sinfo <- sessionInfo()
# Create output dirs for the session and write session info to txt
sdir <- paste("output/", sesname, "/", sep = "") 
sdirg <- paste("output/", sesname, sep = "") 
sdirfigs <- paste("figs/", sesname, "/", sep = "") 
dir.create(sdir)
dir.create(sdirfigs)
# Write session info to txt
writeLines(capture.output(sinfo), paste(sdir, "RSessionInfo.txt", sep=""))
