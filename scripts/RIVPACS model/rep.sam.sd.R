# Function to calculate replicate sampling SD of O/E
rep.sam.sd <- function(occprb.cal, Pc) {
  # Identify taxa above the predicted probability cutoff
  cal.cut <- (occprb.cal >= Pc)  # TRUE/FALSE matrix
  
  # Use occurrence probabilities only for taxa above the cutoff
  E.cal <- apply(occprb.cal * cal.cut, 1, sum)  # Vector of predicted E
  
  # Numerator of site-specific replicate sampling variance
  RS.cal <- apply(occprb.cal * cal.cut, 1, function(x) sum(x * (1 - x)))
  
  # Compute replicate sampling SD
  SDRS <- sqrt(mean(RS.cal / (E.cal^2)))
  
  # Print the result
  print(" ", quote = FALSE)
  print("Replicate sampling SD of O/E: ", quote = FALSE)
  print(SDRS, digits = 4)
}

# Example usage (replace with actual data)
# occprb.cal <- matrix(runif(50, 0, 1), nrow = 10, ncol = 5)  # Simulated probability matrix
# rep.sam.sd(occprb.cal, Pc = 0.5)
