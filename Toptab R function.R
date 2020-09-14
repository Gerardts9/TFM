toptab <- function(fileaut, fileX){
  require (data.table)
  require (qqman)
  library(dplyr)
  
  results <-fread (fileaut)
  resultsNoX <- results[-grep("X", results$MarkerName),] # Sexual chromosome data has to be analysed individually.
  if(nrow(resultsNoX) == 0){resultsNoX <- results} # If autosomic file did not contain X chromosome data, we restore the original data.
  
  resultsX <- fread (fileX)
  resultsTotal <- rbind(resultsNoX,resultsX)
  resultsTotal$MarkerName <- gsub("X",23,resultsTotal$MarkerName)   
  
  FilteredLZ <- resultsTotal[resultsTotal$`P-value` < 0.01,] # Reduce the number of variants to facilitate the creation of the plot.
  FilteredLZ <- na.omit(FilteredLZ)
  
  FilteredLZ <- mutate(FilteredLZ, chr = as.numeric(substr(MarkerName, 4,regexpr(":",MarkerName)-1)), # Create two new columns, chr and BP from existing columns.
                       BP = as.numeric(substr (substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)),1,regexpr(":",substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)))-1)))
  
  Filtered <- FilteredLZ[FilteredLZ$`P-value` < 5e-08,] # Select significant variants.
  
  manhattan(FilteredLZ, chr ="chr", bp ="BP", p ="P-value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
            main="Manhattan_GWAS")
  
  FilteredOrder <- Filtered[order(Filtered$chr,Filtered$BP),]
  
  # Loop to obtain the top SNPs from the results of a meta-analysis for each chromosome:
  toptable <- c()
  for (i in 1:23){
    RESULTSCHR <- FilteredOrder[FilteredOrder$chr == i, ]
    while(nrow(RESULTSCHR) > 0){
      minPv <- min(RESULTSCHR$`P-value`)
      range <- c(RESULTSCHR[RESULTSCHR$`P-value` == minPv , ]$BP - 500000 ,RESULTSCHR[RESULTSCHR$`P-value` == minPv , ]$BP + 500000)
      if (length(range) == 4){ # In case the same p-value is shared between different variants.
        range <- range[c(1,3)]
      } else {range <- range}
      toptable <- rbind(toptable, RESULTSCHR[RESULTSCHR$`P-value` == minPv , ])
      RESULTSCHR <- RESULTSCHR[RESULTSCHR$BP < range[1] | RESULTSCHR$BP > range[2],]
    }
  }
  
  return(list(toptable,FilteredLZ))}
