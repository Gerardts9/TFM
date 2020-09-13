require (data.table)
require (qqman)
library(dplyr)
library(plotrix)

results <-fread ("/home/gerard/Escritorio/MetalScripts/AT_AA/AT_AA_meta_200720_1.txt")
resultsNoX <- results[-grep("X", results$MarkerName),]
if(nrow(resultsNoX) == 0){resultsNoX <- results}

#resultsX <- fread (fileX)
#resultsTotal <- rbind(resultsNoX,resultsX)
#resultsTotal$MarkerName <- gsub("X",23,resultsTotal$MarkerName)

FilteredLZ <- resultsNoX[resultsNoX$`P-value` < 0.01,]
FilteredLZ <- na.omit(FilteredLZ)

FilteredLZ <- mutate(FilteredLZ, chr = as.numeric(substr(MarkerName, 4,regexpr(":",MarkerName)-1)), 
                     BP = as.numeric(substr (substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)),1,regexpr(":",substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)))-1)))

Filtered <- FilteredLZ[FilteredLZ$`P-value` < 5e-08,]

png("/home/gerard/AT_AA/AT_AA_200720.png")
manhattan(FilteredLZ, chr ="chr", bp ="BP", p ="P-value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="AT_AA_Manhattan_GWAS")
dev.off()

FilteredOrder <- Filtered[order(Filtered$chr,Filtered$BP),]

# Loop to obtain the TOPSNPs from the results of a meta-analysis
toptable <- c()
for (i in 1:23){
  RESULTSCHR <- FilteredOrder[FilteredOrder$chr == i, ]
  while(nrow(RESULTSCHR) > 0){
    minPv <- min(RESULTSCHR$`P-value`)
    range <- c(RESULTSCHR[RESULTSCHR$`P-value` == minPv , ]$BP - 500000 ,RESULTSCHR[RESULTSCHR$`P-value` == minPv , ]$BP + 500000)
    if (length(range) == 4){
      range <- range[c(1,3)]
    } else {range <- range}
    toptable <- rbind(toptable, RESULTSCHR[RESULTSCHR$`P-value` == minPv , ])
    RESULTSCHR <- RESULTSCHR[RESULTSCHR$BP < range[1] | RESULTSCHR$BP > range[2],]
  }
}

toptable
toptable[,"ID"] <- c("rs4253315","rs5471")
toptable

fwrite(toptable, file = "/home/gerard/AT_AA/toptable_AT_AA_200720.txt", sep = "\t")



###############
# Multi-Etnic #
###############

require (data.table)
require (qqman)
library(dplyr)

# --> Input File 1 : /home/gerard/Escritorio/MetalScripts/AT_EUROPE/AT/AT_meta_200703_1.TBL
# --> Input File 2 : /home/gerard/Escritorio/MetalScripts/AT_AA/AT_AA_meta_200720_1.txt

AT <-fread ("/home/gerard/Escritorio/MetalScripts/AT_AA/AT_EUR_AA_meta_2007211.txt")

resultsNoX <- results[-grep("X", results$MarkerName),]
if(nrow(resultsNoX) == 0){resultsNoX <- results}


#resultsX <- fread (fileX)
#resultsTotal <- rbind(resultsNoX,resultsX)
#resultsTotal$MarkerName <- gsub("X",23,resultsTotal$MarkerName)      

FilteredLZ <- resultsNoX[resultsNoX$`P-value` < 0.01,]
FilteredLZ <- na.omit(FilteredLZ)

FilteredLZ <- mutate(FilteredLZ, chr = as.numeric(substr(MarkerName, 4,regexpr(":",MarkerName)-1)), 
                     BP = as.numeric(substr (substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)),1,regexpr(":",substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)))-1)))

png("/home/gerard/Escritorio/MetalScripts/AT_AA/AT_EUR_AA_200721.png")
manhattan(FilteredLZ, chr ="chr", bp ="BP", p ="P-value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="AT_EUR_AA_Manhattan_GWAS")
dev.off()

FilteredOrder <- Filtered[order(Filtered$chr,Filtered$BP),]

# Loop to obtain the TOPSNPs from the results of a meta-analysis
toptable <- c()
for (i in 1:23){
  RESULTSCHR <- FilteredOrder[FilteredOrder$chr == i, ]
  while(nrow(RESULTSCHR) > 0){
    minPv <- min(RESULTSCHR$`P-value`)
    range <- c(RESULTSCHR[RESULTSCHR$`P-value` == minPv , ]$BP - 500000 ,RESULTSCHR[RESULTSCHR$`P-value` == minPv , ]$BP + 500000)
    if (length(range) == 4){
      range <- range[c(1,3)]
    } else {range <- range}
    toptable <- rbind(toptable, RESULTSCHR[RESULTSCHR$`P-value` == minPv , ])
    RESULTSCHR <- RESULTSCHR[RESULTSCHR$BP < range[1] | RESULTSCHR$BP > range[2],]
  }
}

toptable
toptable[,"ID"] <- c("rs557150901","rs2227624","rs1964103","rs954475162","rs4665972","rs13244268","rs5471","rs57678445")
toptable

fwrite(toptable, file = "/home/gerard/Escritorio/MetalScripts/AT_AA/toptable_AT_EUR_AA_200721.txt", sep = "\t")