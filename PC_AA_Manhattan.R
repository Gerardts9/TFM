require (data.table)
require (qqman)
library(dplyr)
library(plotrix)
setwd("/home/gerard/Escritorio/")

results <-fread ("./MetalScripts/PC_AA/PC_AA_200720_1.txt")

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

png("/home/gerard/Escritorio/MetalScripts/PC_AA/PC_AA_200720.png")
manhattan(FilteredLZ, chr ="chr", bp ="BP", p ="P-value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PC_AA_Manhattan_GWAS")
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
toptable[,"ID"] <- c("rs200045749","rs116422036","rs867186","rs116524954")
toptable

fwrite(toptable, file = "/home/gerard/Escritorio/MetalScripts/PC_AA/toptable_PC_AA_200720.txt", sep = "\t")

png("/home/gerard/Escritorio/MetalScripts/PC_AA/PC_AA_200720v2.png")
manhattanplot(as.data.frame(FilteredLZ), chr = "chr", pos = "BP", pval = "P-value")
dev.off()

class(FilteredLZ$BP[1])
head(FilteredLZ)



###############
# Multi-Etnic #
###############

require (data.table)
require (qqman)
library(dplyr)

results <-fread ("/home/gerard/PC_AA/PC_ANT_EUR_AA_meta_200722_1.txt")
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

png("/home/gerard/PC_AA/PC_ANT_EUR_AA_200722.png")
manhattan(FilteredLZ, chr ="chr", bp ="BP", p ="P-value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PC_ANT_EUR_AA_Manhattan_GWAS")
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

toptable <- toptable[1:6,]
toptable[,"ID"] <- c("rs12740374","rs1799809","rs1260326","rs116422036","rs55747707","rs867186")
toptable

fwrite(toptable, file = "/home/gerard/PC_AA/toptable_PC_ANT_EUR_AA_200722.txt", sep = "\t")

FilteredLZ <- as.data.frame(FilteredLZ)

png("/home/gerard/PC_AA/PC_ANT_EUR_AA_200722v2.png")
manhattanplot(FilteredLZ, chr = "chr", pos = "BP", pval = "P-value")
dev.off()



###
# PC_AA_EUR:
###
PC <- fread("/home/gerard/Escritorio/MetalScripts/PC_AA/PC_EUR_AA_meta_200722_1.txt")
PC[PC$MarkerName%in%"chr2:38750551:G:C"]

results$`P-value2` <- as.numeric(results$`P-value`)
nrow(results[results$`P-value2` < 5e-08,])


results <-fread ("/home/gerard/PC_AA/PC_EUR_AA_meta_200722_1.txt")
resultsNoX <- results[-grep("X", results$MarkerName),]
if(nrow(resultsNoX) == 0){resultsNoX <- results}
resultsTotal <- resultsNoX

#resultsX <- fread (fileX)
#resultsTotal <- rbind(resultsNoX,resultsX)
#resultsTotal$MarkerName <- gsub("X",23,resultsTotal$MarkerName)      

resultsTotal$`P-value2` <- as.numeric(resultsTotal$`P-value`)

if (sum(resultsTotal$`P-value2` < 5e-324) == 0){
  resultsTotal <- resultsTotal[, !"P-value2"]
} else {stop("Results contain too small p-values")}

sum(resultsTotal$`P-value2` == 0)


res <- resultsTotal[which(resultsTotal$`P-value2` == 0),]

topsnp <- resultsTotal[resultsTotal$MarkerName %in% res[mpfr(res$`P-value`) == min(mpfr(res[["P-value"]]))]$MarkerName]    
topsnp

resultsTotal <- resultsTotal[!resultsTotal$`P-value2` == 0]

FilteredLZ <- resultsTotal[resultsTotal$`P-value2` < 0.01,]
FilteredLZ <- na.omit(FilteredLZ)

FilteredLZ <- mutate(FilteredLZ, chr = as.numeric(substr(MarkerName, 4,regexpr(":",MarkerName)-1)), 
                     BP = as.numeric(substr (substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)),1,regexpr(":",substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)))-1)))

Filtered <- FilteredLZ[FilteredLZ$`P-value2` < 5e-08,]

png("/home/gerard/PC_AA/PC_EUR_AA_200722.png")
manhattan(FilteredLZ, chr ="chr", bp ="BP", p ="P-value2", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="AT_EUR_Manhattan_GWAS")
dev.off()

FilteredOrder <- Filtered[order(Filtered$chr,Filtered$BP),]

# Loop to obtain the TOPSNPs from the results of a meta-analysis
toptable <- c()
for (i in 1:23){
  RESULTSCHR <- FilteredOrder[FilteredOrder$chr == i, ]
  while(nrow(RESULTSCHR) > 0){
    minPv <- min(RESULTSCHR$`P-value2`)
    range <- c(RESULTSCHR[RESULTSCHR$`P-value2` == minPv , ]$BP - 500000 ,RESULTSCHR[RESULTSCHR$`P-value2` == minPv , ]$BP + 500000)
    if (length(range) == 4){
      range <- range[c(1,3)]
    } else {range <- range}
    toptable <- rbind(toptable, RESULTSCHR[RESULTSCHR$`P-value2` == minPv , ])
    RESULTSCHR <- RESULTSCHR[RESULTSCHR$BP < range[1] | RESULTSCHR$BP > range[2],]
  }
}
toptable
topsnp <- mutate(topsnp, chr = as.numeric(substr(MarkerName, 4,regexpr(":",MarkerName)-1)), 
                 BP = as.numeric(substr (substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)),1,regexpr(":",substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)))-1)))

toptable <- rbind(toptable, topsnp)

toptable <- toptable[c(1:5,14),]
toptable[,"ID"] <- c("rs12740374","rs1799809","rs4665972","rs116422036","rs34594435","rs11907011")
toptable

fwrite(toptable, file = "/home/gerard/PC_AA/toptable_PC_EUR_AA_200722.txt", sep = "\t")


FilteredLZ <- as.data.frame(FilteredLZ)
head(FilteredLZ)

png("/home/gerard/PC_AA/PC_EUR_AA_200722v2.png")
manhattanplot(FilteredLZ, chr = "chr", pos = "BP", pval = "P-value")
dev.off()