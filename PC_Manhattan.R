require (data.table)
require (qqman)
library(dplyr)
setwd("/home/gerard/Escritorio/")

###
# ALL
###

#02/07/20:

#PCALL <- toptab("/home/mariasabater/CHARGE/easyQC/PC/meta/PCall_meta_200702_1.TBL","/home/gerard/PC/PC_EUROPE_X1.txt")
#PCALL <- toptab("./MetalScripts/PC_EUROPE/PCall/PCall_meta_200702_1.TBL","./MetalScripts/PC")

require (data.table)
require (qqman)
library(dplyr)
library(Rmpfr)
library(tidyr)

PC <- fread ("/home/gerard/Escritorio/MetalScripts/PC_EUROPE/PCALL/PCall_meta_200702_1.TBL")

PC[PC$MarkerName%in%"chr1:109274968:G:T",]
PC[PC$MarkerName%in%"chr2:127418299:G:A",]
PC[PC$MarkerName%in%"chr2:27375230:T:C",]
PC[PC$MarkerName%in%"chr2:127978818:T:C",]
PC[PC$MarkerName%in%"chr7:73562919:C:T",]
PC[PC$MarkerName%in%"chr20:35179967:C:T",]



#results <- fread ("/home/mariasabater/CHARGE/easyQC/PC/meta/PCall_meta_200702_1.TBL")
resultsNoX <- results[-grep("X", results$MarkerName),]
if(nrow(resultsNoX) == 0){resultsNoX <- results}

#resultsX <- fread ("/home/gerard/PC/PC_EUROPE_X1.txt")
resultsX <- fread ("./MetalScripts/PC_EUROPE/PCALLX/PC_EUROPE_X1.txt")
resultsTotal <- rbind(resultsNoX,resultsX)
resultsTotal$MarkerName <- gsub("X",23,resultsTotal$MarkerName)

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

manhattan(FilteredLZ, chr ="chr", bp ="BP", p ="P-value2", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="AT_EUR_Manhattan_GWAS")

FilteredOrder <- Filtered[order(Filtered$chr,Filtered$BP),]

# Loop to obtain the TOPSNPs from the results of a meta-analysis
toptable <- c()
for (i in c(1:19,21,22,23)){
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


toptable
toptable[,"ID"] <- c("rs599839","rs1799809","rs4665972","rs34594435","rs55707100","rs150070344","rs11907011")
toptable

#fwrite(toptable, file = "/home/gerard/PC/Toptable/PC_All_200702.txt", sep = "\t")
fwrite(toptable, file = "/home/gerard/Escritorio/Locus Zoom + toptable/PC/toptablePC_All_EUR.txt", sep = "\t")

png("/home/gerard/Escritorio/Manhattan Plots/PC/PC_All_EUR.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P-value2", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PC_ALL_EUR_Manhattan_GWAS")
dev.off()


###
# ANT
###

#02/07/20:

#PCANT <- toptab("/home/mariasabater/CHARGE/easyQC/PC/meta/PCant_meta_200702_1.TBL","/home/gerard/PC/PC_ANT_EUROPE_X1.txt")

PCANT <- toptab("./MetalScripts/PC_EUROPE/PCANT/PCant_meta_200702_1.TBL","./MetalScripts/PC_EUROPE/PCANTX/PC_ANT_EUROPE_X1.txt")

toptables <- as.data.frame(PCANT[1])
toptables[,"ID"] <- c("rs1799809","rs1260326","rs17145713","rs867186","rs764598","rs976728","rs17347958","rs541937569","rs561284372","rs193135394")
toptables

#fwrite(toptables, file = "/home/gerard/PC/Toptable/PC_Ant_200702.txt", sep = "\t")

fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/PC/toptablePC_Ant_EUR.txt", sep = "\t")

Manhattan <- as.data.frame(PCANT[2])
png("/home/gerard/Escritorio/Manhattan Plots/PC/PC_Ant_EUR.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PC_ANT_Manhattan_GWAS")
dev.off()


###
# ACT
###

#04/06/20:

#PCACT <- toptab("/home/mariasabater/CHARGE/easyQC/PC/meta/PCact_meta_200604_1.TBL.gz","/home/gerard/PC/PC_ACT_EUROPE_X1.txt")

PCACT <- toptab("./MetalScripts/PC_EUROPE/PCACT/PCact_meta_200604_1.TBL.gz","./MetalScripts/PC_EUROPE/PCACTX/PC_ACT_EUROPE_X1.txt")

toptables <- as.data.frame(PCACT[1])
toptables[,"ID"] <- c("rs1799809","rs4665972","rs11906148","rs764598","rs6088366","rs77914463","rs535571454")
toptables

#fwrite(toptables, file = "/home/gerard/PC/Toptable/PC_Act_200604.txt", sep = "\t")

fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/PC/toptablePSF_Act_EUR.txt", sep = "\t")

Manhattan <- as.data.frame(PCACT[2])
png("/home/gerard/Escritorio/Manhattan Plots/PC/PC_Act_EUR.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PC_ACT_Manhattan_GWAS")
dev.off()



###
# ALL
###

#07/07/20: -CHS-CUSHMAN.

#PCALL <- toptab("/home/mariasabater/CHARGE/easyQC/PC/meta/PCall_meta_200702_1.TBL","/home/gerard/PC/PC_EUROPE_X1.txt")
#PCALL <- toptab("./MetalScripts/PC_EUROPE/PCall/PCall_meta_200702_1.TBL","./MetalScripts/PC")

require (data.table)
require (qqman)
library(dplyr)
library(Rmpfr)
library(tidyr)

results <- fread ("/home/gerard/Escritorio/MetalScripts/PC_EUROPE/PCALL/PCall_meta_200707_1.TBL")
results[results$HetPVal < 1e-6,]
#results <- fread ("/home/mariasabater/CHARGE/easyQC/PC/meta/PCall_meta_200707_1.TBL")
resultsNoX <- results[-grep("X", results$MarkerName),]
if(nrow(resultsNoX) == 0){resultsNoX <- results}

#resultsX <- fread ("/home/gerard/PC/PC_EUROPE_X1.txt")
resultsX <- fread ("./MetalScripts/PC_EUROPE/PCALLX/PC_ALL_200707_X1.txt")
resultsTotal <- rbind(resultsNoX,resultsX)
resultsTotal$MarkerName <- gsub("X",23,resultsTotal$MarkerName)

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

manhattan(FilteredLZ, chr ="chr", bp ="BP", p ="P-value2", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="AT_EUR_Manhattan_GWAS")

FilteredOrder <- Filtered[order(Filtered$chr,Filtered$BP),]

# Loop to obtain the TOPSNPs from the results of a meta-analysis
toptable <- c()
for (i in c(1:19,21,22,23)){
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

toptable
toptable[,"ID"] <- c("rs599839","rs1799809","rs4665972","rs34594435","rs55707100","rs150070344","rs11907011")
toptable

#fwrite(toptable, file = "/home/gerard/PC/Toptable/PC_All_200702.txt", sep = "\t")
fwrite(toptable, file = "/home/gerard/Escritorio/Locus Zoom + toptable/PCNew/toptablePC_All_EUR.txt", sep = "\t")

Manhattan <- FilteredLZ
png("/home/gerard/Escritorio/Manhattan Plots/PCNew/PC_All_EUR.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P-value2", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PC_ALL_EUR_Manhattan_GWAS")
dev.off()


###
# ANT
###

#07/07/20: -CHS-CUSHMAN

#PCANT <- toptab("/home/mariasabater/CHARGE/easyQC/PC/meta/PCant_meta_200702_1.TBL","/home/gerard/PC/PC_ANT_EUROPE_X1.txt")

PCANT <- toptab("./MetalScripts/PC_EUROPE/PCANT/PCant_meta_200707_1.TBL","./MetalScripts/PC_EUROPE/PCANTX/PC_ANT_200707_X1.txt")

toptables <- as.data.frame(PCANT[1])
toptables[,"ID"] <- c("rs1799809","rs1260326","rs17145713","rs867186","rs764598","rs976728","rs17347958","rs541937569","rs561284372","rs142862312")
toptables

#fwrite(toptables, file = "/home/gerard/PC/Toptable/PC_Ant_200702.txt", sep = "\t")

fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/PCNew/toptablePC_Ant_EUR.txt", sep = "\t")

Manhattan <- as.data.frame(PCANT[2])
png("/home/gerard/Escritorio/Manhattan Plots/PCNew/PC_Ant_EUR.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PC_ANT_Manhattan_GWAS")
dev.off()