##################
#### metaUSAT ####
##################

# Load package function:
source("/home/gerard/Descargas/metaUSAT_v1.17.R")
source("/home/gerard/metaUSAT_v1.17.R")
source("C:/Users/Usuario/Downloads/metaUSAT_v1.17.R")

library(devtools)
library(data.table)
library(dplyr)


###
# AT_AA vs PST vs PSF vs PC_AA:
###

################################################
AT <- fread ("/home/gerard/Escritorio/MetalScripts/AT_AA/AT_EUR_AA_meta_2007211.txt")
PST <- fread("/home/gerard/Escritorio/MetalScripts/PS_EUROPE_X/PSTNew/PSTall_meta_200630_1.TBL")
PSF <- fread("/home/gerard/Escritorio/MetalScripts/PS_EUROPE_X/PSF/PSFall_meta_200630_1.TBL")
PC <- fread ("/home/gerard/Escritorio/MetalScripts/PC_AA/PC_EUR_AA_meta_200722_1.txt")

#AT <- fread("/home/mariasabater/CHARGE/easyQC/AT/meta/AT_meta_200703_1.TBL")
#PST <- fread("/home/mariasabater/CHARGE/easyQC/PST/meta/PSTall_meta_200630_1.TBL")
#PSF <- fread("/home/mariasabater/CHARGE/easyQC/PSF/meta/PSFall_meta_200630_1.TBL")
#PC <- fread("/home/mariasabater/CHARGE/easyQC/PC/meta/PCall_meta_200702_1.TBL")

AT <- AT[-grep("X", AT$MarkerName),]
PST <- PST[-grep("X", PST$MarkerName),]
PSF <- PSF[-grep("X", PSF$MarkerName),]
PC <- PC[-grep("X", PC$MarkerName),]

ATX <- fread("/home/gerard/Escritorio/MetalScripts/AT_EUROPE/ATX/AT_EUROPE_200703_X1.txt")
PSTX <- fread("/home/gerard/Escritorio/MetalScripts/PS_EUROPE_X/PSTXNew/PST_EUROPE_X_2006301.txt")
PSFX <- fread("/home/gerard/Escritorio/MetalScripts/PS_EUROPE_X/PSFX/PSF_EUROPE_X_2006301.txt")
PCX <- fread("/home/gerard/Escritorio/MetalScripts/PC_EUROPE/PCALLX/PC_ALL_200707_X1.txt")

#ATX <- fread("/home/gerard/AT/AT_EUROPE_200703_X1.txt")
#PSTX <- fread("/home/gerard/PS/PSTXNew/PST_EUROPE_X_2006301.txt")
#PSFX <- fread("/home/gerard/PS/PSFX/PSF_EUROPE_X_2006301.txt")
#PCX <- fread("/home/gerard/PC/PC_EUROPE_X1.txt")

AT <- rbind(AT,ATX)
PST <- rbind(PST,PSTX)
PSF <- rbind(PSF,PSFX)
PC <- rbind(PC,PCX)

AT$MAC <- ifelse(AT$Freq1 < 0.50, AT$Freq1*AT$N*2, (1-AT$Freq1)*AT$N*2)
PST$MAC <- ifelse(PST$Freq1 < 0.50, PST$Freq1*PST$N*2, (1-PST$Freq1)*PST$N*2)
PSF$MAC <- ifelse(PSF$Freq1 < 0.50, PSF$Freq1*PSF$N*2, (1-PSF$Freq1)*PSF$N*2)
PC$MAC <- ifelse(PC$Freq1 < 0.50, PC$Freq1*PC$N*2, (1-PC$Freq1)*PC$N*2)


# Apply MAC filter:
AT <- AT[AT$MAC > 30,] # 
PST <- PST[PST$MAC > 30,] # 
PSF <- PSF[PSF$MAC > 30,] # 
PC <- PC[PC$MAC > 30,] # 

common <- intersect(intersect(intersect(AT$MarkerName, PST$MarkerName), PSF$MarkerName),PC$MarkerName) # Select common SNPs between the four phenotypes.

# Create files for each chromosome to apply pruning:
for (i in 1:22){
chr <- grep(paste0("chr",i,":"), common, value = TRUE)
chr <- as.list(chr)
fwrite(chr, file = paste0("/media/gerard/CCCOMA_X86FRE_ES-ES_DV9/CHARGE/PLINK/AT_PST_PSF_PC/Chromosomes/chr",i,".txt"), sep = "\t")}


head(common)
chrx <- grep("chrX", common, value = TRUE)
chrx <- as.list(chrx)
fwrite(chrx, file = "/home/gerard/PLINK/Chromosomes/chrX.txt", sep = "\t")

# Pruning bash code:
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr1.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr1.txt --out /home/gerard/Pruning/chr1
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr2.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr2.txt --out /home/gerard/Pruning/chr2
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr3.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr3.txt --out /home/gerard/Pruning/chr3
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr4.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr4.txt --out /home/gerard/Pruning/chr4
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr5.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr5.txt --out /home/gerard/Pruning/chr5
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr6.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr6.txt --out /home/gerard/Pruning/chr6
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr7.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr7.txt --out /home/gerard/Pruning/chr7
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr8.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr8.txt --out /home/gerard/Pruning/chr8
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr9.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr9.txt --out /home/gerard/Pruning/chr9
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr10.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr10.txt --out /home/gerard/Pruning/chr10
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr11.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr11.txt --out /home/gerard/Pruning/chr11
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr12.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr12.txt --out /home/gerard/Pruning/chr12
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr13.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr13.txt --out /home/gerard/Pruning/chr13
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr14.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr14.txt --out /home/gerard/Pruning/chr14
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr15.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr15.txt --out /home/gerard/Pruning/chr15
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr16.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr16.txt --out /home/gerard/Pruning/chr16
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr17.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr17.txt --out /home/gerard/Pruning/chr17
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr18.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr18.txt --out /home/gerard/Pruning/chr18
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr19.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr19.txt --out /home/gerard/Pruning/chr19
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr20.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr20.txt --out /home/gerard/Pruning/chr20
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr21.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr21.txt --out /home/gerard/Pruning/chr21
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr22.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chr22.txt --out /home/gerard/Pruning/chr22
plink --vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr23/chrX.dose.vcf.gz --indep-pairwise 500 'kb' 50 0.1 --extract /home/gerard/PLINK/Chromosomes/chrX.txt --out /home/gerard/Pruning/chrX


AT <- AT[AT$MarkerName %in% common,]
PST <- PST[PST$MarkerName %in% common,]
PSF <- PSF[PSF$MarkerName %in% common,]
PC <- PC[PC$MarkerName %in% common,]

AT$MarkerName <- gsub("X",23,AT$MarkerName)
PST$MarkerName <- gsub("X",23,PST$MarkerName)
PSF$MarkerName <- gsub("X",23,PSF$MarkerName)
PC$MarkerName <- gsub("X",23,PC$MarkerName)

AT <- AT[order(AT$MarkerName),]
PST <- PST[order(PST$MarkerName),]
PSF <- PSF[order(PSF$MarkerName),]
PC <- PC[order(PC$MarkerName),]

identical(AT$MarkerName,PST$MarkerName)
identical(AT$MarkerName,PSF$MarkerName)
identical(AT$MarkerName,PC$MarkerName)

head(AT)
head(PST)
head(PSF)
head(PC)

# Create z-score value, required in the method.
AT$z.score <- AT$Effect/AT$StdErr
PST$z.score <- PST$Effect/PST$StdErr
PSF$z.score <- PSF$Effect/PSF$StdErr
PC$z.score <- PC$Zscore

head(AT)
head(PST)
head(PSF)
head(PC)

zmatrix <- cbind(AT$z.score,PST$z.score, PSF$z.score, PC$z.score)
colnames(zmatrix) <- c("AT","PST","PSF","PC")
pmatrix <- cbind(AT$`P-value`,PST$`P-value`, PSF$`P-value`, as.numeric(PC$`P-value`))
colnames(pmatrix) <-c("AT","PST","PSF","PC")

nrow(zmatrix)
nrow(pmatrix)

# Obtain correlation matrix:
R <- cor.pearson(Z.matrix = zmatrix , P.matrix = pmatrix, p.threshold = 1e-5)
R

common <- AT$MarkerName
length(common)
head(common)

nrow(zmatrix)

#####
# Z-statisic filter:
#zmatrix2 <- zmatrix[zmatrix[,1] < 0.15 & zmatrix[,2] < 0.15 & zmatrix[,3] < 0.15 & zmatrix[,4] < 0.15,]
#nrow(zmatrix2)
#common2 <- common[zmatrix[,1] < 0.15 & zmatrix[,2] < 0.15 & zmatrix[,3] < 0.15 & zmatrix[,4] < 0.15]
#length(common2)
#head(common2)
#####



# Apply the p-value filter:
nrow(pmatrix)
length(common)

pmatrix2 <- pmatrix[pmatrix[,1] < 0.01 | pmatrix[,2] < 0.01 | pmatrix[,3] < 0.01 | pmatrix[,4] < 0.01,]
nrow(pmatrix2)

common2 <- common[pmatrix[,1] < 0.01 | pmatrix[,2] < 0.01 | pmatrix[,3] < 0.01 | pmatrix[,4] < 0.01]
length(common2)

zmatrix2 <- zmatrix[pmatrix[,1] < 0.01 | pmatrix[,2] < 0.01 | pmatrix[,3] < 0.01 | pmatrix[,4] < 0.01,]
nrow(zmatrix2)


# Prepare data to be parallelized:
m1 <- zmatrix2[1:50000,]
m2 <- zmatrix2[50001:100000,]
m3 <- zmatrix2[100001:150000,]
m4 <- zmatrix2[150001:200000,]
m5 <- zmatrix2[200001:250000,]
m6 <- zmatrix2[250001:300000,]
m7 <- zmatrix2[300001:350000,]
m8 <- zmatrix2[350001:386969,]

mlist <- list(m1,m2,m3,m4,m5,m6,m7,m8)

sum(unlist(lapply(mlist,nrow)))

metaUSAT <- function(l){
  pvals <- c()
  for (i in 1:nrow(l)){
  pvals <- c(pvals, metausat(Z = l[i,], R, weights = 1)$p.metausat)}
return(pvals)}

require(parallel); out <- mclapply(mlist, metaUSAT, mc.cores=6)
metaUSAT_AT_PC_PS <- out
save(metaUSAT_AT_PC_PS, file = "/home/gerard/MultiFinal/AT_PC_PS/metaUSAT_AT_PC_PS.Rdata") # Save results.

load("/home/gerard/MultiFinal/AT_PC_PS/metaUSAT_AT_PC_PS.Rdata")
load("/Users/Usuario/Desktop/AT_PC_PS.Rdata")
metaUSAT_AT_PC_PS

p.val <- unlist(metaUSAT_AT_PC_PS)
p.val

min(p.val)

Manhattan <- data.frame(common2, p.val)
head(Manhattan)

#Manhattan <- Manhattan[-grep("X", Manhattan$common2),]
dim(Manhattan)

library(dplyr)
Manhattan <- mutate(Manhattan, chr = as.numeric(substr(common2, 4,regexpr(":",common2)-1)), 
                     BP = as.numeric(substr (substr(common2, regexpr(":",common2)+1 ,nchar(common2)),1,regexpr(":",substr(common2, regexpr(":",common2)+1 ,nchar(common2)))-1)))

head(Manhattan)
Manhattan[Manhattan$common2%in%"chr16:72789996:A:G",]
sum(Manhattan$p.val == 0)
Manhattan[Manhattan$p.val == 0,]
Manhattan <- Manhattan[Manhattan$p.val != 0,]
min(Manhattan$p.val)

library(qqman)
png("/home/gerard/MultiFinal/metaUSAT_AT_PC_PS.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="p.val", snp ="common2",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="metaUSAT_AT_PC_PS")
dev.off()

head(Manhattan)
png("/home/gerard/MultiFinal/AT_PC_PS/metaUSAT_AT_PC_PSv2.png")
manhattanplot(Manhattan, chr = "chr", pos = "BP", pval = "p.val")
dev.off()

head(Manhattan)

toptable <- c()
Manhattan2 <- Manhattan[Manhattan$p.val < 5e-08,]
for (i in 1:22){
  RESULTSCHR <- Manhattan2[Manhattan2$chr == i, ]
  while(nrow(RESULTSCHR) > 0){
    minPv <- min(RESULTSCHR$p.val)
    range <- c(RESULTSCHR[RESULTSCHR$p.val == minPv , ]$BP - 500000 ,RESULTSCHR[RESULTSCHR$p.val == minPv , ]$BP + 500000)
    toptable <- rbind(toptable, RESULTSCHR[RESULTSCHR$p.val == minPv , ])
    RESULTSCHR <- RESULTSCHR[RESULTSCHR$BP < range[1] | RESULTSCHR$BP > range[2],]
  }
}
toptable

fwrite(toptable, file = "/home/gerard/MultiFinal/AT_PC_PS/metaUSAT_AT_PC_PS_toptable.txt", sep = "\t")
