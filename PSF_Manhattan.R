require (data.table)
require (qqman)
library(dplyr)

###
# ALL
###
#PSFALL <- toptab("/home/mariasabater/CHARGE/easyQC/PSF/meta/PSFall_meta_200610_1.TBL.gz","/home/gerard/PSFX/PS/PSF_EUROPE_X1.txt")

setwd("/home/gerard/Escritorio/")
PSFALL <- toptab("./MetalScripts/PS_EUROPE_X/PSF/PSFall_meta_200610_1.TBL","MetalScripts/PS_EUROPE_X/PSFX/PSF_EUROPE_X1.txt")

toptables <- as.data.frame(PSFALL[1])
toptables[,"ID"] <- c("rs566931451","rs528128538","rs121918472","rs558721154","rs150611042")
toptables

fwrite(toptables, file = "/home/gerard/PSFX/toptable.txt", sep = "\t")
fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/toptablePSF_All_EUR.txt", sep = "\t")

Manhattan <- as.data.frame(PSFALL[2])
png("/home/gerard/Escritorio/Manhattan Plots/PSF/PSF_All_EUR.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PSF_ALL_EUR_Manhattan_GWAS")
dev.off()


###
# ANT
###
#PSFANT <- toptab("/home/mariasabater/CHARGE/easyQC/PSF/meta/PSFant_meta_200610_1.TBL.gz","/home/gerard/PS/PSFX/PSF_ANT_EUROPE_X1.txt")

PSFANT <- toptab("./MetalScripts/PS_EUROPE_X/PSF/PSFant_meta_200610_1.TBL","MetalScripts/PS_EUROPE_X/PSFX/PSF_ANT_EUROPE_X1.txt")

toptables <- as.data.frame(PSFANT[1])
toptables[,"ID"] <- c("rs528128538","rs121918472","rs566931451","rs558721154","rs150611042 ")
toptables

#fwrite(toptables, file = "/home/gerard/PSFX/toptable.txt", sep = "\t")

fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/toptablePSF_Ant_EUR.txt", sep = "\t")

Manhattan <- as.data.frame(PSFANT[2])
png("/home/gerard/Escritorio/Manhattan Plots/PSF/PSF_Ant_EUR.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PSF_ANT_EUR_Manhattan_GWAS")
dev.off()


###
# ACT
###
#PSFACT <- toptab("/home/mariasabater/CHARGE/easyQC/PSF/meta/PSFact_meta_200604_1.TBL.gz","/home/gerard/PS/PSFX/PSF_ACT_EUROPE_X1.txt")


PSFACT <- toptab("./MetalScripts/PS_EUROPE_X/PSF/PSFact_meta_200604_1.TBL","MetalScripts/PS_EUROPE_X/PSFX/PSF_ACT_EUROPE_X1.txt")

PSFACT <- fread("./MetalScripts/PS_EUROPE_X/PSF/PSFact_meta_200604_1.TBL")
PSFACT[PSFACT$`P-value` < 5e-08,]
PSFACT[PSFACT$MarkerName%in%("chr3:90431347:C:G")]
PSFACT[PSFACT$MarkerName%in%("chr3:93868695:C:T")]
PSFACT[PSFACT$MarkerName%in%("chr3:93879306:A:G")]
PSFACT[PSFACT$MarkerName%in%("chr3:94976987:C:T")]
PSFACT[PSFACT$MarkerName%in%("chr9:114321523:C:A")]


toptables <- as.data.frame(PSFACT[1])
toptables[,"ID"] <- c("rs150611042")
toptables

#fwrite(toptables, file = "/home/gerard/PSFX/toptable.txt",   sep = "\t")

fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/toptablePSF_Act_EUR.txt", sep = "\t")

Manhattan <- as.data.frame(PSFACT[2])
png("/home/gerard/Escritorio/Manhattan Plots/PSF/PSF_Act_EUR.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PSF_ACT_EUR_Manhattan_GWAS")
dev.off()


###
# ANT
###

# 30/06/20 -HVH1: CAMBIAR AL POWER I EXCEL:

PSFANT <- toptab("/home/gerard/PS/PSF/PSFant_meta_200630_1.TBL","/home/gerard/PS/PSFX/PSF_ANT_EUROPE_X_3006201.txt")

#PSFANT <- toptab("./MetalScripts/PS_EUROPE_X/PSF/PSFant_meta_200630_1.TBL","MetalScripts/PS_EUROPE_X/PSFX/PSF_ANT_EUROPE_X1.txt")

PSFANT <- fread("/home/gerard/Escritorio/MetalScripts/PS_EUROPE_X/PSF/PSFant_meta_200630_1.TBL")
dim(PSFANT)
nrow(PSFANT[PSFANT$`P-value` < 5e-08,])
PSFANT[PSFANT$MarkerName%in%("chr9:114321523:C:A")]


toptables <- as.data.frame(PSFANT[1])
toptables[,"ID"] <- c("rs528128538","rs121918472","rs566931451","rs558721154","rs150611042 ")
toptables

fwrite(toptables, file = "/home/gerard/PS/PSF/toptablePSF_Ant_300620.txt", sep = "\t")

#fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/toptablePSF_Ant_EUR.txt", sep = "\t")

Manhattan <- as.data.frame(PSFANT[2])
png("/home/gerard/PS/PSF/PSF_Ant_EUR_300620.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PSF_ANT_EUR_Manhattan_GWAS")
dev.off()



###
# ALL
###

#30/06/20

PSFALL <- toptab("/home/gerard/PS/PSF/PSFall_meta_200630_1.TBL","/home/gerard/PS/PSFX/PSF_EUROPE_X_2006301.txt")

#PSFALL <- toptab("./MetalScripts/PS_EUROPE_X/PSF/PSFall_meta_200610_1.TBL","MetalScripts/PS_EUROPE_X/PSFX/PSF_EUROPE_X1.txt")

PSFALL <- fread("./MetalScripts/PS_EUROPE_X/PSF/PSFall_meta_200630_1.TBL")
PSFALL[PSFALL$MarkerName%in%"chr2:38750551:G:C"]

nrow(PSFALL[PSFALL$`P-value` < 5e-08,])

toptables <- as.data.frame(PSFALL[1])
toptables[,"ID"] <- c("rs566931451","rs528128538","rs121918472","rs558721154","rs150611042")
toptables

fwrite(toptables, file = "/home/gerard/PS/PSF/toptablePSF_All_300620.txt", sep = "\t")

#fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/toptablePSF_All_EUR.txt", sep = "\t")

Manhattan <- as.data.frame(PSFALL[2])
png("/home/gerard/PS/PSF/PSF_All_EUR_300620.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PSF_ALL_EUR_Manhattan_GWAS")
dev.off()
