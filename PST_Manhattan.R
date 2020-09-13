require (data.table)
require (qqman)
library(dplyr)
setwd("/home/gerard/Escritorio/")

#####
#OLD# 15/05/20
#####

###
# ALL
###
PSTALL <- toptab("/home/mariasabater/CHARGE/easyQC/PST/meta/PSTall_meta_200515_1.TBL.gz","/home/gerard/PS/PSTX/PST_EUROPE_X1.txt")

PSTALL <- toptab("./MetalScripts/PS_EUROPE_X/PST/PSTall_meta_200515_1.TBL.gz","./MetalScripts/PS_EUROPE_X/PSTX/PST_EUROPE_X1.txt")

toptables <- as.data.frame(PSTALL[1])
toptables[,"ID"] <- c("rs1350779","rs150611042")
toptables

#fwrite(toptables, file = "/home/gerard/PSFX/toptable.txt", sep = "\t")

fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/toptablePST_All_EUR.txt", sep = "\t")

Manhattan <- as.data.frame(PSTALL[2])
png("/home/gerard/Escritorio/Manhattan Plots/PST/PST_All_EUR.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PST_ALL_EUR_Manhattan_GWAS")
dev.off()


###
# ANT
###
#PSTANT <- toptab("/home/mariasabater/CHARGE/easyQC/PST/meta/PSTant_meta_200515_1.TBL.gz","/home/gerard/PS/PSTX/PST_ANT_EUROPE_X1.txt")

PSTANT <- toptab("./MetalScripts/PS_EUROPE_X/PST/PSTant_meta_200515_1.TBL","./MetalScripts/PS_EUROPE_X/PSTX/PST_ANT_EUROPE_X1.txt")

toptables <- as.data.frame(PSTANT[1])
toptables[,"ID"] <- c("rs76852077")
toptables

#fwrite(toptables, file = "/home/gerard/PSFX/toptable.txt", sep = "\t")

fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/toptablePST_Ant_EUR.txt", sep = "\t")

Manhattan <- as.data.frame(PSTANT[2])
png("/home/gerard/Escritorio/Manhattan Plots/PST/PST_Ant_EUR.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PST_ANT_EUR_Manhattan_GWAS")
dev.off()


###
# ACT
###
#PSTACT <- toptab("/home/mariasabater/CHARGE/easyQC/PST/meta/PSTact_meta_200515_1.TBL.gz","/home/gerard/PS/PSTX/PST_ACT_EUROPE_X1.txt")

PSTACT <- toptab("./MetalScripts/PS_EUROPE_X/PST/PSTact_meta_200515_1.TBL","MetalScripts/PS_EUROPE_X/PSTX/PST_ACT_EUROPE_X1.txt")

toptables <- as.data.frame(PSTACT[1])
toptables[,"ID"] <- c("rs150611042")
toptables

#fwrite(toptables, file = "/home/gerard/PSFX/toptable.txt", sep = "\t")

fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/toptablePST_Act_EUR.txt", sep = "\t")

Manhattan <- as.data.frame(PSTACT[2])
png("/home/gerard/Escritorio/Manhattan Plots/PST/PST_Act_EUR.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PST_ACT_EUR_Manhattan_GWAS")
dev.off()


#####
#NEW# 16/06/20
#####

###
# ALL
###

#PSTALL <- toptab("/home/mariasabater/CHARGE/easyQC/PST/meta/PSTall_meta_200616_1.TBL","/home/gerard/PS/PSTXNew/PST_ACT_EUROPE_X1.txt")

PSTALL <- toptab("./MetalScripts/PS_EUROPE_X/PSTNew/PSTall_meta_200616_1.TBL","./MetalScripts/PS_EUROPE_X/PSTX/PST_EUROPE_X1.txt")

toptables <- as.data.frame(PSTALL[1])
toptables[,"ID"] <- c("rs150611042")
toptables

fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/toptablePST_All_EUR.txt", sep = "\t")

Manhattan <- as.data.frame(PSTALL[2])
png("/home/gerard/Escritorio/Manhattan Plots/PSTNew/PST_All_EUR.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PST_ALL_EUR_Manhattan_GWAS")
dev.off()



###
# ANT
###
#PSTANT <- toptab("/home/mariasabater/CHARGE/easyQC/PST/meta/PSTant_meta_200616_1.TBL","/home/gerard/PS/PSTXNew/PST_ANT_EUROPE_X1.txt")

PSTANT <- toptab("./MetalScripts/PS_EUROPE_X/PSTNew/PSTant_meta_200616_1.TBL","./MetalScripts/PS_EUROPE_X/PSTXNew/PST_ANT_EUROPE_X1.txt")

toptables <- as.data.frame(PSTANT[1])
toptables[,"ID"] <- c("rs10213226")
toptables

#fwrite(toptables, file = "/home/gerard/PSFX/toptable.txt", sep = "\t")

fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/toptablePST_Ant_EUR.txt", sep = "\t")

Manhattan <- as.data.frame(PSTANT[2])
png("/home/gerard/Escritorio/Manhattan Plots/PSTNew/PST_Ant_EUR.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PST_ANT_EUR_Manhattan_GWAS")
dev.off()


###
# ACT
###
#PSTACT <- toptab("/home/mariasabater/CHARGE/easyQC/PST/meta/PSTact_meta_200616_1.TBL","/home/gerard/PS/PSTXNew/PST_ACT_EUROPE_X1.txt")

PSTACT <- toptab("./MetalScripts/PS_EUROPE_X/PSTNew/PSTact_meta_200616_1.TBL","MetalScripts/PS_EUROPE_X/PSTXNew/PST_ACT_EUROPE_X1.txt")

PSTACT <- fread("./MetalScripts/PS_EUROPE_X/PSTNew/PSTact_meta_200616_1.TBL")
PSTACT[PSTACT$`P-value` < 5e-08,]
PSTACT[PSTACT$MarkerName%in%("chr9:114321523:C:A")]


toptables <- as.data.frame(PSTACT[1])
toptables[,"ID"] <- c("rs116994374")
toptables

fwrite(toptables, file = "/home/gerard/PSFX/toptable.txt", sep = "\t")

fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/toptablePST_Act_EUR.txt", sep = "\t")

Manhattan <- as.data.frame(PSTACT[2])
png("/home/gerard/Escritorio/Manhattan Plots/PSTNew/PST_Act_EUR.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PST_ACT_EUR_Manhattan_GWAS")
dev.off()


#####
#NEW2# 30/06/20
#####

###
# ANT
###

# 30/06/20: Sense HVH1 cambiar-ho al pdf i excel:

PSTANT <- toptab("/home/gerard/PS/PST/PSTant_meta_200630_1.TBL","/home/gerard/PS/PSTXNew/PST_ANT_EUROPE_X_3006201.txt")

PSTANT <- fread("/home/gerard/Escritorio/MetalScripts/PS_EUROPE_X/PSTNew/PSTant_meta_200630_1.TBL")
PSTANT[PSTANT$HetPVal < 1e-06]
PSTANT[PSTANT$MarkerName%in%("chr9:114321523:C:A"),]

#PSTANT <- toptab("./MetalScripts/PS_EUROPE_X/PSTNew/PSTant_meta_200616_1.TBL","./MetalScripts/PS_EUROPE_X/PSTXNew/PST_ANT_EUROPE_X1.txt")

toptables <- as.data.frame(PSTANT[1])
toptables[,"ID"] <- c("rs10213226")
toptables

fwrite(toptables, file = "/home/gerard/PS/toptablePST_Act_300620.txt", sep = "\t")

#fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/toptablePST_Ant_EUR.txt", sep = "\t")

Manhattan <- as.data.frame(PSTANT[2])
png("/home/gerard/PS/PST/PST_Ant_EUR_300620.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PST_ANT_EUR_Manhattan_GWAS")
dev.off()


###
# ALL
###

# 30/06/20: Sense HVH1 cambiar-ho al pdf i excel:

#PSTALL <- toptab("/home/gerard/PS/PST/PSTall_meta_200630_1.TBL","/home/gerard/PS/PSTXNew/PST_EUROPE_X_2006301.txt")

PSTALL <- toptab("/home/gerard/Escritorio/MetalScripts/PS_EUROPE_X/PSTNew/PSTall_meta_200630_1.TBL","/home/gerard/Escritorio/MetalScripts/PS_EUROPE_X/PSTXNew/PST_EUROPE_X_2006301.txt")


PSTALL <-fread ("/home/gerard/Escritorio/MetalScripts/PS_EUROPE_X/PSTNew/PSTall_meta_200630_1.TBL")
PSTALL[PSTALL$MarkerName%in%"chr2:38750551:G:C",]

toptables <- as.data.frame(PSTALL[1])
toptables[,"ID"] <- c("rs150611042")
toptables

fwrite(toptables, file = "/home/gerard/PS/toptablePST_All_300620.txt", sep = "\t")

fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/toptablePST_All_EUR.txt", sep = "\t")


Manhattan <- as.data.frame(PSTALL[2])
png("/home/gerard/PS/PST/PST_All_EUR_300620.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="PST_ALL_EUR_Manhattan_GWAS")
dev.off()
