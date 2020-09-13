require (data.table)
require (qqman)
library(dplyr)
setwd("/home/gerard/Escritorio/")

###
# ALL
###

# 20/07/02

AT <- toptab("/home/mariasabater/CHARGE/easyQC/AT/meta/AT_meta_200702_1.TBL","/home/gerard/AT/AT_EUROPE_200702_X1.txt")
#AT <- toptab("./MetalScripts/AT_EUROPE_X/AT/AT_meta_200702_1.TBL","./MetalScripts/AT_EUROPE_X/ATX/AT_EUROPE_200702_X1.txt")

toptables <- as.data.frame(AT[1])
toptables[,"ID"] <- c("rs557150901","rs2227624","rs1964103","rs954475162","rs11127048","rs13244268")
toptables

fwrite(toptables, file = "/home/gerard/AT/Toptable/toptable_AT_200702.txt", sep = "\t")

#fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/AT/toptable_AT_200702.txt", sep = "\t")

Manhattan <- as.data.frame(AT[2])
png("/home/gerard/AT/Manhattan/AT_200702.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="AT_EUR_Manhattan_GWAS")
dev.off()



###
# ALL
###

# 20/07/03

AT <- toptab("/home/mariasabater/CHARGE/easyQC/AT/meta/AT_meta_200703_1.TBL","/home/gerard/AT/AT_EUROPE_200703_X1.txt")


AT <- fread("./MetalScripts/AT_EUROPE/AT/AT_meta_200703_1.TBL")
AT[AT$MarkerName%in%"chr1:173914872:A:T",]

AT <- toptab("./MetalScripts/AT_EUROPE/AT/AT_meta_200703_1.TBL","./MetalScripts/AT_EUROPE/ATX/AT_EUROPE_200703_X1.txt")

toptables <- as.data.frame(AT[1])
toptables[,"ID"] <- c("rs557150901","rs2227624","rs1964103","rs954475162","rs11127048","rs13244268")
toptables

fwrite(toptables, file = "/home/gerard/AT/Toptable/toptable_AT_200703.txt", sep = "\t")

#fwrite(toptables, file = "/home/gerard/Escritorio/Locus Zoom + toptable/AT/toptable_AT_200703.txt", sep = "\t")

Manhattan <- as.data.frame(AT[2])
png("/home/gerard/AT/Manhattan/AT_200703.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="P.value", snp ="MarkerName",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="AT_EUR_Manhattan_GWAS")
dev.off()
