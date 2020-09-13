install.packages("haploR", dependencies = TRUE)
library(haploR)
require(openxlsx)

#AT:
AT <- queryHaploreg(query=c("rs2227624","rs4665972","rs13244268","rs5471"))
AT
write.xlsx(x=AT, file="/Users/Usuario/Desktop/AT.haploR.xlsx")

SERPINC1 <- queryHaploreg(query=c("rs2227624"), ldThresh = 0.15)
write.xlsx(x=SERPINC1, file="/Users/Usuario/Desktop/SERPINC1.AT.haploR.xlsx")


#PC:
PC <- queryHaploreg(query=c("rs12740374","rs1799809","rs4665972","rs116422036","rs34594435","rs11907011"))
PC

write.xlsx(x=PC, file="/Users/Usuario/Desktop/PC.haploR.xlsx")



#PST:
PST <- queryHaploreg(query=c("rs150611042"), ldThresh = 0.5)
write.xlsx(x=PST, file="/Users/Usuario/Desktop/PST.haploR.xlsx")



#PSF:
PSF <- queryHaploreg(query=c("rs121918472","rs150611042"), ldThresh = 0.5)
write.xlsx(x=PSF, file="/Users/Usuario/Desktop/PSF.haploR.xlsx")



gsub("^([^:]*:[^:]*):.*$", "\\1", AT$MarkerName)
/home/pferran/Others/sfwr/locuszoom/bin/locuszoom --metal LocusZoomAt.txt --build hg38 --pop EUR --source 1000G_Nov2014 --refsnp chr1:173914872 --flank 250kb --pvalcol P.value --ld-vcf /home/pferran/bdd/ref/ldref/proc/2imp/chr1.dose.vcf.gz
chr1:173914872