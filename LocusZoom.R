require (data.table)
library(dplyr)

# AT_AA:
AT_AA <-fread ("/home/gerard/Escritorio/MetalScripts/AT_AA/AT_EUR_AA_meta_2007211.txt")
AT_AA <- AT_AA[AT_AA$`P-value` < 0.01,]
AT_AA <- mutate(AT_AA, chr = as.numeric(substr(MarkerName, 4,regexpr(":",MarkerName)-1)), 
                     BP = as.numeric(substr (substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)),1,regexpr(":",substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)))-1)))
AT_AA <- na.omit(AT_AA)

#CHR1:
chr1 <- as.data.frame(AT_AA[AT_AA$chr == 1, ])
chr1$`#CHROM` <- "chr1"
chr1$BEGIN <- chr1$BP
chr1 <- chr1[order(chr1$BP),]
head(chr1)
chr1 <- chr1[c("#CHROM", "BEGIN", "MarkerName", "Allele1", "Allele2","Freq1","FreqSE","MinFreq","MaxFreq","Effect","StdErr","P-value","Direction","HetISq","HetChiSq","HetDf","N")]
fwrite(chr1, sep = "\t", file = "/home/gerard/Escritorio/LocusZoom/AT_Ethnic/rs2227624")
bgzip rs2227624 -f; tabix -b 2 -e 2 rs2227624.gz -f

#CHR2:
chr2 <- as.data.frame(AT_AA[AT_AA$chr == 2, ])
chr2$`#CHROM` <- "chr2"
chr2$BEGIN <- chr2$BP
chr2 <- chr2[order(chr2$BP),]
head(chr2)
chr2 <- chr2[c("#CHROM", "BEGIN", "MarkerName", "Allele1","Allele2","Freq1","FreqSE","MinFreq","MaxFreq","Effect","StdErr","P-value","Direction","HetISq","HetChiSq","HetDf","N")]
fwrite(chr2, sep = "\t", file = "/home/gerard/Escritorio/LocusZoom/AT_Ethnic/rs4665972")
bgzip rs4665972 -f; tabix -b 2 -e 2 rs4665972.gz -f

#CHR7:
chr7 <- as.data.frame(AT_AA[AT_AA$chr == 7, ])
chr7$`#CHROM` <- "chr7"
chr7$BEGIN <- chr7$BP
chr7 <- chr7[order(chr7$BP),]
head(chr7)
chr7 <- chr7[c("#CHROM", "BEGIN", "MarkerName", "Allele1","Allele2","Freq1","FreqSE","MinFreq","MaxFreq","Effect","StdErr","P-value","Direction","HetISq","HetChiSq","HetDf","N")]
fwrite(chr7, sep = "\t", file = "/home/gerard/Escritorio/LocusZoom/AT_Ethnic/rs13244268")
bgzip rs13244268 -f; tabix -b 2 -e 2 rs13244268.gz -f

#CHR16:
chr16 <- as.data.frame(AT_AA[AT_AA$chr == 16, ])
chr16$`#CHROM` <- "chr16"
chr16$BEGIN <- chr16$BP
chr16 <- chr16[order(chr16$BP),]
head(chr16)
chr16 <- chr16[c("#CHROM", "BEGIN", "MarkerName", "Allele1","Allele2","Freq1","FreqSE","MinFreq","MaxFreq","Effect","StdErr","P-value","Direction","HetISq","HetChiSq","HetDf","N")]
fwrite(chr16, sep = "\t", file = "/home/gerard/Escritorio/LocusZoom/AT_Ethnic/rs5471")
bgzip rs5471 -f; tabix -b 2 -e 2 rs5471.gz -f

###
#PC_AA:
###

PC_AA <- fread("/home/gerard/Escritorio/MetalScripts/PC_AA/PC_EUR_AA_meta_200722_1.txt")
PC_AA$`P-value`  <- as.numeric(PC_AA$`P-value`)
PC_AA <- PC_AA[PC_AA$`P-value` < 0.01,]
PC_AA <- mutate(PC_AA, chr = as.numeric(substr(MarkerName, 4,regexpr(":",MarkerName)-1)), 
                BP = as.numeric(substr (substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)),1,regexpr(":",substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)))-1)))
PC_AA <- na.omit(PC_AA)

#CHR1:
chr1 <- as.data.frame(PC_AA[PC_AA$chr == 1, ])
chr1$`#CHROM` <- "chr1"
chr1$BEGIN <- chr1$BP
chr1 <- chr1[order(chr1$BP),]
head(chr1)
chr1 <- chr1[c("#CHROM", "BEGIN", "MarkerName", "Allele1","Allele2","Freq1","FreqSE","MinFreq","MaxFreq","Weight","Zscore","P-value","Direction","HetISq","HetChiSq","HetDf","N")]
fwrite(chr1, sep = "\t", file = "/home/gerard/Escritorio/LocusZoom/PC_Ethnic/rs12740374")
bgzip rs12740374 -f; tabix -b 2 -e 2 rs12740374.gz -f


#CHR2:
chr2 <- as.data.frame(PC_AA[PC_AA$chr == 2, ])
chr2$`#CHROM` <- "chr2"
chr2$BEGIN <- chr2$BP
chr2 <- chr2[order(chr2$BP),]
head(chr2)
chr2 <- chr2[c("#CHROM", "BEGIN", "MarkerName", "Allele1","Allele2","Freq1","FreqSE","MinFreq","MaxFreq","Weight","Zscore","P-value","Direction","HetISq","HetChiSq","HetDf","N")]
fwrite(chr2, sep = "\t", file = "/home/gerard/Escritorio/LocusZoom/PC_Ethnic/rs1799809")
bgzip rs1799809 -f; tabix -b 2 -e 2 rs1799809.gz -f
class(chr2$`P-value`)

#CHR7:
chr7 <- as.data.frame(PC_AA[PC_AA$chr == 7, ])
chr7$`#CHROM` <- "chr7"
chr7$BEGIN <- chr7$BP
chr7 <- chr7[order(chr7$BP),]
head(chr7)
chr7 <- chr7[c("#CHROM", "BEGIN", "MarkerName", "Allele1","Allele2","Freq1","FreqSE","MinFreq","MaxFreq","Weight","Zscore","P-value","Direction","HetISq","HetChiSq","HetDf","N")]
fwrite(chr7, sep = "\t", file = "/home/gerard/Escritorio/LocusZoom/PC_Ethnic/rs34594435")
bgzip rs34594435 -f; tabix -b 2 -e 2 rs34594435.gz -f


#CHR20:
chr20 <- as.data.frame(PC_AA[PC_AA$chr == 20, ])
chr20$`#CHROM` <- "chr20"
chr20$BEGIN <- chr20$BP
chr20 <- chr20[order(chr20$BP),]
chr20 <- chr20[chr20$`P-value` != 0 ,]
head(chr20)
chr20 <- chr20[c("#CHROM", "BEGIN", "MarkerName", "Allele1","Allele2","Freq1","FreqSE","MinFreq","MaxFreq","Weight","Zscore","P-value","Direction","HetISq","HetChiSq","HetDf","N")]
fwrite(chr20, sep = "\t", file = "/home/gerard/Escritorio/LocusZoom/PC_Ethnic/rs11907011")
bgzip rs11907011 -f; tabix -b 2 -e 2 rs11907011.gz -f


###
#PST:
###

PSTALL <-fread ("/home/gerard/Escritorio/MetalScripts/PS_EUROPE_X/PSTNew/PSTall_meta_200630_1.TBL")
PSTALL <- PSTALL[PSTALL$`P-value` < 0.01,]
PSTALL <- mutate(PSTALL, chr = as.numeric(substr(MarkerName, 4,regexpr(":",MarkerName)-1)), 
                BP = as.numeric(substr (substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)),1,regexpr(":",substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)))-1)))
PSTALL <- na.omit(PSTALL)

#CHR9:
chr9 <- as.data.frame(PSTALL[PSTALL$chr == 9, ])
chr9$`#CHROM` <- "chr9"
chr9$BEGIN <- chr9$BP
chr9 <- chr9[order(chr9$BP),]
head(chr9)
chr9 <- chr9[c("#CHROM", "BEGIN", "MarkerName", "Allele1", "Allele2","Freq1","FreqSE","MinFreq","MaxFreq","Effect","StdErr","P-value","Direction","HetISq","HetChiSq","HetDf","N")]
fwrite(chr9, sep = "\t", file = "/home/gerard/Escritorio/LocusZoom/PSTALL/rs150611042")
bgzip rs150611042 -f; tabix -b 2 -e 2 rs150611042.gz -f


###
#PSFALL:
###

PSFALL <- fread("/home/gerard/Escritorio/MetalScripts/PS_EUROPE_X/PSF/PSFall_meta_200630_1.TBL")
PSFALL <- PSFALL[PSFALL$`P-value` < 0.01,]
PSFALL <- mutate(PSFALL, chr = as.numeric(substr(MarkerName, 4,regexpr(":",MarkerName)-1)), 
                 BP = as.numeric(substr (substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)),1,regexpr(":",substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)))-1)))
PSFALL <- na.omit(PSFALL)

#CHR3:
chr3 <- as.data.frame(PSFALL[PSFALL$chr == 3, ])
chr3$`#CHROM` <- "chr3"
chr3$BEGIN <- chr3$BP
chr3 <- chr3[order(chr3$BP),]
head(chr3)
chr3 <- chr3[c("#CHROM", "BEGIN", "MarkerName", "Allele1", "Allele2","Freq1","FreqSE","MinFreq","MaxFreq","Effect","StdErr","P-value","Direction","HetISq","HetChiSq","HetDf","N")]
fwrite(chr3, sep = "\t", file = "/home/gerard/Escritorio/LocusZoom/PSFALL/rs121918472")
bgzip rs121918472 -f; tabix -b 2 -e 2 rs121918472.gz -f


###
#AT_EU:
###

AT <- fread("/home/gerard/Escritorio/MetalScripts/AT_EUROPE/AT/AT_meta_200703_1.TBL")
AT <- AT[AT$`P-value` < 0.01,]
AT <- mutate(AT, chr = as.numeric(substr(MarkerName, 4,regexpr(":",MarkerName)-1)), 
                 BP = as.numeric(substr (substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)),1,regexpr(":",substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)))-1)))
AT <- na.omit(AT)
class(AT$`P-value`)

#CHR2:
chr2 <- as.data.frame(AT[AT$chr == 2, ])
chr2$`#CHROM` <- "chr2"
chr2$BEGIN <- chr2$BP
chr2 <- chr2[order(chr2$BP),]
head(chr2)
chr2 <- chr2[c("#CHROM", "BEGIN", "MarkerName", "Allele1", "Allele2","Freq1","FreqSE","MinFreq","MaxFreq","Effect","StdErr","P-value","Direction","HetISq","HetChiSq","HetDf","N")]
fwrite(chr2, sep = "\t", file = "/home/gerard/Escritorio/LocusZoom/AT_EU/rs11127048")
bgzip rs11127048 -f; tabix -b 2 -e 2 rs11127048.gz -f


###
#PC_EU:
###
PC <- fread ("/home/gerard/Escritorio/MetalScripts/PC_EUROPE/PCALL/PCall_meta_200702_1.TBL")
PC$`P-value`  <- as.numeric(PC$`P-value`)
PC <- PC[PC$`P-value` < 0.01,]
PC <- mutate(PC, chr = as.numeric(substr(MarkerName, 4,regexpr(":",MarkerName)-1)), 
                BP = as.numeric(substr (substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)),1,regexpr(":",substr(MarkerName, regexpr(":",MarkerName)+1 ,nchar(MarkerName)))-1)))
PC <- na.omit(PC)

#CHR1:
chr1 <- as.data.frame(PC[PC$chr == 1, ])
chr1$`#CHROM` <- "chr1"
chr1$BEGIN <- chr1$BP
chr1 <- chr1[order(chr1$BP),]
head(chr1)
chr1 <- chr1[c("#CHROM", "BEGIN", "MarkerName", "Allele1","Allele2","Freq1","FreqSE","MinFreq","MaxFreq","Weight","Zscore","P-value","Direction","HetISq","HetChiSq","HetDf","N")]
fwrite(chr1, sep = "\t", file = "/home/gerard/Escritorio/LocusZoom/PC_EU/rs599839")
bgzip rs599839 -f; tabix -b 2 -e 2 rs599839.gz -f

#CHR15:
chr15 <- as.data.frame(PC[PC$chr == 15, ])
chr15$`#CHROM` <- "chr15"
chr15$BEGIN <- chr15$BP
chr15 <- chr15[order(chr15$BP),]
head(chr15)
chr15 <- chr15[c("#CHROM", "BEGIN", "MarkerName", "Allele1","Allele2","Freq1","FreqSE","MinFreq","MaxFreq","Weight","Zscore","P-value","Direction","HetISq","HetChiSq","HetDf","N")]
fwrite(chr15, sep = "\t", file = "/home/gerard/Escritorio/LocusZoom/PC_EU/rs150070344")
bgzip rs150070344 -f; tabix -b 2 -e 2 rs150070344.gz -f





