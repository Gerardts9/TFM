##################
###### MPAT ######
##################

setwd("/home/gerard/Escritorio/")
library(MPAT)
library(data.table)


#####
data(lipids)
head(lipids)

lipids_zscore = as.matrix(lipids[,c("Zscore.HDL","Zscore.LDL","Zscore.TG","Zscore.TC")])
Sigma = cor(lipids_zscore )
Sigma
mixAda(Z.vec=lipids_zscore[1,],Sigma) # p-value for the first SNP


###
# AT vs PST vs PSF vs PC: 
###

AT <- fread("/home/gerard/Escritorio/MetalScripts/AT_EUROPE/AT/AT_meta_200703_1.TBL")
PST <- fread("/home/gerard/Escritorio/MetalScripts/PS_EUROPE_X/PSTNew/PSTall_meta_200630_1.TBL")
PSF <- fread("/home/gerard/Escritorio/MetalScripts/PS_EUROPE_X/PSF/PSFall_meta_200630_1.TBL")
PC <- fread("/home/gerard/Escritorio/MetalScripts/PC_EUROPE/PCALL/PCall_meta_200702_1.TBL")

#AT <- fread("/home/mariasabater/CHARGE/easyQC/AT/meta/AT_meta_200703_1.TBL")
#PST <- fread("/home/gerard/PS/PST/PSTall_meta_200630_1.TBL")
#PSF <- fread("/home/gerard/PS/PSF/PSFall_meta_200630_1.TBL")
#PC <- fread("/home/mariasabater/CHARGE/easyQC/PC/meta/PCall_meta_200702_1.TBL")

common <- intersect(intersect(intersect(AT$MarkerName, PST$MarkerName), PSF$MarkerName),PC$MarkerName) # Select common SNPs in the four phenotypes.
length(common)

AT <- AT[AT$MarkerName %in% common,]
PST <- PST[PST$MarkerName %in% common,]
PSF <- PSF[PSF$MarkerName %in% common,]
PC <- PC[PC$MarkerName %in% common,]


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
Sigma <- cor(zmatrix)
Sigma

common <- AT$MarkerName
length(common)


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
m8 <- zmatrix2[350001:400000,]
m9 <- zmatrix2[400001:450000,]
m10 <- zmatrix2[450001:458475,]

mlist <- list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)

sum(unlist(lapply(mlist,nrow)))

multiada <- function(l){
  pvals <- c()
  for (i in 1:nrow(l)){
    pvals <- c(pvals, mixAda(Z.vec=l[i,], Sigma))}
  return(pvals)
}

mixAda(Z.vec=zmatrix2[1,], Sigma)

require(parallel); out <- mclapply( mlist, multiada, mc.cores=6)
mixADA_AT_PC_PS <- out
save(mixADA_AT_PC_PS, file ="/home/gerard/mixADA_AT_PC_PS.Rdata")

load("/home/gerard/mixADA_AT_PC_PS.Rdata")

p.val <- unlist(mixADA_AT_PC_PS)
min(p.val)
length(p.val)

Manhattan <- data.frame(common2, p.val)
head(Manhattan)

Manhattan <- Manhattan[-grep("X", Manhattan$common2),]

library(dplyr)
Manhattan <- mutate(Manhattan, chr = as.numeric(substr(common2, 4,regexpr(":",common2)-1)), 
                    BP = as.numeric(substr (substr(common2, regexpr(":",common2)+1 ,nchar(common2)),1,regexpr(":",substr(common2, regexpr(":",common2)+1 ,nchar(common2)))-1)))

head(Manhattan)

library(qqman)
png("/home/gerard/mixADA_AT_PC_PS.png")
manhattan(Manhattan, chr ="chr", bp ="BP", p ="p.val", snp ="common2",suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
          main="mixADA_AT_PC_PS")
dev.off()


png("/home/gerard/mixADA_AT_PC_PSv2.png")
manhattanplot(Manhattan, chr = "chr", pos = "BP", pval = "p.val")
dev.off()


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

length(common)
length(common2)


# BAZ1B Gene common in PC and AT:
mixAda(Z.vec=m3[match("chr7:73490480:C:T", common2),], Sigma) # p-value = 2.373534e-10
mixAda(Z.vec=zmatrix2[match("chr7:73471480:G:A", common2),], Sigma) # p-value = 2.443714e-10


# GCKR: Common gene in PC and AT:
mixAda(Z.vec=zmatrix[match("chr2:27508073:T:C", common),], Sigma) # p-value = 2.151763e-10
mixAda(Z.vec=zmatrix[match("chr15:44601447:G:C", common),], Sigma) # p-value = 2.151763e-10


mixAda(Z.vec=zmatrix[match("chr2:127291098:A:C", common),], Sigma) # p-value = 2.151763e-10