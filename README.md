# TFM

### Easy QC:

EasyQC was applied to each original dataset obtained from the different participating cohorts, separating between X and autosomal chromosomes. AT EasyQC contains the script used to apply the quality control step to the chromosome X cohorts (ARIC, GAIT2, retroveCASES, retroveCONTROLS, CHS, HVH1 and MARTHA). All the scripts used were created from the script given by the creators of the package on their website (https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/software/) and following the EasyQC manual (https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf).


### METAL Scripts

Different METAL scripts were used to meta-analyse each phenotype 


### Manhattan.R

Manhattan.R scripts were used to create the Manhattan Plots for each of the phenotypes analysed along with the table that include the significant SNPs in 1MB regions. A custom function called toptab() is used in some of them. This function takes as input the autosomal and X chromosome results from METAL, and returns a table with the significant variants plus a table to be used in the creation of Regional Plots. 


### Toptab R function

Contains the toptab() custom function used in scripts used to create Mahnattan plots and significant variants tables.


### Regional plots:

LocusZoom.R contains the script used to create the significant variants' regional plots obtained in the meta-analyses. 


### Supplementary Table 1:

Autosomal Antithrombin quality control step output.


### Supplementary Tables 2-5:

HaploR complete query results, from AT, PC, PST and PSF phenotypes, respectively.
