# TFM

### Easy QC:

EasyQC was applied to each original dataset obtained from the different participating cohorts, separating between X and autosomal chromosomes. contains the autosomal Antithrombin script used to appy the quality control to the ten cohorts  

### METAL Scripts

Different METAL scripts were used to meta-analyse each phenotype 


### Manhattan.R

Manhattan.R scripts were used to create the Manhattan Plots for each of the phenotypes analysed along with the table that include the significant SNPs in 1MB regions. A custom function called toptab() is used in some of them. This function takes as input the autosomal and X chromosome results from METAL, and returns a table with the significant variants plus a table to be used in the creation of Regional Plots. 


### Toptab R function

Contains the toptab() custom function used in scripts used to create Mahnattan plots and significant variants tables.
