## Lollipop-Plots for somatic variants:
## 
## We will be using the maftools package
## This package contains a ready2use function to create a lollipop plot
## If your requirements are not met, we can work on specific tasks individually
## 
## 01/05/2022


#' Firstly, we need to install the appropriate package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("maftools")

#' load the library via
library(maftools) #' if that doesn't work --> type: .rs.restartR() in the console and press enter

#' preparing some dummy example data

# input.raw = read.csv('~/Documents/MSKCC/00_Data/Annotated_MAF_oncokb_hugoified.txt', sep = '\t')
# subset_raw = input.raw[sample(x = nrow(input.raw), size = 100, replace = T), ]
# write.table(subset_raw, file = '~/Desktop/dummy_maf.txt', sep = '\t', row.names = F, quote = F)

#' read the inuput data via:
x = maftools::read.maf('~/Desktop/dummy_maf.txt') #' <your_path_to_file>

#' make the visualization
lollipopPlot(maf = x, 
             gene = 'TP53', 
             labelPos = 'all', 
             labPosSize = 0.7, 
             cBioPortal = T)

#' out