## Investigate AR-deletions in Prostate adenocarcinoma
## Joint-project with Kyrie Pappas
## 
## Start: 07/16/2021


library(pctGCdata)
x = facets::readSnpMatrix('~/Desktop/mnt/ATMcountdata/countsMerged____P-0002273-T01-IM3_P-0002273-N01-IM3.dat.gz')
x1 = facets::preProcSample(rcmat = x, gbuild = 'hg19')
x2 = facets::procSample(x1)
x3 = facets::emcncf(x2)

#' we further have the gene-level matrix where we can directly assess the AR gene (IMPACT data, 40K)
#' AR-loss can promote resistance to PARPi (mouse models)
#' technical options are available
#' however, consider that AR is near the centromere; only one copy (outside of PAR1 and PAR2)

