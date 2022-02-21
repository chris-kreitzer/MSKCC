## Querying AR in Prostate cancer cohort
## Specifically asking for deletions (and deep deletions)
## 02-21-2022

x = read.csv('~/Documents/MSKCC/dmp-2021/dmp_022022/data_mutations_extended.txt', sep = '\t', skip = 1)
cohort = read.csv('~/Documents/MSKCC/dmp-2021/dmp_022022/ProstateCohort_022022.tsv', sep = '\t')
cna = read.csv('~/Documents/MSKCC/dmp-2021/dmp_022022/data_CNA.txt', sep = '\t')
Prostate = x[which(x$Tumor_Sample_Barcode %in% cohort$Sample.ID),, drop = F]


AR = cna[which(cna$Hugo_Symbol == 'AR'), ]
colnames(AR) = gsub(pattern = '\\.', replacement = '-', colnames(AR))
AR = AR[,which(colnames(AR) %in% Prostate$Tumor_Sample_Barcode)]

colnames(AR)[which(AR == -2, arr.ind = TRUE)]






# Two samples show a AR deep deletion:
# - 'P-0003101-T02-IM5' (metastasis in Lung 2015 and Kidney (20% purity) 2016: called from cBIOa)
# - 'P-0004482-T01-IM5' (metastasis in Liver (60% purity) 2015 added)
#' both are FacetsQC TRUTH
#' both are non-WGD

msk_annotated = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/msk_impact_facets_annotated.cohort.txt', sep = '\t')
msk_annotated$sample_id = substr(msk_annotated$sample_id, start = 1, stop = 17)
msk_annotated = msk_annotated[which(msk_annotated$sample_id %in% c('P-0003101-T02-IM5', 'P-0004482-T01-IM5')), ]
msk_gene = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/msk_impact_facets_annotated.gene_level.txt', sep = '\t')
msk_gene$sample = substr(msk_gene$sample, start = 1, stop = 17)
msk_gene = msk_gene[which(msk_gene$sample %in% c('P-0003101-T02-IM5', 'P-0004482-T01-IM5')), ]
msk_gene = msk_gene[which(msk_gene$gene == 'AR'), ]
msk_gene
