## PIK3R1 Re-Analysis and addressing comments of 
## Reviewers: PIK3R1 Rebuttal
## 
## 01/24/2022
## chris-kreitzer
## 

clean()
setwd('~/Documents/MSKCC/PIK3R1/')
library(ggplot2)

genomics = load('Data/1.0.genomic_data_all.Rdata')
clinical = load('Data/1.1.data_samples.Rdata')
facets = read.csv('Data/facets_qc_TRUE.txt', sep = '\t')
facets$sample_id = substr(facets$sample_id, start = 1, stop = 17)
samples = read.csv('Data/60bfa9dfe4b0f0abac8b43f5_clinical_data.tsv', sep = '\t')
ccf = read.csv('Data/msk_impact_facets_annotated.ccf.maf.gz', sep = '\t')


## PIK3R1 clonality status:
samplesQC = samples[which(samples$Sample.ID %in% facets$sample_id), ]
ccf_prostate = ccf[which(ccf$Tumor_Sample_Barcode %in% samplesQC$Sample.ID), ]

PIK = ccf_prostate[which(ccf_prostate$Hugo_Symbol == 'PIK3R1'), ]
mets = unique(samples$Sample.ID[which(samples$Sample.Type == 'Metastasis')])

#' Sample-Type
PIK$type = NA
for(i in 1:nrow(PIK)){
  if(PIK$Tumor_Sample_Barcode[i] %in% mets){
    PIK$type[i] = 'Metastasis'
  } else
    PIK$type[i] = 'Primary'
}


boxplot(PIK$ccf_expected_copies ~ PIK$type, xlab = '', ylab = 'cancer cell fraction', 
        main = 'CCF of PIK3R1 across different sampling sites')
t.test(PIK$ccf_expected_copies ~ PIK$type)


#' Display the VAF evolution in metastatic samples
ggplot(PIK, aes(x = t_var_freq, color = type, fill = type)) +
  geom_density(alpha = 0.2, size = 1.1) +
  scale_color_manual(values = c('Primary' = '#3b685f',
                                'Metastasis' = '#85510d'), name = '') +
  scale_fill_manual(values = c('Primary' = '#3b685f',
                                'Metastasis' = '#85510d'), name = '') +
  labs(x = 'Variant allele frequency', y = 'density') +
  scale_y_continuous(expand = c(0, 0.01), limits = c(0, 6)) +
  scale_x_continuous(expand = c(0, 0.01), limits = c(0, 1)) +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12, face = 'bold', color = 'black'),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        legend.position = c(0.85, 0.8))

#' Are PIK3R1 clonality estimates independent of PTEN and PIK3CA
GOI = c('PTEN', 'PIK3R1', 'PIK3CA')

ccf_prostate_GOI = ccf_prostate[which(ccf_prostate$Hugo_Symbol %in% GOI), ]
co_table = as.matrix(t(xtabs(~ccf_prostate_GOI$Hugo_Symbol + ccf_prostate_GOI$clonality,
                             subset = ccf_prostate_GOI$clonality != 'INDETERMINATE')))

#' fisher's test
fisher.test(co_table, simulate.p.value = T, B = 1e5)


## BradleyTerry Model for ranking:





length(unique(samples$Sample.ID))
head(ccf$Tumor_Sample_Barcode)
head(samples$Sample.ID)
