## PIK3R1 Re-Analysis and addressing comments of 
## Reviewers: PIK3R1 Rebuttal
## 
## 01/24/2022
## chris-kreitzer
## 

clean()
setwd('~/Documents/MSKCC/PIK3R1/')
library(ggplot2)
library(BradleyTerry2)


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
GOI = c('PTEN', 'PIK3R1', 'PIK3CA', 'MTOR')

ccf_prostate_GOI = ccf_prostate[which(ccf_prostate$Hugo_Symbol %in% GOI), ]
co_table = as.matrix(t(xtabs(~ccf_prostate_GOI$Hugo_Symbol + ccf_prostate_GOI$clonality,
                             subset = ccf_prostate_GOI$clonality != 'INDETERMINATE')))

#' fisher's test
fisher.test(co_table, simulate.p.value = T, B = 1e5)


## BradleyTerry Model for ranking:
test = ccf_prostate_GOI[, c('Hugo_Symbol', 'Tumor_Sample_Barcode', 'ccf_expected_copies', 'clonality')]

a = subset(test, duplicated(test$Tumor_Sample_Barcode))
a = unique(a$Tumor_Sample_Barcode)
test = test[which(test$Tumor_Sample_Barcode %in% a), ]

allout = data.frame()
for(i in unique(test$Tumor_Sample_Barcode)){
  dat = test[which(test$Tumor_Sample_Barcode == i), ]
  if(nrow(dat) == 2 & !any(is.na(dat$ccf_expected_copies))){
    out = data.frame(winner = dat$Hugo_Symbol[which.max(dat$ccf_expected_copies)],
                     loser = dat$Hugo_Symbol[which.min(dat$ccf_expected_copies)],
                     ID = i)
    
    allout = rbind(allout, out)
  } else next
  
}
  
#' manually create matrix based on allout
k = matrix(c(0, 7, 6, 3, 0, 0,  1, 0, 0), 
           nrow = 3, byrow = T,
           dimnames = list(c('PTEN', 'PIK3CA', 'PIK3R1'),
                           c('PTEN', 'PIK3CA', 'PIK3R1')))

m = BradleyTerry2::countsToBinomial(k)

PIK_model = BTm(cbind(win1, win2), player1, player2, ~player, id = "player", data = m)

#' Visualization of the BT model outpt
plot.in = data.frame(item = c('PTEN', 'PIK3CA', 'PIK3R1'),
                            estimate = c(0, -0.8473, -1.7918),
                            SE = c(0, 0.6901, 1.08))

helper.lines.end = length(unique(plot.in$item)) + 1
plot.in$category = ifelse(plot.in$estimate > 0, 'clonal', 'subclonal')

cut_line = plot.in[plot.in$estimate >=0, ]
cut_line = cut_line$item[which.min(cut_line$estimate)]
cut_point = which(plot.in$item == cut_line)
cut_point = (length(unique(plot.in$item)) - cut_point) + 0.5

cancer_type = 'Prostate Cancer'
cancer_type_sample_size = length(unique(test$Tumor_Sample_Barcode))


BT.MLE.model = ggplot(plot.in,
                      aes(x = reorder(item, estimate), 
                          y = estimate, 
                          ymin = estimate - SE, 
                          ymax = estimate + SE)) +
  
  geom_pointrange(size = 0.5) +
  
  geom_linerange() +
  
  coord_flip() +
  
  geom_vline(xintercept = seq(0.5, helper.lines.end, 1), 
             size = 0.2, 
             linetype = 'dashed', 
             color = 'grey') +
  
  scale_y_continuous(expand = c(0, 0.01),
                     limits = c(-3, 1),
                     sec.axis = dup_axis()) +
  
  theme(panel.background = element_blank(),
        axis.text.y = element_text(size = 12, 
                                   face = 'bold',
                                   color = 'black'),
        axis.ticks.y = element_line(lineend = 'round', color = 'black', size = 0.6),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(size = 0.6, lineend = 'round'),
        panel.border = element_rect(fill = NA, color = 'black', size = 1.4)) +
  
  labs(x = '', y = '', title = paste0(cancer_type, ' (n=', cancer_type_sample_size, ')'))

BT.MLE.model

