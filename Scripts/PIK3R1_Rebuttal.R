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
GOI = c('PTEN', 'PIK3R1', 'PIK3CA')

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



## Question 2: Focality of PIK3R1 alteration:
#' Are these PIK3R1 losses specific or always in the context of wider 5q deletions? 
#' How often are these homozygous vs mono allelic losses? 
#' How often do they co-occur with PTEN/PIK3CA alterations?

gene_call = read.csv('msk_impact_facets_annotated.gene_level.txt.gz', sep = '\t')
gene_call$sample = substr(gene_call$sample, start = 1, stop = 17)
genes_prostate = gene_call[which(gene_call$sample %in% ccf_prostate$Tumor_Sample_Barcode), ]
genes_prostate_goi = genes_prostate[which(genes_prostate$gene %in% GOI), ]  
genes_prostate_goi = genes_prostate_goi[, c('sample', 'gene', 'chrom', 'seg_start', 'seg_end', 'cn_state')]
genes_prostate_goi$length = genes_prostate_goi$seg_end - genes_prostate_goi$seg_start

#' every sample should contain 3 calls (for every gene)
ii = which(table(genes_prostate_goi$sample) != 3)
genes_prostate_goi = genes_prostate_goi[!genes_prostate_goi$sample %in% names(ii), ]


#' Assess the focality based on arm-length
chromo_arm = data.frame(target = c('PIK3CA', 'PIK3R1', 'PTEN'),
                        chromosome = c(3, 5, 10),
                        chromo.arm = c('3q', '5q', '10q'),
                        start = c(90504854, 46405641, 39254935),
                        end = c(198022430, 180915260, 135534747))
chromo_arm$length = chromo_arm$end - chromo_arm$start


gene_out = data.frame()
for(i in 1:nrow(genes_prostate_goi)){
  if(genes_prostate_goi$length[i] >= chromo_arm$length[which(chromo_arm$target == genes_prostate_goi$gene[i])]){
    out = data.frame(target = genes_prostate_goi$gene[i],
                     breadth = paste0('>', chromo_arm$chromo.arm[which(chromo_arm$target == genes_prostate_goi$gene[i])]),
                     call = genes_prostate_goi$cn_state[i],
                     sample = genes_prostate_goi$sample[i])
  } else {
    out = data.frame(target = genes_prostate_goi$gene[i],
                     breadth = paste0(genes_prostate_goi$length[i] / chromo_arm$length[which(chromo_arm$target == genes_prostate_goi$gene[i])], '%'),
                     call = genes_prostate_goi$cn_state[i],
                     sample = genes_prostate_goi$sample[i])
  }
  gene_out = rbind(gene_out, out)
}


## PIK3R1 alteration:
#' Amplification
amp_state = c('AMP', 
              'GAIN', 
              'CNLOH AFTER', 
              'LOSS & GAIN', 
              'LOSS AFTER')

#' Het_category
het_state = c('CNLOH', 
              'CNLOH & GAIN', 
              'CNLOH BEFORE', 
              'CNLOH BEFORE & LOSS', 
              'HETLOSS',
              'LOSS BEFORE',
              'LOSS BEFORE & AFTER')

#' unspecified
unspecified = c('AMP (many states)', 
                'DIPLOID or CNLOH', 
                'GAIN (many states)', 
                'INDETERMINATE', 
                'LOSS (many states)',
                'LOSS BEFORE or DOUBLE LOSS AFTER',
                'TETRAPLOID or CNLOH BEFORE')

homo_state = c('HOMDEL', 'DOUBLE LOSS AFTER')


## PIK3R1:
PIK3R1 = gene_out[which(gene_out$target == 'PIK3R1'), ]
PIK3R1 = PIK3R1[which(PIK3R1$call %in% homo_state | 
                        PIK3R1$call %in% het_state), ]

PIK_vec = c()
for(i in 1:nrow(PIK3R1)){
  if(substr(x = PIK3R1$breadth[i], 
            start = nchar(PIK3R1$breadth[i]), 
            stop = nchar(PIK3R1$breadth[i])) == '%'){
    val = round(as.numeric(substr(PIK3R1$breadth[i], start = 1, stop = 15)), 4)
    
  } else next
  PIK_vec = c(PIK_vec, val)
}

## PIK3CA:
PIK3CA = gene_out[which(gene_out$target == 'PIK3CA'), ]
PIK3CA = PIK3CA[which(PIK3CA$call %in% homo_state | 
                        PIK3CA$call %in% het_state), ]

PIK3CA_vec = c()
for(i in 1:nrow(PIK3CA)){
  if(substr(x = PIK3CA$breadth[i], 
            start = nchar(PIK3CA$breadth[i]), 
            stop = nchar(PIK3CA$breadth[i])) == '%'){
    val = round(as.numeric(substr(PIK3CA$breadth[i], start = 1, stop = 15)), 4)
    
  } else next
  PIK3CA_vec = c(PIK3CA_vec, val)
}


## PTEN:
PTEN = gene_out[which(gene_out$target == 'PTEN'), ]
PTEN = PTEN[which(PTEN$call %in% homo_state | 
                        PTEN$call %in% het_state), ]

PTEN_vec = c()
for(i in 1:nrow(PTEN)){
  if(substr(x = PTEN$breadth[i], 
            start = nchar(PTEN$breadth[i]), 
            stop = nchar(PTEN$breadth[i])) == '%'){
    val = round(as.numeric(substr(PTEN$breadth[i], start = 1, stop = 15)), 4)
    
  } else next
  PTEN_vec = c(PTEN_vec, val)
}


## Are PIK3R1 deletions focal?
PIK_comp = data.frame(target = 'PIK3R1',
                      category = c('<25%', '25-50%', '50-75%', '>75%'),
                      value = c(length(PIK_vec[which(PIK_vec <= 0.25)]),
                                length(PIK_vec[which(PIK_vec > 0.25 & PIK_vec <= 0.5)]),
                                length(PIK_vec[which(PIK_vec > 0.5 & PIK_vec <= 0.75)]),
                                length(PIK_vec[which(PIK_vec > 0.75)])))

#' number > 5q
PIK_comp$value[which(PIK_comp$category == '>75%')] = PIK_comp$value[which(PIK_comp$category == '>75%')] + length(grep('^>.', x = PIK3R1$breadth))
PIK_comp$category = factor(PIK_comp$category, levels = c('>75%', '50-75%', '25-50%', '<25%'))
colours = RColorBrewer::brewer.pal(n = 4, name = 'Blues')

ggplot(PIK_comp, aes(x = target, y = value, fill = category)) +
  geom_bar(stat = 'identity', position = 'fill', width = 0.3) +
  scale_fill_manual(values = colours, name = '% of 5q') +
  coord_flip() +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Fraction', title = 'PIK3R1 losses (segment length) relative to 5q length')


ggplot(as.data.frame(PIK_vec), aes(x = PIK_vec)) +
  geom_density(alpha = 0.2, size = 1.1) +
  labs(x = 'Fraction of 5q arm', y = 'density') +
  scale_y_continuous(expand = c(0, 0.01), limits = c(0, 2)) +
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



## Co-occurrence with PTEN
co = gene_out[which(gene_out$call %in% homo_state | 
             gene_out$call %in% het_state), ]

co_all = data.frame()
for(i in unique(co$sample)){
  da = co[which(co$sample == i), ]
  if(nrow(da) == 1) next
  else if(nrow(da) > 1 & 'PIK3R1' %in% da$target) {
    a_targets = as.character(da$target)
    a_targets = a_targets[!a_targets %in% 'PIK3R1']
    o = data.frame(gene1 = 'PIK3R1',
                   gene2 = a_targets,
                   id = i)
  }
  co_all = rbind(co_all, o)
}

co_all %>%
  group_by(gene1, gene2) %>%
  summarize(count = n())










