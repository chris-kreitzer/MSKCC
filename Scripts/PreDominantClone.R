PRAD = read.csv('~/Documents/MSKCC/PIK3R1/PIK3R1_Facets_Calls', sep = '\t')
sampleIDs = PRAD$SAMPLE[which(PRAD$FACETS_QC == 'TRUE')]

mutations = read.csv('~/Documents/MSKCC/PIK3R1/Data/msk_impact_facets_annotated.ccf.maf.gz', sep = '\t')
mutations = mutations[which(mutations$Tumor_Sample_Barcode %in% sampleIDs), ]
View(mutations)

mut_high = mutations[which(mutations$Tumor_Sample_Barcode == 'P-0057414-T01-IM6'), ]
plot(density(mut_high$ccf_Mcopies))
plot(density(mut_high$t_var_freq))
hist(mut_high$ccf_Mcopies, nclass = 20, right = T)
lines(density(mut_high$ccf_Mcopies))


#' example ccf plot:
par(mfrow = c(2,2),
    mar = c(4,2,3,2))

#' TP53 mutations
hist(mutations$ccf_Mcopies[which(mutations$Hugo_Symbol == 'TP53')], 
     nclass = 50, col = 'black', border = 'white',
     las = 2,
     xaxt = 'n', 
     main = NULL,
     xlab = '',
     yaxt = 'n',
     ylab = '')
axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1),
     labels = c(0, 0.25, 0.5, 0.75, 1))
box(lwd = 2)
mtext(side = 2, text = "Frequency", line = 0.5, font = 2, adj = 0)
mtext(side = 1, text = "CCF", line = 2, font = 2, adj = 0)
title(main = 'TP53 mutations')

#' RB1
hist(mutations$ccf_Mcopies[which(mutations$Hugo_Symbol == 'RB1')], 
     nclass = 50, col = 'black', border = 'white',
     las = 2,
     xaxt = 'n', 
     main = NULL,
     xlab = '',
     yaxt = 'n',
     ylab = '')
axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1),
     labels = c(0, 0.25, 0.5, 0.75, 1))
box(lwd = 2)
mtext(side = 2, text = "Frequency", line = 0.5, font = 2, adj = 0)
mtext(side = 1, text = "CCF", line = 2, font = 2, adj = 0)
title(main = 'RB1 mutations')

#' PIK3R1
hist(mutations$ccf_Mcopies[which(mutations$Hugo_Symbol == 'PIK3CA')], 
     nclass = 50, col = 'black', border = 'white',
     las = 2,
     xaxt = 'n', 
     main = NULL,
     xlab = '',
     yaxt = 'n',
     ylab = '')
axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1),
     labels = c(0, 0.25, 0.5, 0.75, 1))
box(lwd = 2)
mtext(side = 2, text = "Frequency", line = 0.5, font = 2, adj = 0)
mtext(side = 1, text = "CCF", line = 2, font = 2, adj = 0)
title(main = 'PIK3CA mutations')

#' PIK3R1
#' PIK3R1
hist(mutations$ccf_Mcopies[which(mutations$Hugo_Symbol == 'PIK3R1')], 
     nclass = 50, col = 'black', border = 'white',
     las = 2,
     xaxt = 'n', 
     main = NULL,
     xlab = '',
     yaxt = 'n',
     ylab = '')
axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1),
     labels = c(0, 0.25, 0.5, 0.75, 1))
box(lwd = 2)
mtext(side = 2, text = "Frequency", line = 0.5, font = 2, adj = 0)
mtext(side = 1, text = "CCF", line = 2, font = 2, adj = 0)
title(main = 'PIK3R1 mutations')


##-----------------------------------------------------------------------------
mutations = read.csv('~/Documents/MSKCC/datahub/public/prad_su2c_2019/data_mutations.txt', sep = '\t')
mutations$VAF = mutations$t_alt_count / (mutations$t_alt_count + mutations$t_ref_count)

#' print VAF values
hist(mutations$VAF[which(mutations$Hugo_Symbol == 'TP53')], 
     nclass = 50, col = 'black', border = 'white',
     las = 2,
     xaxt = 'n', 
     main = NULL,
     xlab = '',
     yaxt = 'n',
     ylab = '')
axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1),
     labels = c(0, 0.25, 0.5, 0.75, 1))
box(lwd = 2)
mtext(side = 2, text = "Frequency", line = 0.5, font = 2, adj = 0)
mtext(side = 1, text = "CCF", line = 2, font = 2, adj = 0)
title(main = 'PIK3R1 mutations')


