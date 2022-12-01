##----------------+
## LUAD freqCNA plot;
## for Subhi 
##----------------+

library(copynumber)

luad = read.csv('~/Documents/MSKCC/MSKCC_SlackFiles/LUAD_ascat.seg.txt', sep = '\t')
luad = luad[with(luad, order(ID, chrom, loc.start)), ]
row.names(luad) = NULL

luad$chrom = as.integer(luad$chrom)
colnames(luad) = c('sampleID', 'chrom', 'start.pos', 'end.pos', 'n.probes', 'mean')
luad$start.pos = as.integer(luad$start.pos)
luad$end.pos = as.integer(luad$end.pos)
luad$n.probes = as.numeric(luad$n.probes)


for(i in unique(luad$sampleID)){
  print(i)
  segs = luad[which(luad$sampleID == i), ]
  copynumber::plotFreq(segments = segs, thres.gain = 0.2)
  #' sample: ALCH-ABW9-TTP1-A failed (due to missing chromosome 22)
  #' 
}

luad_new = luad[!luad$sampleID %in% c('ALCH-ABW9-TTP1-A'), ]

#' make the plot
plotFreq(segments = luad_new,
         thres.gain = 0.2,
         thres.loss = -0.2, 
         col.gain = adjustcolor("#b2182b"), 
         col.loss = adjustcolor("darkblue"),
         ylim = c(-75,75))

