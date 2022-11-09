segs = read.csv('~/Documents/MSKCC/Subhi/B2M/mskimpact_segments - 2022-11-08T112923.114.seg', sep = '\t')
segs = segs[which(segs$chrom == 15), ]

breakpoints = c()
for(i in unique(segs$ID)){
  n = length(segs$seg.mean[which(segs$ID == i)])
  m = n-1
  if(m == 0) next
  loc.end = segs$loc.end[which(segs$ID == i)][seq(1,m,1)]
  breakpoints = c(breakpoints, loc.end)
}


plot(density(breakpoints), 
     ylab = '',
     xlab = '',
     xaxt = 'n',
     yaxt = 'n',
     main = '', lwd = 2)
axis(side = 1, at = c(2e+07, 4e+07, 6e+07, 8e+07, 1e+08), labels = c('20Mb', '40Mb', '60Mb', '80Mb', '100Mb'))
mtext(text = 'Genomic Position', side = 1, line = 1.8)
mtext(text = 'Density', side = 2, line = 0.8)
mtext(text = 'Breakpoint Density; chr. 15', side = 3, line = 0.8, adj = 0, cex = 1.5)
text(x = 4.9e+07, y = 0.2e-08, label = '44.2Mb\n(FRMD5)')
box(lwd = 2)

x = density(breakpoints)
x$x[which.max(x$y)]

abline(v = x$x[which.max(x$y)], lty = 'dashed', lwd = 0.5)
