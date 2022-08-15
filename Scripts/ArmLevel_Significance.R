## Endometrial Cancer: 

## checking if copy-number alterations (gain or losses) differ among
## primary and metastatic samples (chromosome-arm-wise).
## Subhi request
## 
## start: 08/12/2022
## revision: 08/15/2022
## chris-kreitzer

clean()
gc()
.rs.restartR()
setup()


#' Data-import
hg19 = '~/Profile'
hg19$chrom[which(hg19$chrom == 23)] = 'X'
cohort = read.csv('~/Documents/MSKCC/Subhi/Endometrial_Cancer/CN-H (Mets vs Primary).txt', sep = '\t')
alterations = read.csv('~/Documents/MSKCC/Subhi/Endometrial_Cancer/mskimpact_segments - 2022-08-08T161649.492.seg', sep = '\t')
cytoband = read.csv('~/Documents/MSKCC/05_IMPACT40K/ascets/cytoband_coordinates_hg19.txt', sep = '\t')
# write.table(cytoband, file = '~/Documents/MSKCC/dmp-2021/Genomics/cytoband_coordiantes_hg19.txt', sep = '\t', row.names = F, quote = F)


#' assign whether called fragment belong to q or p arm of chromosome
#' I am using ASCETS for this task:
#' Importantly, we can use more granular cytobands if we want to:
#' In this approach, we are just splitting the chromosomes in p and q arm!

chromo_arm = data.frame()
for(i in 1:nrow(hg19)){
  chrom = hg19$chrom[i]
  out = data.frame(arm = c(paste0(chrom, 'p'), paste0(chrom, 'q')),
                   chrom = chrom,
                   start = c(0, hg19$centromere[i]),
                   end = c(hg19$centromere[i], hg19$size[i]))
  
  chromo_arm = rbind(chromo_arm, out)
  
}

chromo_arm$arm[which(chromo_arm$chrom == 24)] = c('Yp', 'Yq')
chromo_arm$chrom[which(chromo_arm$chrom == 24)] = c('Y', 'Y')


## RUNNING ASCETS -
source('~/Documents/MSKCC/05_IMPACT40K/ascets/ascets_resources.R')
colnames(alterations) = c('sample', 'chrom', 'segment_start', 'segment_end', 'num_mark', 'log2ratio')

#' run algorithm
ascets_output = ascets(cna = alterations, 
                       cytoband = chromo_arm, 
                       name = 'Endometrial_Cancer', 
                       min_boc = 0.7)
write_outputs_to_file(ascets_output, location = "~/Documents/MSKCC/Subhi/Endometrial_Cancer/")


##' Assign significance - is there a difference between primary and mets
##' stratified by cytoband
all_out = data.frame()
for(i in 2:ncol(ascets_output$calls)){
  try({
    #' primaries
    primary_gain = table(ascets_output$calls[which(ascets_output$calls$sample %in% cohort$Sample_ID[which(cohort$Sample.Type == 'Primary')]), i])['AMP'][[1]]
    primary_gain = ifelse(is.na(primary_gain) | primary_gain == 0, 0, primary_gain)
    primary_loss = table(ascets_output$calls[which(ascets_output$calls$sample %in% cohort$Sample_ID[which(cohort$Sample.Type == 'Primary')]), i])['DEL'][[1]]
    primary_loss = ifelse(is.na(primary_loss) | primary_loss == 0, 0, primary_loss)
    primary_neutral = table(ascets_output$calls[which(ascets_output$calls$sample %in% cohort$Sample_ID[which(cohort$Sample.Type == 'Primary')]), i])['NEUTRAL'][[1]]
    #' mets
    mets_gain = table(ascets_output$calls[which(ascets_output$calls$sample %in% cohort$Sample_ID[which(cohort$Sample.Type == 'Metastasis')]), i])['AMP'][[1]]
    mets_gain = ifelse(is.na(mets_gain) | mets_gain == 0, 0, mets_gain)
    mets_loss = table(ascets_output$calls[which(ascets_output$calls$sample %in% cohort$Sample_ID[which(cohort$Sample.Type == 'Metastasis')]), i])['DEL'][[1]]
    mets_loss = ifelse(is.na(mets_loss) | mets_loss == 0, 0, mets_loss)
    mets_neutral = table(ascets_output$calls[which(ascets_output$calls$sample %in% cohort$Sample_ID[which(cohort$Sample.Type == 'Metastasis')]), i])['NEUTRAL'][[1]]
    
    #' chi-squared test
    print(c(i, matrix(c(primary_gain, mets_gain, primary_neutral, mets_neutral), ncol = 2)))
    chisq_gain = chisq.test(matrix(c(primary_gain, mets_gain, primary_neutral, mets_neutral), ncol = 2), correct = F)$p.value
    chisq_loss = chisq.test(matrix(c(primary_loss, mets_loss, primary_neutral, mets_neutral), ncol = 2), correct = F)$p.value
    
    out = data.frame(arm = colnames(ascets_output$calls)[i],
                     gain = chisq_gain,
                     loss = chisq_loss)
    
    all_out = rbind(all_out, out)
  })
}



#' assign_arm = function(seg.file, genomic_coordinates){
#'   all_out = data.frame()
#'   try({
#'     for(i in unique(seg.file$chrom)){
#'       centromere = genomic_coordinates$centromere[which(genomic_coordinates$chrom == i)]
#'       seg_file = seg.file[which(seg.file$chrom == i), ]
#'       seg_file$arm = NA
#'       #' assign arm
#'       seg_file$arm = ifelse(seg_file[, 'loc.start'] <= centromere &
#'                               seg_file[, 'loc.end'] <= centromere, 'p',
#'                             ifelse(seg_file[, 'loc.start'] >= centromere &
#'                                      seg_file[, 'loc.end'] >= centromere, 'q', 'centromere'))
#'       
#'       all_out = rbind(all_out, seg_file)
#'     }
#'     
#'     
#'     #' chose one alteration/sample
#'     #' if there are more than 5 segments per chromosome-arm; 
#'     
#'     
#'   })
#'   
#'   return(all_out)
#' }
#' 
#' annotation = assign_arm(seg.file = alterations, genomic_coordinates = hg19)
