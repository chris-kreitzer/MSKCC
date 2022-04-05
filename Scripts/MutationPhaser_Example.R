##-----------------------------------------------------------------------------
## MutationPhaser()
## Test whether two mutations are on the same allele (cis) or opposite (trans)
## 04/05/2022
## chris-kreitzer


clean()

setwd('~/Documents/MSKCC/dmp-2021/MutationPhaser_Test/')
library(MutationPhaser)

#' working on one example: P-0003433-T01-IM5
mutations = read.csv('~/Documents/GitHub/PARP-Prostate-Cancer/Data/data_mutations_extended.oncokb.txt.gz', sep = '\t')
MAF = mutations[which(mutations$Tumor_Sample_Barcode == 'P-0003433-T01-IM5'), ]
MAF = MAF[, c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 
              'End_Position', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2')]
MAF = as.data.table(MAF)
write.table(MAF, 
            file = '~/Documents/MSKCC/dmp-2021/MutationPhaser_Test/gorlick.maf.txt', 
            sep = '\t', quote = F, row.names = F)

#' BAM MAP
BAM_MAP = data.table(Tumor_Sample_Barcode = 'P-0003433-T01-IM5',
                     bam = '~/Documents/MSKCC/dmp-2021/MutationPhaser_Test/YQ584241-T.bam')
write.table(BAM_MAP, file = '~/Documents/MSKCC/dmp-2021/MutationPhaser_Test/gorlick.bam.txt', sep = '\t', row.names = F)


#' phasing
x = MutationPhaser::phase(maf_file = '~/Documents/MSKCC/dmp-2021/MutationPhaser_Test/gorlick.maf.txt', 
                          map_file = '~/Documents/MSKCC/dmp-2021/MutationPhaser_Test/gorlick.bam.txt', 
                          ref_genome = '~/Documents/MSKCC/dmp-2021/Genomics/gr37.fasta',
                          max_reference_region = 1e4,
                          cpus = 4)

annotated = MutationPhaser::annotate_phase(d = x)



