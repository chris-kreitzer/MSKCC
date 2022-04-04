##-----------------------------------------------------------------------------
## Mutation Phasing: Investigate whether a mutation occur in cis or trans
## Library from Ed Redznik:
## 
## 04/04/2022
## chris-kreitzer
## 

clean()

## Dependencies
#' install.packages('bedr')
#' install.packages('Hmisc')

# install.packages('~/Documents/MSKCC/MutationPhaser/', repos = NULL, type = 'source')
library(MutationPhaser)



#' working on one example: P-0015132-T01-IM6 
#' look whether BRCA2 Germline and BRCA2 somatic alteration occur on the same allele or not (cis VS trans)  
MAF = mutations_PRAD[which(mutations_PRAD$Tumor_Sample_Barcode == 'P-0015132-T01-IM6'), ]
MAF = MAF[, c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 
              'End_Position', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2')]
MAF = as.data.table(MAF)
write.table(MAF, 
            file = '~/Documents/MSKCC/dmp-2021/maf.test.txt', 
            sep = '\t', quote = F, row.names = F)

#' BAM MAP
BAM_MAP = data.table(Tumor_Sample_Barcode = 'P-0015132-T01-IM6',
                     bam = '~/Documents/MSKCC/dmp-2021/BAMs/P-0015132-T01-IM6/JH330573-T.bam')

write.table(BAM_MAP, file = '~/Documents/MSKCC/dmp-2021/bam.test.txt', sep = '\t', row.names = F)


#' phasing
MutationPhaser::phase(maf_file = '~/Documents/MSKCC/dmp-2021/maf.test.txt', 
                      map_file = '~/Documents/MSKCC/dmp-2021/bam.test.txt', 
                      ref_genome = '~/Documents/MSKCC/dmp-2021/Genomics/gr37.fasta',
                      phased_output = '~/Desktop/tmp/new.txt', 
                      max_reference_region = 1e4,
                      cpus = 4)




