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
x = MutationPhaser::phase(maf_file = '~/Documents/MSKCC/dmp-2021/maf.test.txt', 
                      map_file = '~/Documents/MSKCC/dmp-2021/bam.test.txt', 
                      ref_genome = '~/Documents/MSKCC/dmp-2021/Genomics/gr37.fasta',
                      max_reference_region = 1e4,
                      cpus = 4)





##-----------------------------------------------------------------------------
## Manually extract aligned read at BRCA2 locus

partner = 1
x = MutationPhaser::phase() #' see above
tmp = x[1, ]
tmp$bam = '~/Documents/MSKCC/dmp-2021/BAMs/P-0015132-T01-IM6/JH330573-T.bam'




get_variant_barcodes <- function(tmp,ref_genome,max_reference_region,partner=1,return_barcodes=F) { 
  out <- tryCatch({            
    
    ## get fields for the relevant composite mutation partner
    start_field <- paste('Start_Position',partner,sep='.')
    end_field <- paste('End_Position',partner,sep='.')
    ref_field <- paste('Reference_Allele',partner,sep='.')
    alt_field <- paste('Tumor_Seq_Allele2',partner,sep='.')
    type_field <- paste('Variant_Type',partner,sep='.')
    
    ## bam file for sample
    
    variant_start <- tmp[[start_field]]             
    ref_allele <- tmp[[ref_field]]
    alt_allele <- tmp[[alt_field]]
    variant_type <- tmp[[type_field]]
    
    ## path to bam 
    bamfile <- BamFile(tmp$bam)
    
    ## what param: fields we wamt from the bam
    what <- scanBamWhat()
    
    ## which param: regions for this composite_mutation as GRangesList
    regions <- tmp[, c('Chromosome', start_field,end_field), with=F]
    names(regions) <- c('chr','start','end')
    regions$start <- regions$start
    regions$end <- regions$end
    which <- makeGRangesFromDataFrame(regions)
    
    ## flag param: samtools flags to use when extracting reads
    flags <- scanBamFlag( isNotPassingQualityControls=F, isDuplicate=F)
    
    ## combine into params object; load the bam; create 'dat' with the relevant data
    right_padding <- max(c(nchar(ref_allele),nchar(alt_allele)))
    if(is.na(right_padding)) right_padding <- 0
    what <- c("qname","seq","pos","cigar","qwidth","rname","mrnm","mpos","mate_status") #"rname","qwidth","qual")
    param <- ScanBamParam(which=which, mapqFilter=10, flag=flags, what=what)
    bam <- scanBam(bamfile, param=param)[[1]]        
    region <- paste0(tmp$Chromosome,':',min(bam$pos),'-',max(bam$pos+width(bam$seq)+right_padding))
    bp_spanned_in_region <- abs(max(bam$pos+width(bam$seq)+right_padding) - min(bam$pos))
    
    if(variant_type!='Fusion' & bp_spanned_in_region <= max_reference_region) {
      invisible(capture.output(seq <- get.fasta(x=region,fasta=ref_genome,check.chr=F,verbose=F)$sequence))
      #message('using cdna for p53!'); seq <- readRDS('~/lab/projects/composite_mutations/data/processed_data/figure5/tcga/p53_cdna.fa.rds')
      dat <- data.table(qname=bam$qname,pos=bam$pos,seq=as.character(bam$seq),cigar=bam$cigar,qwidth=bam$qwidth,
                        rname=bam$rname,mpos=bam$mpos,mrnm=as.character(bam$mrnm))
      
      ## get reference sequence for each alignment's SEQ
      fullrefs <- rep(as.character(seq),nrow(dat))
      dat$ref_start <- dat$pos-min(dat$pos)+1
      dat$ref_end <- dat$ref_start + bam$qwidth - 1 + right_padding
      dat$ref <- substr(fullrefs,dat$ref_start,dat$ref_end)                       
      
      ## re-align the sequence of each alignment to the reference sequence for its respective genomic region
      alignments <- recnw(dat$seq,dat$ref,gap_penalty=8,match=5,mismatch=-4,free_hgap_1=T,free_hgap_2=T,free_tgap_1=F,free_tgap_2=F) ## recnw.cpp 
      alignments <- rbindlist(lapply(alignments,function(s) {as.list(strsplit(s,'[|]')[[1]])}))
      names(alignments) <- c('seq1','seq2','score')
    } else if(variant_type!='Fusion' & bp_spanned_in_region > max_reference_region) {
      stop(paste0('Region spanned by SAM alignments exceeds maximum expected (',max_reference_region,'bp)'))
      
    } else if(variant_type=='Fusion'){ 
      ## don't query the reference sequence for fusions, just focus on where the SAM alignments map to
      dat <- data.table(qname=bam$qname,pos=bam$pos,seq=as.character(bam$seq),cigar=bam$cigar,qwidth=bam$qwidth,
                        rname=bam$rname,mpos=bam$mpos,mrnm=as.character(bam$mrnm))
    }
    
    if(variant_type=='DEL') {
      ## Say read has the alternate (deleted) allele if:
      ##  - the alignment returns a deletion at the appropriate locus
      ##  - AND the cigar for the read described a deletion
      ## Say read has the wild-type allele if 
      ##  - the cigar had no D or I
      ## Otherwise:
      ##  - assign allele = 'NA'
      dat$allele <- as.character(NA) ## intialize the allele called by each read as a default 'NA' value
      dat$call <- '.'
      dat$alignment <- alignments$seq1              
      dat$vstart <- variant_start - dat$pos + 1
      dat$vend <- dat$vstart + nchar(ref_allele) - 1     
      dat[vstart < 1,vstart:=1]
      dat[vend < 1,allele:=NA]
      dat[vstart > 0 & vend > 0,allele:=substr(alignment,vstart,vend)]
      dat[grepl('-',allele)==T & grepl('D',cigar)==T,call:=alt_allele] 
      dat[grepl('-',allele)==F,call:=ref_allele] 
      
    } else if(variant_type=='INS') {
      ## Say read has the alternate (inserted) allele if:
      ##  - the alignment returns a insertion at the appropriate locus
      ##  - AND the cigar for the read described a insertion
      ## Say read has the wild-type allele if 
      ##  - the cigar had no D or I
      ## Otherwise:
      ##  - assign allele = 'NA'
      dat$allele <- as.character(NA) ## intialize the allele called by each read as a default 'NA' value
      dat$call <- '.'
      dat$allele <- as.character(NA) ## intialize the allele called by each read as a default 'NA' value
      dat$alignment <- alignments$seq2
      dat$vstart <- variant_start - dat$pos + 2
      dat$vend <- dat$vstart + nchar(alt_allele) - 1
      dat[vstart < 1,vstart:=1]
      dat[vend < 1,allele:=NA]
      dat[vstart > 0 & vend > 0,allele:=substr(alignment,vstart,vend)]
      dat[grepl('-',allele)==T & grepl('I',cigar)==T,call:=alt_allele] 
      dat[grepl('-',allele)==F,call:=ref_allele] 
      
    } else if(variant_type %in% c('SNP','DNP','ONP')) {           
      ## S/D/ONPs, after re-aligning SEQ to the ref, simply check for the ref_allele or 
      ## alt_allele at the expected position
      dat$allele <- as.character(NA) ## intialize the allele called by each read as a default 'NA' value
      dat$call <- '.'
      dat$alignment <- alignments$seq1
      dat$vstart <- variant_start - dat$pos + 1
      dat$vend <- dat$vstart + nchar(alt_allele) - 1
      dat[vstart < 1,vstart:=1]
      dat[vend < 1,allele:=NA]
      dat[vstart > 0 & vend > 0,call:=substr(alignment,vstart,vend)]
      dat[call %nin% c(ref_allele,alt_allele),call:='.'] 
      
    } else if(variant_type %in% 'Fusion') {
      ## experimental!
      
      ## first, check if the fusion includes a gene from a different chromosome from the wild-type allele
      ## here, we make the assumptions:
      ## 1. Either the alt-allele's chromosome is the same as the ref-allele's chromosome,
      ##    or the fusion-allele's chromosome is the second-most common chromosome for SAM reads/mates
      ## 2. Fusions were only called/in MAF if supported by 2+ reads
      ## 3. Reads mapping to the fusion's different chromosome is enough to say they called the fusion (I'm not checking
      ##    that the position at the alternate chromosome is where is it supposed to be for the fusion-partner gene)
      chromosome_table <- table(dat$mrnm)
      chromosome_table <- chromosome_table[order(chromosome_table,decreasing=T)]
      chromosome_table <- chromosome_table[chromosome_table >= 2]
      chromosome_table <- chromosome_table[1:2]
      ref_chromosome <- names(chromosome_table)[1]
      alt_chromosome <- names(chromosome_table)[2]
      if(is.na(alt_chromosome)) alt_chromosome <- ref_chromosome
      
      if(alt_chromosome != ref_chromosome) {
        dat$fusion_pair <- dat$rname %in% alt_chromosome | dat$mrnm %in% alt_chromosome
      } else {
        ## the fusion-partner is on the same chromosome as the referene allele, so now check for reads
        ## with impossibly-large genomic distances from their mate. This will only work for paired-end reads.
        dat[rname %in% ref_chromosome & mrnm %in% ref_chromosome,position_difference:=abs(pos-mpos)]
        dat$distances <- as.integer(NA)
        dat[rname %in% ref_chromosome & mrnm %in% ref_chromosome,distances:=abs(pos-mpos)]
        dat$fusion_pair <- dat$distances >= 2e3 ## hard-coded threshold for maximum short, paired-end read start:start distance    
      }
      ref_allele <- 'no fusion'
      alt_allele <- 'fusion'                
      dat$call <- ifelse(dat$fusion_pair==T,alt_allele,ref_allele)
    } else {
      stop('variant type is not supported!')
    }
    ref_barcodes <- dat$qname[dat$call==ref_allele]
    alt_barcodes <- dat$qname[dat$call==alt_allele]
    list(ref_barcodes=ref_barcodes,alt_barcodes=alt_barcodes,error='',dat=dat)
    
  }, error=function(e) {
    ref_barcodes <- as.character(NA)
    alt_barcodes <- as.character(NA)
    list(ref_barcodes=ref_barcodes,alt_barcodes=alt_barcodes,error=as.character(e),dat='foo')
  })        
  
  if(return_barcodes==F) {
    list(ref_barcodes=out$ref_barcodes,alt_barcodes=out$alt_barcodes,error=out$error)
  } else {
    out$dat  
  }
}

out$ref_barcodes
out$error
View(out$dat)









