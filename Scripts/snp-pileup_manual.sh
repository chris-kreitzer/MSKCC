#! /bin/bash/

## First, `snp-pileup` was used to get count matrices for respective .bam files (see excel sheet for sample paths)

#' example run:
I used the following vcf file (see below); stored locally
Parameters used for snp-pileup based on recommondation of Ahmet Zehir

`snp-pileup -v -A -q15 -Q20 -r10 -g /juno/home/kreitzec/Data/dbsnp_137.hg19__RmDupsClean__plusPseudo50__DROP_SORT.vcf.gz <output-folder> <normal.bam> <tumor.bam>`

#' I simply submitted this job via 
`bash <filename>`; since bsub was not working in batch mode


Running snp-pileup function (08/30/2022): inspired by Reliable Analysis of Clinical Tumor-Only Whole-Exome Sequencing Data; Markus Riester)

snp-pileup -A -P50 -q15 -Q20 -r10,0 /juno/work/ci/resources/genomes/GRCh37/facets_snps/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf.gz <<output_name>> <<NORMAL_BAM>> <<TUMOR_BAM>>

A: count_orphans
P50: count read-depth every 50 bp postion (even if there is no SNP location from dbsnp
q15: minimum mapping Quality
Q20: minimum base quality
r10,0: there needs to be at least 10 reads at every position in the first file provided (in this case the normal); so other positions that have less than 10 reads are skipped

