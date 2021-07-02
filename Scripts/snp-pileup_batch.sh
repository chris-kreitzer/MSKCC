#!/bin/bash

#BSUB -J WES
#BSUB -n 4
#BSUB -R rusage[mem=4]
#BSUB -W 01:59
#BSUB -o s.log
#BSUB -eo s.err

module load samtools/1.9

snp-pileup -v -A -q15 -Q20 -r10 -g /juno/home/kreitzec/Data/dbsnp_137.hg19__RmDupsClean__plusPseudo50__DROP_SORT.vcf.gz /juno/home/kreitzec/test2.gz /juno/home/kreitzec/normal.bam /juno/home/kreitzec/tumor.bam

