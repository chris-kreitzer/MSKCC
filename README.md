## Miscellaneous scripts, processing pipelines and analysis steps for colleagues

#### annotate-maf-wrapper.R
example script, where annotate-maf-wrapper.R is located in `${HOME}`   
* `Rscript <script_location> -m <maf_location> -f <facets_output.Rdata_location> -o <output_directory>`   
* `Rscript /juno/home/kreitzec/annotate-maf-wrapper.R -m /juno/work/schultz/subhi/mafs_for_mafAnno/F12-6-p1--F12-blood.union.v4.annotated.maf -f /juno/work/subhi/Ingo_WES_Facets/s_F12-6-p1.gz.Rdata -o /juno/work/schultz/subhi/F12-6-p1_annotated.txt`   



## Data access various formats and important links (MSK related); VPN required    
[MSK-knowledgeSystems](https://github.mskcc.org/)  access only via VPN   
[Delphi-Sample-Tracker](https://delphi.mskcc.org/sample-tracker/home)   
[High Performance Computing](http://mskcchpc.org/) Questions related to juno and server access   
[CVR Molecular Diagnostic service](https://cvr.mskcc.org:8083/search) - open @Chrome: mskcc Christoph    
[DMP-2022 DMP Data Directory](https://github.mskcc.org/knowledgesystems/dmp-2022) #' Page 404 not found / need to login @Enterprise (msk initials)    
[OncoKB annotated Mutations](https://github.mskcc.org/knowledgesystems/oncokb-annotated-msk-impact) daily updates with OncoKB annotations / need login as well    
[Darwin](https://ddp.mskcc.org/search) Clinical Interpretation tool; all clinical data associated with tumor specimen; specifically pathology and radiology events and treatment data.  
[TCGA_SCNA] (/juno/work/ccs/tcga_facets/)   


## Survival analysis; be cautions about left-truncation and length bias
Read this handy [article](https://towardsdatascience.com/how-well-do-you-really-know-survival-analysis-221885b89a8e#:~:text=Left%20truncation%20occurs%20when%20data,study)    
Use this [Github](https://github.com/slb2240/delayed_entry_clin_genom_studies/blob/main/crc_stage_iv_os_dx.R) for R   
And some other useful [publication](https://pubmed.ncbi.nlm.nih.gov/34734967/)  


## BAM(s) intercourse; sequencing reads   
`xjuno`  
`module load samtools/1.9`   

### samtools idxstats – reports alignment summary statistics
- Retrieve and print stats in the index file corresponding to the input file.   
- Before calling idxstats, the input BAM file should be indexed by samtools index.   
- OUTPUt to stdout: reference sequence name, sequence length, # mapped read-segments and # unmapped read-segments

#### example:
`samtools idxstats <sample.bam>`


#### samtools view -h <sample.bam> chr20 > <out>.sam/bam
- With no options or regions specified, prints all alignments in the specified input alignment file in SAM format  
  
`chr1`   
Output all alignments mapped to the reference sequence named `chr1' (i.e. @SQ SN:chr1).

`chr2:1000000`   
The region on chr2 beginning at base position 1,000,000 and ending at the end of the chromosome.


## BigWig (bedGraph) format:
- The bedGraph (or bigwig) format is always the same: `chr-start-end-value`.  
- the score can be anything, e.g. an average read coverage
  - typical file extension: .bw, .bigwig
  - binary version of a bedGraph or wig file   


### snp-pileup -- countmatrix(facets) -- IGV
Importantly. If we are running `snp-pileup -v -A -g <vcf> <normal> <tumor.bam>` with no specified filters
(e.g. read quality etc.) we exactly get the raw alignment counts as displayed in IGV (see P-0014995 bam / IGV / countmatrix_unfiltered)
Moreover, we only get read count information on those positions which are provided in the vcf file supplied (so not every base postion is covered from the IGV). Facets uses the VCF file `/juno/work/ci/resources/genomes/GRCh37/facets_snps/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf` to assess the counts at specfic loci.   
![image](https://user-images.githubusercontent.com/58993468/166894021-8f51a751-c07b-48eb-8639-51be064c184c.png)

  


### BAM/SAM INFO important!
- **FLAG** = bit-wise Flag: Infos about duplicated read or NOT
- Example: Make sure you're only looking at alignment records that represent mapped reads. The -F 0x4 option says to filter records where the 0x4 flag (read unmapped) is 0, resulting it only mapped reads being output.
- **MAPQ** = mapping quality of reads; depends on alignment tool; generally MAPQ >> 30 == UNIQUE mapping (good)
- **CIGAR** = exact describtion of the alignment (whether read matches/mismatches - insertion (I) - deletion (D) or (S) == soft clipping)
- Example: Analyzing the CIGAR string for indels

### Read information: Alignment (mostly bowtie2)
- If the mapping quality score (MAPQ) of a read is zero, then it means that the read maps to multiple positions equally well. 
- MAPQ is a property of single reads, and is not changed by pairing. 
- MAPQ = mapping quality = uniquness of read mapping to genome; MAPQ >> 30 unique in the genome; MAPQ << 15 random match
  (−10 log10 Pr{mapping position is wrong})
- **Cigar 100M**: 
  - M match/mismatch
  - N indicates spliced. 
  - I is insertion into the reference and 
  - D is deletion from the reference.  

e.g. cigar of 100M: indicates that the read starts at 1002 and is mapped continuously for the next 100 bases. 
Note that still there might be mismatches. Suppose your cigar was 20M 1D 20M 60N 30M 1I 20M 100N 9M => read length = sum of all M/I/S operations. 
There is no S here. So, read length = 20+20+30+1+20+9=100. Suppose your read started at 1002, then, let's see how the read is for this cigar string.

![CIGAR](https://user-images.githubusercontent.com/58993468/154679668-5ad990ec-36e6-43a1-a9bf-f403be00ecfe.png)

#### Examples (https://wikis.utexas.edu/display/CoreNGSTools/Filtering+with+SAMTools)
Count reads that mappedi with indels   
`samtools view -F 0x4 yeast_pe.sort.bam | cut -f 6 | grep -P '[ID]' | wc -l`

#### Count all mapped reads
`samtools view -c -F 0x4 yeast_pe.sort.bam`

#### Filtering by location range
##### count the number of reads mapped to chromosome 2 (chrII)
`samtools view -c -F 0x4 yeast_pe.sort.bam chrII`

##### count the number of reads mapped to chromosomes 1 that overlap coordinates 1000-2000
`samtools view -c -F 0x4 yeast_pe.sort.bam chrI:1000-2000`   
 
#### counting only mapped (primary aligned) reads
`samtools view -c -F 260 SAMPLE.bam` [-F 260 (excludes unmapped and secondary alignments)]      
 
options   
  -c  count reads and print the total number

  -f bitcode  output reads that fulfill the checked 'bitcode' criteria, see SAM bitcode fields

  -F bitcode  exclude reads that match one or more checked 'bitcode' criteria, see SAM bitcode fields

  -F 260  output primary aligned mapped reads (no secondary alignments included)
  
## Sequencing coverage and breadth of coverage   
Sequencing coverage or depth (coverage and depth are used interchangeably) determines the number of times sequenced nucleotide bases covered the target genome. For example, if genome size is 100 Mbp and you have sequenced 5 M reads of 100 bp size, then sequencing coverage at genome level would be 5X. 

  
### Median insert size
samtools stats DE840153-T.bam |grep ^SN | cut -f 2-

#### get read depth for each position on chromosome [use -a parameter to get read depth for all positions]  
samtools depth PC14_L001_R1.bam > read_depth.txt

#### get overall read depth
awk '{sum+=$3;} END {print sum/NR;}' read_depth.txt
- $3 means read depth at each position of chromosome (third column from read_depth.txt)
- NR means total rows in read_depth.txt [total mapped chromosome size (with -a option, you will get whole chromosome size)]   
  
#### IMPORTANT Bioinformatic TOOLs:
MutationPhaser (https://github.com/reznik-lab/MutationPhaser/blob/master/R/func.R)   
ASCAT.sc (https://github.com/VanLoo-lab/ASCAT.sc)   
TcellExtract (https://github.com/McGranahanLab/TcellExTRECT)   
refphase (ASCAT tutorial): mirrored subclonal imbalance: https://bitbucket.org/schwarzlab/refphase/src/master/   
GeneticAnnotations: https://atlasgeneticsoncology.org/teaching/30067/nomenclature-for-the-description-of-mutations-and-other-sequence-variations   
Genome Nexus (Ino): https://www.genomenexus.org/    
  
PingID (enrollment): new device:  https://thespot.mskcc.org/esc/?id=kb_article&sys_id=bc010b9d1b0bf0d006ec0e94604bcb2f&table=kb_knowledge
You must be connected via VPN

Test example to have green button! not allowed

  
  

