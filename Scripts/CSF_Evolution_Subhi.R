###############################################################################
## Copy-Number Evolution in CSF samples from primary Tumors
## 
## 06/20/2022
## chris-kreitzer
## 

clean()
gc()
.rs.restartR()

setwd('~/Documents/MSKCC/Subhi/CSF/')

PatientID_a = readxl::read_excel(path = 'CSF_Cohort_Evolution.xlsx', sheet = 1)
PatientID_b = readxl::read_excel(path = 'CSF_Cohort_Evolution.xlsx', sheet = 2)
FacetsAnnotation = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')

folders = list.files(path = getwd(), pattern = '^0.*', full.names = T)
all_files = c()
for(i in 1:length(folders)){
  files = list.files(path = folders[i], full.names = T)
  all_files = c(all_files, files)
}


## IM-GBM-204:
all_files[grep('s_C_K6YHWN_S001_d', all_files)]


#' load all count-matrices simultanously
#' working on one sample as an example!
input = lapply(list.files(path = 'IM-GBM-204/', full.names = T), function(x) facets::readSnpMatrix(filename = x))

comp = data.frame()
for(countmatrix in 1:length(input)){
  facets_processed = facetsSuite::run_facets(read_counts = input[[countmatrix]],
                                             cval = 100,
                                             genome = 'hg19')
  
  QC = facets_fit_qc(facets_processed)$facets_qc
  QC_waterfall = facets_fit_qc(facets_processed)$waterfall_filter_pass
  dipLogR = facets_processed$dipLogR
  
  for(gene in c('CDKN2A',
                'CDK4',
                'CDK6',
                'PTEN',
                'EGFR',
                'PDGFRA',
                'KIT',
                'KDR',
                'MET',
                'MDM2',
                'MDM4',
                'RB1',
                'NF1',
                'TP53')){
    
    snps = extract_genes_rd(gene = gene,
                            countmatrix = facets_processed$snps,
                            type = 'snps',
                            padding = 0)
    median_cnlr = median(snps$cnlr)
    nhet = length(snps$het[which(snps$het == 1)])
    gene_change = facetsSuite::gene_level_changes(facets_output = facets_processed,
                                                  genome = 'hg19')
    gene_cn = gene_change[which(gene_change$gene == gene), 'cn_state']
    alteration = ifelse(abs(dipLogR - median_cnlr) >= 0.4 &
                          median_cnlr >= 0.2, 'AMP',
                        ifelse(abs(dipLogR - median_cnlr) >= 0.4 &
                                 median_cnlr <= -0.2, 'DEL', 'Diploid'))
    
    out = data.frame(Sample = 'IM-GBM-204',
                     order = countmatrix,
                     gene = gene,
                     gene_cn = gene_cn,
                     alteration = alteration,
                     median_cnlr = median_cnlr,
                     nhet = nhet,
                     QC = QC,
                     waterfall_filter_pass = QC_waterfall,
                     dipLogR = dipLogR)
    comp = rbind(comp, out)
  }
}

ggplot(comp, aes(x = order, y = gene, fill = median_cnlr, color = alteration)) +
  geom_tile(aes(color = alteration, width = 0.8, height = 0.8), size = 0.45) +
  scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red', name = 'CnLR') +
  scale_color_manual(values = c('AMP' = 'red',
                                'DEL' = 'blue',
                                'DIPLOID' = 'grey95'), name = '') +
  scale_x_continuous(breaks = 1:3,
                     labels = c("primary", "CSF1", 'CSF2'),
                     sec.axis = dup_axis()) +
  theme_bw() +
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0),
        axis.text.x = element_blank(),
        axis.text = element_text(color = 'black')) +
  labs(x = '', y = '', title = 'IM-GBM-204') +
  
  coord_equal()


