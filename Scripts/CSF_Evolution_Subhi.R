###############################################################################
## Copy-Number Evolution in CSF samples from primary Tumors
## 
## 06/20/2022
## chris-kreitzer


clean()
gc()
.rs.restartR()

setwd('~/Documents/MSKCC/Subhi/CSF/')


PatientCohort = readxl::read_excel(path = 'Data/Chris_Data_Facets.xlsx', sheet = 1)
FacetsAnnotation = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')

#' Strive through all CSF count-matrices from juno cluster
CSF_folders = list.files(path = getwd(), pattern = '^0.*', full.names = T)
CSF_files = c()
for(i in 1:length(CSF_folders)){
  files = list.files(path = CSF_folders[i], full.names = T)
  CSF_files = c(CSF_files, files)
}

IMPACT_files = list.files(path = paste0(getwd(), '/CSF_IMPACT_COUNTMATRIX/'),  full.names = T)

#' merge all file paths:
all_files = c(CSF_files, IMPACT_files)


#' add path variable to table
PatientCohort$path = NA
for(i in 1:nrow(PatientCohort)){
  if(any(grepl(pattern = PatientCohort$Sample_ID[i], all_files))){
    PatientCohort$path[i] = grep(pattern = PatientCohort$Sample_ID[i], x = all_files, value = T)
  } else {
    PatientCohort$path[i] = NA
  }
}


## manually fetch the files from xjuno:: tedious work; Subhi
#' fetch_files = PatientCohort$path[which(grepl('juno', x = PatientCohort$path))]
#' IMPACT_countfiles = list.files(path = 'CSF_IMPACT_COUNTMATRIX/', full.names = T, recursive = T)

#' create folders and fetch all the files within the folders
for(id in unique(PatientCohort$Sample_Tag)){
  ids = PatientCohort[which(PatientCohort$Sample_Tag == id), ]
  dir.create(path = paste0(getwd(), '/', id))
  
  for(file in 1:nrow(ids)){
    if(any(grepl(pattern = ids$Sample_ID[file], x = all_files))){
      filename = grep(pattern = ids$Sample_ID[file], x = all_files, value = T)
      system(command = paste0('cp ', filename, ' ', paste0(getwd(), '/', id, '/')))
    }
  }
}



#' load all count-matrices simultanously
#' working on one sample as an example!
IDs = list.files(path = getwd(), pattern = 'IM-GBM', full.names = T, all.files = T)

for(sampleid in 1:length(IDs)){
  input = lapply(list.files(path = IDs[sampleid], full.names = T), 
                 function(x) facets::readSnpMatrix(filename = x))
  input_names = list.files(path = IDs[sampleid])
  names(input) = input_names
  print(IDs[sampleid])
  
  comp = data.frame()
  for(countmatrix in 1:length(input)){
    facets_processed = facetsSuite::run_facets(read_counts = input[[countmatrix]],
                                               cval = 100,
                                               genome = 'hg19')
    
    QC = facets_fit_qc(facets_processed)$facets_qc
    QC_waterfall = facets_fit_qc(facets_processed)$waterfall_filter_pass
    dipLogR_flag = facets_fit_qc(facets_processed)$dipLogR_flag
    dipLogR = facets_processed$dipLogR
    purity = facets_processed$purity
    ploidy = facets_processed$ploidy
    gene_level = facetsSuite::gene_level_changes(facets_processed, genome = 'hg19')
    
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
      
      median_cnlr_observed = gene_level[which(gene_level$gene == gene), 'median_cnlr_seg']
      median_cnlr_observed = median_cnlr_observed - dipLogR
      cf.em = gene_level[which(gene_level$gene == gene), 'cf.em']
      #' assign clonality
      clonality = ifelse(cf.em >= (purity * 0.8), 'clonal', 'subclonal')
      gene_cn = gene_level[which(gene_level$gene == gene), 'cn_state']
      
      Alteration = ifelse(abs(median_cnlr_observed - dipLogR) >= 0.4 &
                            median_cnlr_observed >= 0.2, 'Amplification',
                          ifelse(abs(median_cnlr_observed - dipLogR) >= 0.4 &
                                   median_cnlr_observed <= -0.2, 'Deletion', 'Diploid'))
      
      out = data.frame(Sample = names(input)[countmatrix],
                       order = countmatrix,
                       gene = gene,
                       gene_cn = gene_cn,
                       Alteration = Alteration,
                       median_cnlr_observed = median_cnlr_observed,
                       cf.em,
                       clonality = clonality,
                       purity = purity,
                       QC = QC,
                       waterfall_filter_pass = QC_waterfall,
                       dipLogR_flag = dipLogR_flag,
                       dipLogR = dipLogR)
      comp = rbind(comp, out)
    }
  }
  write.table(comp, file = paste0(getwd(), '/', basename(IDs[sampleid]), '/', basename(IDs[sampleid]), '_summary.txt'), sep = '\t', row.names = F)
  rm(comp, out, gene_level, input)
}

a = read.csv('IM-GBM-126/IM-GBM-126_summary.txt', sep = '\t')

## Visualization:
IDs = list.files(path = getwd(), pattern = 'IM-GBM', full.names = T, all.files = T)

for(sampleid in 1:length(IDs)){
  filename = list.files(path = IDs[sampleid], pattern = 'summary.txt', full.names = T)
  file = read.csv(file = filename, sep = '\t')
  file$subclonal = ifelse(file$cf.em >= (file$purity * 0.8), 'clonal', 'subclonal')
  print(filename)
  
  QC = any(file$waterfall_filter_pass == F)
  QC = ifelse(QC == F, T, F)
  dip = any(file$dipLogR_flag == F)
  
  #' change the names:
  breaks = length(unique(file$Sample))
  labels = unique(file$Sample)
  for(i in 1:length(labels)){
    if(grepl(pattern = 'counts', x = labels[i])){
      labels[i] = substr(x = labels[i], start = 17, stop = 33)
    } else {
      labels[i] = substr(x = labels[i], start = 1, stop = 25)
    }
  }
  
  #' plot
  p = ggplot(file, aes(x = order, y = gene, fill = median_cnlr_observed, color = Alteration)) +
    geom_tile(aes(color = Alteration, width = 0.8, height = 0.8), size = 0.45) +
    geom_text(aes(label = ifelse(subclonal == 'clonal', '', '*')), col = 'black') +
    scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red', name = 'CnLR') +
    scale_color_manual(values = c('Amplification' = 'red',
                                  'Deletion' = 'blue',
                                  'Diploid' = 'grey95'), name = '') +
    scale_x_continuous(breaks = seq(1, breaks, 1),
                       labels = labels,
                       sec.axis = dup_axis()) +
    theme_bw() +
    theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0),
          axis.text.x = element_blank(),
          axis.text = element_text(color = 'black')) +
    labs(x = '', y = '', title = paste0(basename(IDs[sampleid]), '_', QC)) +
    
    coord_equal()
  
  ggsave(filename = paste0(IDs[sampleid], '/', basename(IDs[sampleid]), '_plot.pdf'),
         plot = p, device = 'pdf', width = 8, height = 8, units = 'in')
}



#### fetch plot and put them in one folder
IDs = list.files(path = getwd(), pattern = 'IM-GBM', full.names = T, all.files = T)

for(sampleid in 1:length(IDs)){
  filename = list.files(path = IDs[sampleid], pattern = '.pdf', full.names = T)
  system(command = paste0('cp ', filename, ' ', '/Users/chriskreitzer/Documents/MSKCC/Subhi/CSF/plots_06222022/'))
}


#' out