##-----------------------------------------------------------------------------
## Facets Re-Analysis of solid and CSF samples
##-----------------------------------------------------------------------------
##
## start: 05/04/2022
## chris-kreitzer
##

clean()
gc()
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')

library(facets)
library(patchwork)


IMPACT505 = vroom::vroom('data_gene_panel_impact505.txt', delim = '\t')
IMPACT505 = as.character(colnames(IMPACT505))
IMPACT468 = vroom::vroom('~/Documents/GitHub/PARP-Prostate-Cancer/Data/data_gene_panel_impact468.txt', delim = '\t')
IMPACT468 = as.character(colnames(IMPACT468))
genes_IMPACT = unique(union(IMPACT468, IMPACT505))



## Folder-wise:
folders = list.files(path = getwd(), pattern = '^0.*', full.names = T, recursive = F)

all_out = data.frame()
for(folder in unique(folders)){
  try({
    dir.create(path = paste0(folder, '/FacetsFits'), mode = 777)
    countFiles = list.files(path = folder, full.names = T)
    countFiles = countFiles[-grep('FacetsFits', countFiles)]
    
    sample_summary = data.frame()
    for(files in unique(countFiles)){
      print(basename(files))
      countmatrix = facetsSuite::read_snp_matrix(input_file = files)
      countmatrix = countmatrix[, c(1,2,3,5,4,6)]
      
      #' normal Facets
      fit = facetsSuite::run_facets(read_counts = countmatrix,
                                    cval = 100, 
                                    genome = 'hg19', 
                                    seed = 100, 
                                    snp_nbhd = 250)
      
      #' QC measures
      qc = facets_fit_qc(fit)
      qc_fail = names(grep(pattern = 'FALSE', qc, value = T))
      qc_fail = qc_fail[!qc_fail %in% c('wgd', 'facets_qc')]
      
      #' gene-level calls
      gene_call = facetsSuite::gene_level_changes(facets_output = fit, 
                                                  genome = 'hg19')
      
      gene_call = gene_call[which(gene_call$gene %in% genes_IMPACT), ]
      gene_call = gene_call[!gene_call$filter %in% c('suppress_segment_too_large'), ]
      
      sub_summary = data.frame(sample = basename(files),
                               qc = qc$facets_qc,
                               notes = paste(qc_fail, collapse = ';'))
      
      sample_summary = rbind(sample_summary, sub_summary)
      
      #' merge all objects
      out = list(facets_fit = fit, 
                 gene_level = gene_call,
                 qc = qc)
      saveRDS(object = out, file = paste0(folder, '/FacetsFits/', basename(files), '_fitted.rds'))
      rm(countmatrix, fit, qc, qc_fail, gene_call, out)
    }
    all_out = rbind(all_out, sample_summary)
  })
} 

write.table(all_out, file = 'SampleSummary_QC.txt', sep = '\t', quote = F, row.names = F)

