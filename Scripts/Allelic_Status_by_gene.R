##----------------+
## Allelic-status for
## arbitrary genes
##----------------+

## start: 11/08/2022
## revision: 11/10/2022
## chris-kreitzer

allelic_status = function(samples = NULL,
                          gene = NULL,
                          copy_number_data = NULL,
                          mutation_data = NULL,
                          cna_filter = c('PASS', 'RESCUE')){
  
  message('\nPlease provide samples as characters (e.g. c(, , ))\n
          Please provide gene as single string (e.g. TP53)\n
          Please provide mutation in OncoKB-annotated format\n
          Please provide CNAs from Facets and at the GENE-LEVEL')
  
  gene_summary_out = data.frame()
  samples = as.character(samples)
  mutations = mutation_data[which(mutation_data$Tumor_Sample_Barcode %in% samples), ]
  cna_sample_match = sapply(unique(copy_number_data$sample),
                            function(x) if (grepl(paste(samples, collapse = '|'),  x, ignore.case = TRUE)) {"Present"} else {"Absent"})
  cna_sample_names = names(cna_sample_match)[which(cna_sample_match == 'Present')]
  cna_s = copy_number_data[which(copy_number_data$sample %in% cna_sample_names), ]
  
  #---------------+
  # concentrate GOI
  #---------------+
  if(!is.null(gene)){
    gene = gene
    mutation_goi = mutations[which(mutations$Hugo_Symbol == gene), ]
    cna_goi = cna_s[which(cna_s$gene == gene), ]
    cna_goi$sample = substr(cna_goi$sample, start = 1, stop = 17)
  } else {
    return(list(mutations = mutations,
                cna = cna_s))
  }
  
  #---------------+
  # loop through samples
  #---------------+
  for(i in unique(samples)){
    try({
      if(i %in% mutation_goi$Tumor_Sample_Barcode){
        muts = mutation_goi[which(mutation_goi$Tumor_Sample_Barcode == i), ]
        
        #' single; or multiple mutations?
        mutation_specific = paste0(gene, ';', muts$HGVSp)
        mutation_specific_p = ifelse(length(mutation_specific) > 1, 
                                   paste(mutation_specific, collapse = '|'), mutation_specific)
        
        #' Germline or Somatic
        mutation_status = muts$Mutation_Status
        mutation_status = ifelse(length(mutation_status) > 1, 
                                 paste(mutation_status, collapse = ';'), mutation_status)
        
        #' oncogenicity?
        mutation_effect = muts$ONCOGENIC
        mutation_effect = ifelse(length(mutation_effect) > 1, 
                                 paste(mutation_effect, collapse = ';'), mutation_effect)
        
        #' cis-trans mutations?
        composite_mutation = ifelse(length(mutation_specific) > 1, 'check: cis/trans mutations?', 'no')
        
      } else {
        mutation_specific_p = 'none'
        mutation_status = NA
        mutation_effect = NA
        composite_mutation = NA
      }
      
      #-------------+
      # CNA data
      #-------------+
      if(i %in% cna_goi$sample){
        cna = cna_goi[which(cna_goi$sample == i), ]
        tcn = cna$tcn
        lcn = cna$lcn
        cn_state = cna$cn_state
        filter = cna$filter
      } else {
        tcn = NA
        lcn = NA
        cn_state = NA
        filter = NA
      }
      
      #-------------+
      # gather info and return
      #-------------+
      out = data.frame(id = i,
                       gene = gene,
                       mutation = mutation_specific_p,
                       n_mutations = ifelse(mutation_specific_p == 'none', 0, nrow(muts)),
                       status = mutation_status,
                       effect = mutation_effect,
                       composite = composite_mutation,
                       tcn = tcn,
                       lcn = lcn,
                       cn_state = cn_state,
                       filter = filter)
      gene_summary_out = rbind(gene_summary_out, out)
    })
  }

  #---------------+
  # Allelic-status
  #---------------+
  gene_summary_out$cna_AI = NA
  gene_summary_out$cna_AI_n = NA
  gene_summary_out = gene_summary_out[order(gene_summary_out$tcn, decreasing = T), ]
  
  #---------------+
  # Criteria(s)
  #---------------+
  
  # CNA balanced
  cna_balanced = function(data){
    tcn = data[1]
    lcn = data[2]
    # arguments
    arg1 = all(!is.na(tcn) & !is.na(lcn))
    arg2 = all(tcn > 0 & !is.na(lcn) & lcn != 0 & abs(diff(c(tcn, lcn))) %in% c(1,2,3))
    return(c(arg1, arg2))
  }
  
  args_cna = apply(gene_summary_out[, c('tcn', 'lcn')], 1, cna_balanced)
  args_short = apply(args_cna, 2, all)
  
  for(i in seq_along(args_short)){
    if(args_short[i]){
      gene_summary_out$cna_AI[i] = 'balanced'
      gene_summary_out$cna_AI_n[i] = 0
    } else {
      tcn = gene_summary_out$tcn[i]
      tcn = ifelse(is.na(tcn) | tcn == 0, 0, tcn)
      lcn = gene_summary_out$lcn[i]
      lcn = ifelse(is.na(lcn) | lcn == 0, 0, lcn)
      n_alleles = tcn - lcn
      gene_summary_out$cna_AI[i] = paste0('imbalanced')
      gene_summary_out$cna_AI_n[i] = n_alleles
    }
  }
  
  gene_summary_out$allelic_call = NA
  gene_summary_out$allelic_call = ifelse(is.na(gene_summary_out$tcn), 'check Facets fit', NA)

  #---------------+
  # subset clear cases
  #---------------+
  failed_samples = gene_summary_out[which(gene_summary_out$allelic_call == 'check Facets fit'), ]
  
  if(nrow(failed_samples) != 0){
    gene_summary_out = gene_summary_out[!gene_summary_out$id %in% unique(failed_samples$id), ]
    for(id in 1:nrow(gene_summary_out)){
      if(gene_summary_out$filter[id] %in% cna_filter){
        a_call = ifelse(gene_summary_out$cna_AI[id] == 'balanced' &
                          gene_summary_out$cna_AI_n[id] %in% c(0,2,3) &
                          gene_summary_out$mutation[id] == 'none', 'wildtype',
                        ifelse(gene_summary_out$cna_AI[id] == 'balanced' &
                                 gene_summary_out$cna_AI_n[id] %in% c(0,2,3) &
                                 gene_summary_out$mutation[id] != 'none', 'monoallelic',
                               ifelse(gene_summary_out$cna_AI[id] == 'imbalanced' &
                                        gene_summary_out$cna_AI_n[id] == 0, 'biallelic',
                                      ifelse(gene_summary_out$cna_AI[id] == 'imbalanced' &
                                               gene_summary_out$cna_AI_n[id] == 1 &
                                               gene_summary_out$mutation[id] == 'none', 'monoallelic',
                                             ifelse(gene_summary_out$cna_AI[id] == 'imbalanced' &
                                                      gene_summary_out$cna_AI_n[id] == 1 &
                                                      gene_summary_out$mutation[id] != 'none', 'biallelic',
                                                    ifelse(gene_summary_out$cna_AI[id] == 'imbalanced' &
                                                             gene_summary_out$cna_AI_n[id] > 1 &
                                                             gene_summary_out$mutation[id] == 'none', 'monoallelic',
                                                           ifelse(gene_summary_out$cna_AI[id] == 'imbalanced' &
                                                                    gene_summary_out$cna_AI_n[id] > 1 &
                                                                    gene_summary_out$mutation[id] != 'none', 'check: expected_alt_copies', NA)))))))
        
        gene_summary_out$allelic_call[id] = a_call
      } else {
        gene_summary_out$allelic_call[id] = 'ambiguous:FacetsFilter'
      }
    }
    gene_summary_out = rbind(gene_summary_out, failed_samples)
    return(gene_summary_out)
  } else {
    for(id in 1:nrow(gene_summary_out)){
      if(gene_summary_out$filter[id] %in% cna_filter){
        a_call = ifelse(gene_summary_out$cna_AI[id] == 'balanced' &
                          gene_summary_out$cna_AI_n[id] %in% c(0,2,3) &
                          gene_summary_out$mutation[id] == 'none', 'wildtype',
                        ifelse(gene_summary_out$cna_AI[id] == 'balanced' &
                                 gene_summary_out$cna_AI_n[id] %in% c(0,2,3) &
                                 gene_summary_out$mutation[id] != 'none', 'monoallelic',
                               ifelse(gene_summary_out$cna_AI[id] == 'imbalanced' &
                                        gene_summary_out$cna_AI_n[id] == 0, 'biallelic',
                                      ifelse(gene_summary_out$cna_AI[id] == 'imbalanced' &
                                               gene_summary_out$cna_AI_n[id] == 1 &
                                               gene_summary_out$mutation[id] == 'none', 'monoallelic',
                                             ifelse(gene_summary_out$cna_AI[id] == 'imbalanced' &
                                                      gene_summary_out$cna_AI_n[id] == 1 &
                                                      gene_summary_out$mutation[id] != 'none', 'biallelic',
                                                    ifelse(gene_summary_out$cna_AI[id] == 'imbalanced' &
                                                             gene_summary_out$cna_AI_n[id] > 1 &
                                                             gene_summary_out$mutation[id] == 'none', 'monoallelic',
                                                           ifelse(gene_summary_out$cna_AI[id] == 'imbalanced' &
                                                                    gene_summary_out$cna_AI_n[id] > 1 &
                                                                    gene_summary_out$mutation[id] != 'none', 'check: expected_alt_copies', NA)))))))
        
        gene_summary_out$allelic_call[id] = a_call
        
      } else {
        gene_summary_out$allelic_call[id] = 'ambiguous:FacetsFilter'
      }
    }
    return(gene_summary_out)
  }
}


#' out