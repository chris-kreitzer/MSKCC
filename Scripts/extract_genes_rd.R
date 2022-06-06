extract_genes_rd = function(gene, 
                            countmatrix,
                            padding = 50){
  
  message(paste('Left/Right padding of', padding, 'bp considered'))
  
  selected_gene = genes_hg19[which(genes_hg19$gene == gene), ]
  countmatrix = countmatrix
  countmatrix_chromosome = countmatrix[which(countmatrix$Chromosome == selected_gene$chrom), ]
  gene_rd = countmatrix_chromosome[which(countmatrix_chromosome$Position >= (as.numeric(selected_gene$start) - as.numeric(padding)) &
                                           countmatrix_chromosome$Position <= (as.numeric(selected_gene$end) + as.numeric(padding))), ]
  
  return(gene_rd)
}