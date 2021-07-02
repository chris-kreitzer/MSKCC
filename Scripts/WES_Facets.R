library(facetsSuite)
library(dplyr)

all.files = list.files(path = '~/Subhi_WES', pattern = '*.gz', full.names = T)

for(i in unique(all.files)){
data.in = facetsSuite::read_snp_matrix(i)
data.out = facetsSuite::run_facets(data.in, cval = 1000, ndepth = 35, genome = 'hg19', seed = 1000)

create_legacy_output =  function(facets_output, directory, sample_id) {
    
    #sample_id = ifelse(sel_run_type == '', sample_id, paste0(sample_id, '_', sel_run_type))
    output_prefix = paste0(directory, '/', sample_id)
    
    ### create cncf.txt
    facets_output$segs %>%
        mutate(ID=sample_id,
               loc.start = start,
               loc.end = end) %>%
        select(ID, chrom, loc.start, loc.end, seg, num.mark, nhet, cnlr.median, mafR, 
               segclust, cnlr.median.clust, mafR.clust, cf, tcn, lcn, 
               cf.em, tcn.em, lcn.em) %>%
        write.table(file=paste0(output_prefix, '.cncf.txt'), quote=F, row.names=F, sep='\t')
    
    ### create .Rdata
    out =
        list(
            jointseg = facets_output$snps,
            out = facets_output$segs %>% 
                select(chrom, seg, num.mark, nhet, cnlr.median, mafR, 
                       segclust, cnlr.median.clust, mafR.clust, cf, tcn, lcn),
            nX=23,
            chromlevels = c(1:22, "X"),
            dipLogR = facets_output$dipLogR,
            alBalLogR = facets_output$alBalLogR,
            IGV = NULL
        )
    fit = 
        list(
            loglik = facets_output$loglik,
            purity = facets_output$purity,
            ploidy = facets_output$ploidy,
            dipLogR = facets_output$dipLogR,
            seglen = -1,
            cncf = facets_output$segs,
            emflags = facets_output$em_flags
        )
    save(out, fit, file = paste0(output_prefix, ".Rdata"), compress=T)
}

create_legacy_output(facets_output = data.out, directory = '~/Subhi_WES', sample_id = basename(i))
}

