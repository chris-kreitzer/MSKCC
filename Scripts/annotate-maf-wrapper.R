#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(facetsSuite)
})

args = commandArgs(TRUE)
if (length(args) == 0) {
    message(paste(c('Run annotate-maf-wrapper.R --help for list of input arguments.',
                    'Note: For multiple samples, use --sample-mapping, for a single sample --facets-output argument can be used'), collapse = '\n'))
    quit()
}

parser = ArgumentParser(description = 'Annotate MAF with local copy number and CCF estimates.')
parser$add_argument('-m', '--maf-file', required = TRUE,
                    help = 'Mutations in MAF format')
parser$add_argument('-s', '--sample-mapping', required = FALSE,
                    help = 'Tab-separated file with header where first column contains sample ID and second column path to: .rds or .Rdata or .cncf.txt (the latter two from legacy) output files from facetsSuite')
parser$add_argument('-f', '--facets-output', required = FALSE,
                    help = 'A single facetsSuite .rds or .Rdata output file')
parser$add_argument('-a', '--facets-algorithm', required = FALSE,
                    default = 'em', help = 'Which FACETS algorithm to use [default %(default)s]')
parser$add_argument('-o', '--output', required = FALSE,
                    help = 'Output file [default input.ccf.maf]')
parser$add_argument('-p', '--parallel', required = FALSE,
                    default = FALSE,  help = 'Parallelize [default %(default)s]')
args = parser$parse_args()

output = ifelse(is.null(args$output),
                paste0(gsub('\\.[a-z]+$', '', args$maf_file), '.ccf.maf'),
                args$output)

if (is.null(args$sample_mapping) & is.null(args$facets_output)) {
    stop('Provide either a sample mapping or single-sample Facets output file.', call. = F)
}
if (!is.null(args$sample_mapping) & !is.null(args$facets_output)) {
    stop('Provide only one of --sample-mapping or --facets-output.', call. = F)
}

# Match sample IDs in input ---------------------------------------------------------------------------------------

maf = fread(args$maf_file, header = T)

if (!is.null(args$sample_mapping)) {
    
    # Read and check sample map overlap with input MAF
    sample_map = fread(args$sample_mapping, header = TRUE, col.names = c('sample', 'file'))
    if (any(duplicated(sample_map$sample))) {
        stop('Sample map contains duplicate sample names', call. = FALSE)
    }
    maf_unique_samples = length(setdiff(unique(maf$Tumor_Sample_Barcode), sample_map$sample))
    sample_map_unique_samples = length(setdiff(sample_map$sample, unique(maf$Tumor_Sample_Barcode)))
    common_samples = intersect(unique(maf$Tumor_Sample_Barcode), sample_map$sample)
    
    # Report non-overlapping samples
    if (maf_unique_samples > 0) {
        message(paste(maf_unique_samples, 'samples in input MAF not in sample map.'))
    }
    if (sample_map_unique_samples > 0) {
        message(paste(sample_map_unique_samples, 'samples in sample map not in input MAF.'))
    }
    
    message(paste('Annotating', length(common_samples), 'sample(s).'))
    
} else if (!is.null(args$facets_output)) {
    if (length(unique(maf$Tumor_Sample_Barcode)) != 1) {
        stop("Single-sample Facets input provided, but MAF contains multiple samples. Can't match on sample name.", call. = F)
    }
    common_samples = unique(maf$Tumor_Sample_Barcode)
    sample_map = list(file = args$facets_output, sample = common_samples)
    
    message(paste0('Annotating sample ', common_samples, '.'))
}

# Annotate --------------------------------------------------------------------------------------------------------
annotate_sample = function(sample_id) {
    sample_maf = maf[maf$Tumor_Sample_Barcode == sample_id,]
    facets_seg_file = sample_map$file[which(sample_map$sample == sample_id)]
    if ( grepl(".cncf(.edited)?.txt", facets_seg_file)) {
        rdata_pfx = basename(gsub(".cncf(.edited)?.txt","",facets_seg_file))
        rdata_files = list.files(dirname(facets_seg_file), pattern = paste0(rdata_pfx,'.(Rdata$|rds)$'), ignore.case=T, full.names = T)
        if (length(rdata_files) == 0) {
            rdata_files = list.files(dirname(facets_seg_file), pattern = 'rdata$|rds$', ignore.case=T, full.names = T)
        } 
        sample_purity = load_facets_output(rdata_files[1])$purity
        ccf_annotate_maf_legacy(sample_maf, facets_seg_file, sample_purity, args$facets_algorithm)
    } else {
        sample_facets = load_facets_output(facets_seg_file)
        ccf_annotate_maf(sample_maf, sample_facets$segs, sample_facets$purity, args$facets_algorithm)
    }
}

if (args$parallel == TRUE) {
    library(parallel)
    output_maf = mclapply(common_samples, annotate_sample, mc.cores = detectCores())
} else {
    output_maf = lapply(common_samples, annotate_sample)
}
output_maf = rbindlist(output_maf)

# Add back samples that were missing in sample map
if (any(!maf$Tumor_Sample_Barcode %in% sample_map$sample)) {
    output_maf = rbindlist(output_maf,
                           maf[!which(maf$Tumor_Sample_Barcode %in% sample_map$sample), ],
                           use.names = TRUE, fill = TRUE)
}

write.table(output_maf, output, quote = F, sep = '\t', col.names = T, row.names = F)