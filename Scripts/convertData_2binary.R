##-----------------------------------------------------------------------------
## Convert alteration data into binary matrix
## ----------------------------------------------------------------------------
## 
#' @name mutation_matrix
#' 
#' @description Take the alteration data (mutations) in the underlying maf format
#' (github) and create a binary matrix. Within the binary matrix mutations the
#' Variant Classifications (missense, nonsense, etc.) are added for respective genes; 
#' required for subsequent processing (eg Oncoprint - complexheatmap()). Note that this
#' function is a stand-alone function. Full integration can be found in function - onco_matrix
#' 
#' @export
#' @param path_data_repository Absolute path where data is stored. Mutations
#' CNA, gene panel
#' @param sampleIDs character vector indicating which samples to analyse. If NULL
#' all the samples in the input mutation/cna matrix will be converted. This may take long
#' 
#' @return list object containing two data frames
#' @return alteration_matrix: as described above
#' @return data.cna: input cna object subsetted for sampleIDs
#' 
#' @author chris-kreitzer
#' 
#' start: 03/17/2022


   
#' @import needed packages
library(vroom)
library(tidyr)


##-----------------------------------------------------------------------------
## Convert mutation data
## ----------------------------------------------------------------------------
.mutation_matrix = function(path_data_repository, sampleIDs = NULL){
  
  message('------------------------------\n')
  message('Checking if neccessary files are available\n')
  
  #' data objects
  data_objects = list.files(path = path_data_repository, full.names = T)
  
  stopifnot('Mutation data is not present in data repository' = 
              any(grepl(pattern = '*.mutations.*', x = data_objects)))
  stopifnot('CNA data is not present in data repository' = 
              any(grepl(pattern = '*CNA.*', x = data_objects)))
  stopifnot('Gene panel data not present in data repository' = 
              any(grepl(pattern = '*gene_panel.*', x = data_objects)))
  
  #' load data
  mutations = read.csv(file = data_objects[grepl(pattern = '*.mutations.*', data_objects)], sep = '\t')
  mutations = as.data.frame(mutations)
  cna = read.csv(file = data_objects[grepl(pattern = '*.CNA.*', data_objects)], sep = '\t')
  cna = as.data.frame(cna)
  gene_panel_505 = suppressMessages(vroom::vroom(file = data_objects[grepl(pattern = '*.gene_panel_impact505.*', data_objects)]))
  gene_panel_505 = colnames(gene_panel_505)
  if(length(gene_panel_505) != 505){
    message('Please install <vroom>[install.packages("vroom")] to load the gene panel data properly')
  }
  
  gene_panel_468 = suppressMessages(vroom::vroom(file = data_objects[grepl(pattern = '*.gene_panel_impact468.*', data_objects)]))
  gene_panel_468 = colnames(gene_panel_468)
  gene_panel = data.frame(genes = c('Sample.ID', unique(c(gene_panel_468, gene_panel_505))))
  
  #' subset mutations if SampleIDs are provided
  if(!is.null(sampleIDs)){
    ids = sampleIDs
    mutations = mutations[which(mutations$Tumor_Sample_Barcode %in% ids),, drop = F]
    cna = cna[which(cna$SAMPLE_ID %in% ids),, drop = F]
  }
  
  
  #' Alteration matrix
  alteration_matrix = setNames(data.frame(matrix(ncol = length(gene_panel$genes), 
                                                 nrow = 0)), gene_panel$genes)
  
  
  message('------------------------------------------------\n')
  message('Starting to create binary matrix...\n')
  
  counter = 1
  
  for(i in unique(mutations$Tumor_Sample_Barcode)){
    data.mut.sub = mutations[mutations$Tumor_Sample_Barcode == i, ]
    
    if(nrow(data.mut.sub) != 0){
      
      for(j in unique(data.mut.sub$Hugo_Symbol)){
        if(j %in% colnames(alteration_matrix)){
          alteration_matrix[counter, j] = data.mut.sub$Variant_Classification[which(data.mut.sub$Hugo_Symbol == j)][1]
          
        }
      }
    }
    
    alteration_matrix[counter, 'Sample.ID'] = i
    counter = counter + 1
  }
  
  return(list(alteration_matrix = alteration_matrix,
              data_cna = cna))
}


##
#' test = mutation_matrix(path_data_repository = 'Data/', sampleIDs = unique(sam$SAMPLE_ID))



##-----------------------------------------------------------------------------
## onco_matrix
## ----------------------------------------------------------------------------

#' @name onco_matrix
#' @description Take the binary mutation matrix from @mutation_matrix and further
#' add information fetched from cna and fusion data. 
#'
#' @export
#' @param path_data_repository Absolute path where data is stored. Mutations
#' CNA, gene panel
#' @param sampleIDs character vector indicating which samples to analyse. If NULL
#' all the samples in the input mutation/cna matrix will be converted. This may take long
#'
#' @return Dataframe with alterations
#' 
#' @author chris-kreitzer


onco_matrix = function(path_data_repository, sampleIDs = NULL){
  
  message('---------------------------------\n')
  message('Currently all mutations (either oncogenic or VUS) are included.\n
          If you want to exclusively work with oncogenic ones, please subset the input data\n
          OR: contact me to update the code chunk!\n')
  
  #' convert mutations; see @aliases .mutation_matrix()
  data_out = .mutation_matrix(path_data_repository = path_data_repository, 
                                   sampleIDs = sampleIDs)
  
  #' add cna data
  message('----------------------------\n')
  message('Adding CNA data...\n')
  
  cna = as.data.frame(data_out$data_cna)
  mut = as.data.frame(data_out$alteration_matrix)
  
  for(i in unique(cna$SAMPLE_ID)){
    if(i %in% mut$Sample.ID){
      sub = cna[which(cna$SAMPLE_ID == i), ]
      for(j in unique(sub$HUGO_SYMBOL)){
        if(j %in% colnames(mut)){
          mut[which(mut$Sample.ID == i), which(colnames(mut) == j)] = paste(mut[which(mut$Sample.ID == i), which(colnames(mut) == j)],
                                                                  sub$ALTERATION[which(sub$HUGO_SYMBOL == j)], sep = ';')
        } else next
      }
      rm(sub)
    } else next
  }
  
  
  # rm NA values
  # remove NA values
  k = mut
  for(i in 1:ncol(k)){
    k[ , i] = ifelse(grepl('NA', k[, i]), 
                     sub(pattern = 'NA;', replacement = '', k[,i]), k[,i])
  }
  
  
  # change mutation classes
  truncating.classes = c('Nonsense_Mutation','Splice_Site','Frame_Shift_Del','Frame_Shift_Ins')
  message('Currently: ', paste(truncating.classes, collapse = ', '), 'are considered as truncating mutations\n')
  
  for(i in 1:ncol(k)){
    k[ , i] = ifelse(grepl('Nonsense_Mutation', k[, i]), sub(pattern = 'Nonsense_Mutation', replacement = 'Truncating_Mutation', k[,i]),
                     ifelse(grepl('Splice_Site', k[, i]), sub(pattern = 'Splice_Site', replacement = 'Truncating_Mutation', k[,i]),
                            ifelse(grepl('Frame_Shift_Del', k[, i]), sub(pattern = 'Frame_Shift_Del', replacement = 'Truncating_Mutation', k[,i]),
                                   ifelse(grepl('Frame_Shift_Ins', k[, i]), sub(pattern = 'Frame_Shift_Ins', replacement = 'Truncating_Mutation', k[,i]),
                                          ifelse(grepl('In_Frame_Ins', k[, i]), sub(pattern = 'In_Frame_Ins', replacement = 'Inframe_Mutation', k[,i]),
                                                 ifelse(grepl('In_Frame_Del', k[, i]), sub(pattern = 'In_Frame_Del', replacement = 'Inframe_Mutation', k[,i]),
                                                        ifelse(grepl('3\'Flank', k[, i]), sub(pattern = '3\'Flank', replacement = 'VUS', k[,i]),
                                                               ifelse(grepl('5\'UTR', k[, i]), sub(pattern = '5\'UTR', replacement = 'VUS', k[,i]),
                                                                      ifelse(grepl('Intron', k[, i]), sub(pattern = 'Intron', replacement = 'VUS', k[,i]),
                                                                             ifelse(grepl('Splice_Region', k[, i]), sub(pattern = 'Splice_Region', replacement = 'VUS', k[,i]),
                                                                                    ifelse(grepl('Nonstop_Mutation', k[, i]), sub(pattern = 'Nonstop_Mutation', replacement = 'VUS', k[,i]),
                                                                                           ifelse(grepl('Translation_Start_Site', k[, i]), sub(pattern = 'Translation_Start_Site', replacement = 'VUS', k[,i]),
                                                                                                  ifelse(grepl('5\'Flank', k[, i]), sub(pattern = '5\'Flank', replacement = 'VUS', k[,i]), k[,i])))))))))))))
  }
  
  Somatic_alteration_matrix = k
  
  # look 
  for(i in 2:ncol(Somatic_alteration_matrix)){
    for(j in 1:nrow(Somatic_alteration_matrix)){
      Somatic_alteration_matrix[j, i] = ifelse(is.na(Somatic_alteration_matrix[j, i]), "", Somatic_alteration_matrix[j, i])
    }
  }
  
  
  #' add fusion data
  message('----------------------------\n')
  message('Adding Fusion data...\n')
  
  data_objects = list.files(path = path_data_repository, full.names = T)
  
  stopifnot('Fusion data is not present in data repository' = 
              any(grepl(pattern = '*.fusions.*', x = data_objects)))
  
  #' load data
  Fusions = suppressWarnings(read.csv(file = data_objects[grepl(pattern = '*.fusions.*', data_objects)], sep = '\t'))
  data.fusion = as.data.frame(Fusions)
  
  for(i in unique(data.fusion$Tumor_Sample_Barcode)){
    fusion_sub = data.fusion[data.fusion$Tumor_Sample_Barcode == i, ]
    
    if(i %in% Somatic_alteration_matrix$Sample.ID){
      for(j in unique(fusion_sub$Hugo_Symbol)){
        if(j %in% colnames(Somatic_alteration_matrix)[-1]){
          Somatic_alteration_matrix[Somatic_alteration_matrix$Sample.ID == i, j] = paste(Somatic_alteration_matrix[Somatic_alteration_matrix$Sample.ID == i, j], 'Fusion', sep = ';')
        }
      }
    }
  }
  
  return(Somatic_alteration_matrix)
}


#' test
#' x = onco_matrix(path_data_repository = 'Data/', sampleIDs = unique(samp$SAMPLE_ID))

#' out



