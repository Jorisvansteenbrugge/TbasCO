
#' Pre-processing of RNAseq data
#' @export
Pre_process_input <- function(filepath, normalize = TRUE){
  RNAseq.table <- read.csv2(filepath)

  # To be sure the corresponding collumns are converted to their correct type
  # TODO: also for sample columns
  RNAseq.table$Bin <- as.character(RNAseq.table$Bin)


  # Identify some basic features of the table
  RNAseq.features <- Get_matrix_features(RNAseq.table)

  # Combine the table and the features in one object
  RNAseq.data <- list("table"    = RNAseq.table,
                      "features" = RNAseq.features)


  return(RNAseq.data)


  #RNAseq.table <- Normalization(RNAseq.table)
}



Get_matrix_features <- function(RNAseq.table){
  .Get_sample_columns <- function(column.names){
    column_contains_Sample <- sapply(column.names,
                                     function(x) grepl("Sample", x, fixed = T))

    return( as.numeric(which(column_contains_Sample == TRUE)) )
  }

  sample_columns <- .Get_sample_columns(colnames(RNAseq.table))

  return(list(
    "bins"                = sort(unique(RNAseq.table$Bin)),
    "sample_columns"      = sample_columns,
    "rank_columns"        = ncol(RNAseq.table)+1:length(sample_columns),
    "KO_presense_absense" = NA
  ))
}



Normalization <- function(RNAseq.table, normalize_function){

  if(missing(normalize_function)){
    # Perform the default normalization
  }
}

#
# edgeRmethods <- function(SS, SE, method_name, RNAseq_Annotated_Matrix){
#   library(edgeR)
#   RNAseq_Annotated_Matrix1<-DGEList(counts=RNAseq_Annotated_Matrix[, SS:SE],lib.size=library_size)
#   norm_factors <- calcNormFactors(RNAseq_Annotated_Matrix1, method=method_name)
#   RNAseq_Annotated_Matrix[, SS:SE] <- t(t(as.matrix(norm_factors)) / norm_factors$samples[,3])
#   return(RNAseq_Annotated_Matrix)
# }



# setClass('General_features', representation(name                = 'character',
#                                             bins                = 'vector',
#                                             sample_columns      = 'vector',
#                                             rank_columns        = 'vector',
#                                             KO_presence_absence = 'matrix'),
#          prototype(name = 'matrix_features'))
