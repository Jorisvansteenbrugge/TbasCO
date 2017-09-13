
#' Pre-processing of RNAseq data
#' @export
Pre_process_input <- function(filepath, normalize.var = FALSE){
  RNAseq.table <- read.csv2(filepath)

  # To be sure the corresponding collumns are converted to their correct type
  # TODO: also for sample columns
  RNAseq.table$Bin <- as.character(RNAseq.table$Bin)


  # Identify some basic features of the table
  RNAseq.features <- .Get_matrix_features(RNAseq.table)

  # Combine the table and the features in one object
  RNAseq.data <- list("table"    = RNAseq.table,
                      "features" = RNAseq.features)


  return(RNAseq.data)


  if(! normalize.var){
    # Skip
  } else{
    RNAseq.data$table <- .Normalization(RNAseq.data$table, normalize.var)
  }

  print("filter")
  RNAseq.data$table <- .Remove_sd_zero(RNAseq.data)
}



.Get_matrix_features <- function(RNAseq.table){

  .Get_sample_columns <- function(column.names){

    column_contains_Sample <- sapply(column.names,
                                     function(x) grepl("Sample", x, fixed = T) )

    return( as.numeric(which(column_contains_Sample == TRUE)) )
  }

  sample_columns <- .Get_sample_columns(colnames(RNAseq.table))
  bins <- sort(unique(RNAseq.table$Bin))
  return(list(
    "bins"                        = bins,
    "sample_columns"              = sample_columns,
    "rank_columns"                = ncol(RNAseq.table)+1:length(sample_columns),
    "annotation_presence_absence" = .Get_annotation_presence_absence(RNAseq.table, bins)
  ))
}

# Calculate which KO is present for each Bin
.Get_annotation_presence_absence <- function(RNAseq.table, bins){
  annotation_terms <- unique(RNAseq.table$Annotation)

  .Lookup_bin <- function(bin){

    present_bins <- RNAseq.table[which(RNAseq.table$Bin == bin), "Annotation"]

    return(annotation_terms %in% present_bins)
  }

  annotation_presence_absence <- sapply(bins, .Lookup_bin)


  colnames(annotation_presence_absence) <- bins
  rownames(annotation_presence_absence) <- annotation_terms

  return(annotation_presence_absence)
}

.Normalization <- function(RNAseq.table, normalize.var){
  if(typeof(normalize.var) == "closure"){
    normalize.var(RNAseq.table)
  } else{
    #default normalization
  }
}


.Remove_sd_zero <- function(RNAseq.data){
  sample_cols <- RNAseq.data$features$sample_columns

  data_matrix <- matrix(apply(RNAseq.data$table[, sample_cols], 2, as.numeric),
                        ncol = length(sample_cols) )

  stdev_vector <- Calc_stdev_rows(data_matrix)

  RNAseq.data$table <- cbind(RNAseq.data$table, stdev_vector)
  print("filter baby")
  return( RNAseq.data$table[which(RNAseq.data$table[, ncol(RNAseq.data$table)] != 0), ] )


}
