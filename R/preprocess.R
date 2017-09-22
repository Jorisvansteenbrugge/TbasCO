
#' Pre-processing of RNAseq data
#' @param filepath RNAseq data file
#' @param normalize.method Variable
#' @export
#' @author JJM van Steenbrugge
Pre_process_input <- function(filepath, normalize.method = FALSE,
                              filter.method = "stdev"){
  RNAseq.table <- read.csv2(filepath)
  # Identify some basic features of the table
  RNAseq.features <- .Get_matrix_features(RNAseq.table)

  # To be sure the corresponding collumns are converted to their correct type
  RNAseq.table$Bin <- as.character(RNAseq.table$Bin)

  for(sample_col in RNAseq.features$sample_columns){
    RNAseq.table[, sample_col] <- as.numeric(RNAseq.table[, sample_col])
  }

  if(normalize.method != FALSE){
    RNAseq.table <- .Normalization(RNAseq.table,
                                   RNAseq.features,
                                   normalize.method)
  }
  RNAseq.table   <- .filter(RNAseq.table,
                            RNAseq.features$sample.columns,
                            filter.method)

  # Combine the table and the features in one object
  return(list("table"    = RNAseq.table,
              "features" = RNAseq.features))
}

#' Aumatically identify features based on column names
#' @param RNAseq.table
#' @return A list containing a vector of unique bins, a vector with the sample column
#' numbers, a vector containing the rank column numbers and a annotation_presence_absence matrix
#' @author JJM van Steenbrugge
.Get_matrix_features <- function(RNAseq.table){

  .Get_sample_columns <- function(column.names){

    column_contains_Sample <- sapply(column.names,
                                     function(x) grepl("Sample", x, fixed = T) )

    return( as.numeric(which(column_contains_Sample == TRUE)) )
  }

  sample.columns <- .Get_sample_columns(colnames(RNAseq.table))
  bins <- sort(unique(RNAseq.table$Bin))
  return(list(
    "bins"                        = bins,
    "sample.columns"              = sample.columns,
    "rank.columns"                = ncol(RNAseq.table)+1:length(sample.columns),
    "annotation_presence_absence" = .Get_annotation_presence_absence(RNAseq.table, bins)
  ))
}

#' Calculate which KO is present for each Bin
#' @param RNAseq.table
#' @param bins a vector with all unique bins
#' @author JJM van Steenbrugge
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

#' Normalize RNAseq data
#' @param RNAseq.table
#' @param normalize.method variable that is either a string (default method) or a function to be applied
#' for normalization of RNAseq data.
#' @author JJM van Steenbrugge
#' @return
.Normalization <- function(RNAseq.table, RNAseq.features, normalize.method){
  if(typeof(normalize.method) == "closure"){
    normalize.method(RNAseq.table)
  } else{
    require(edgeR)

    sample.columns <- RNAseq.features$sample.columns

    DGEobj <- DGEList(counts=RNAseq.table[, sample.columns])

    norm_factors      <- calcNormFactors(DGEobj, method='TMM')

    RNAseq.table[, sample.columns] <- t(t(as.matrix(norm_factors)) /
                                            norm_factors$samples[,3])
    return(RNAseq.table)
  }
}

#' Remove data rows that have a stdev of zero.
#' @param RNAseq.table
#' @param sample.columns a vector with the collumns that contain RNAseq data.
#' @author JJM van Steenbrugge
.filter <- function(RNAseq.table, sample.columns,
                    filter.method){
  if(filter.method == "stdev"){

    stdevs <- apply(RNAseq.table[, sample.columns], 1, sd)

    return(RNAseq.table[which(stdevs > 0),])

  } else if(filter.method == 'MAD'){

    .Get.MAD <- function(row){

      avg <- mean(row)
      mean.abs.deviations <- sapply(row, function(x) abs(x - avg))

      return( sum(mean.abs.deviations) / length(row) )

    }
    MADs <- apply(RNAseq.table[, sample.columns], 1, .Get.MAD)
    return( RNAseq.table[which(MADs > 1 ), ] ) #maybe a bigger number than 1
  } else if(typeof(filter.method) == 'closure') {
    filter.method(RNAseq.table)
  }
}

.Create.Rank.Columns <- function(RNAseq.table, RNAseq.features){
  header <- sapply( 1: length(RNAseq.features$rank.columns),
                    function(x) paste("rank", x, sep = "") )


  for ( sample_col in RNAseq.features$sample.columns) {
    for( bin in RNAseq.features$bins ) {
      current_bin <-  RNAseq.table[ which( RNAseq.table$Bin == bin ), ]
    }
  }

  for (sample_col in RNAseq.features$sample.columns) {
    for (s in 1: length(all_bins)) {
      matrix_bins <- RNAseq.table[which(RNAseq.table$Bin == all_bins[s]), ]
      matrix_bins[, matrix_features@bin_column_index + i]<- rank(-matrix_bins[, sample_col + 1],
                                                                   na.last     = TRUE,
                                                                   ties.method = "random") /
                                                            max(rank(-matrix_bins[, sample_col + 1],
                                                                   na.last     = TRUE,
                                                                   ties.method = "random"))

      RNAseq.table[which(RNAseq.table$Bin == all_bins[s]), ] <- matrix_bins
    }
  }
}
