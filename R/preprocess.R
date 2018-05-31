
#' Pre-processing of RNAseq data
#' @description Wrapper of preprocess functions, including: Identifying Matrix
#' features, Normalization, Filtering for non informative lines
#' @param filepath RNAseq data file
#' @param normalize.method User defined method to normalize data
#' @param filter.low.coverage boolean expression to wether we should filter out
#' genomes with a low overall expression level.
#' @export
#' @examples
#' Pre_process_input(file.path, normalize.method = FALSE, filter.method = FALSE) #No normalization or filtering
#' Pre_process_input(file.path, normalize.method = TRUE, filter.method = 'MAD') # Default Normalization, and MAD filtering
#' Pre_process_input(file.path, normalize.method = custom_function_name, filter.method = FALSE) # Custom Normalization, no filtering
#' @return RNAseq.data A list object containin the RNAseq table with ranked expression collumns possible filtered and normalized,
#' together with some features (position of expression collumns, rank collumns, annotation database, a list of the genomes, and
#' a presence absence table of all annotations)
#' @author JJM van Steenbrugge
Pre_process_input <- function(file.path, annotation.db.path, normalize.method = FALSE,
                              filter.method = "stdev",
                              filter.low.coverage = T){

  RNAseq.table     <- read.csv2(file.path)

  # Identify some basic features of the table
  annotation.db    <- read.table(annotation.db.path, sep = "\t",
                                 quote = "", stringsAsFactors = F,
                                 header = T)
  annotation.db    <- Create.Module.groups(annotation.db)
  RNAseq.features  <- Get_matrix_features(RNAseq.table,annotation.db)
  RNAseq.features$annotation.db <- annotation.db

  # To be sure the corresponding collumns are converted to their correct type
  RNAseq.table$Bin <- as.character(RNAseq.table$Bin)

  for(sample.col in RNAseq.features$sample.columns){
    RNAseq.table[, sample.col] <- as.numeric(RNAseq.table[, sample.col])
  }

  if (normalize.method != FALSE){
    RNAseq.table <- Normalization(RNAseq.table,
                                  RNAseq.features,
                                  normalize.method)
  }
  if (filter.method    != FALSE){
    RNAseq.table <- Filter(RNAseq.table,
                           RNAseq.features$sample.columns,
                           filter.method)
  }

  RNAseq.table     <- Create.Rank.Columns(RNAseq.table, RNAseq.features)

  if (filter.low.coverage) {
    RNAseq.data <- Filter.Low.Coverage(list("table"    = RNAseq.table,
                                            "features" = RNAseq.features))
  } else {
    RNAseq.data <- list("table"    = RNAseq.table,
                        "features" = RNAseq.features)
  }


  # Combine the table and the features in one object
  return(RNAseq.data)
}

#' Automatically identify features based on column names
#' @name Get_matrix_features
#' @description Robust method to infer information from matrix column names
#' @param RNAseq.table
#' @return A list containing a vector of unique bins, a vector with the sample column
#' numbers, a vector containing the rank column numbers and a annotation_presence_absence matrix
#' @author JJM van Steenbrugge
Get_matrix_features <- function(RNAseq.table,annotation.db){

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
    "annotation_presence_absence" = Get_annotation_presence_absence(RNAseq.table, bins,annotation.db)
  ))
}

#' Calculate which KO is present for each Bin
#' @name Get_annotation_presence_absence
#' @description Calculate which KO is present for each Bin
#' @param RNAseq.table
#' @param bins a vector with all unique bins
#' @author JJM van Steenbrugge
Get_annotation_presence_absence <- function(RNAseq.table, bins,annotation.db){

  a <- unique(annotation.db$`all annotations in a module`)
  b <- unique(as.character(RNAseq.table$Annotation))
  annotation_terms <- unique(c(a,b))

  .Lookup_bin <- function(bin){

    present_bins <- RNAseq.table[which(RNAseq.table$Bin == bin), "Annotation"]

    return(annotation_terms %in% present_bins)
  }

  annotation_presence_absence <- sapply(bins, .Lookup_bin)


  colnames(annotation_presence_absence) <- bins
  rownames(annotation_presence_absence) <- annotation_terms

  return( annotation_presence_absence )
}

#' Normalize RNAseq data
#' @param RNAseq.table
#' @param normalize.method variable that is either a string (default method) or a function to be applied
#' for normalization of RNAseq data.
#' @author JJM van Steenbrugge
Normalization <- function(RNAseq.table, RNAseq.features, normalize.method){
  if(typeof(normalize.method) == "closure"){
    return(normalize.method(RNAseq.table))
  }
  return(RNAseq.table)
}

#' Remove data rows that have a stdev of zero.
#' @param RNAseq.table
#' @param sample.columns a vector with the collumns that contain RNAseq data.
#' @author JJM van Steenbrugge
Filter <- function(RNAseq.table, sample.columns,
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

    return(
      RNAseq.table[which(MADs > 1 ), ] ) #maybe a bigger number than 1

  } else if(typeof(filter.method) == 'closure') {
    filter.method(RNAseq.table)
  }

}
#' Add a ranking collumn for each sample that holds the negative rank of transcripts
#' assigned to each gene, per bin. Normalized on the highest negative rank.
#' @param RNAseq.table
#' @param RNAseq.features
#' @author JJM van Steenbrugge
#' @author BO Oyserman
Create.Rank.Columns <- function(RNAseq.table, RNAseq.features){
  header <- sapply( 1: length(RNAseq.features$rank.columns),
                    function(x) paste("Rank", x, sep = "") )


  for (sample.idx in 1:length(RNAseq.features$sample.columns)) {
    sample.col <- RNAseq.features$sample.columns[sample.idx]
    rank.col   <- rep(NA, length(sample.col))


    for (bin in RNAseq.features$bins) {
      sample.bin.expr <- RNAseq.table[which(RNAseq.table$Bin == bin), sample.col]

      sample.bin.rank <- rank(-sample.bin.expr,
                              na.last     = TRUE,
                              ties.method = "random")

      sample.bin.norm.rank <- sample.bin.rank / max(sample.bin.rank)
      rank.col[which(RNAseq.table$Bin == bin)] <- sample.bin.norm.rank
    }

    RNAseq.table <- cbind(RNAseq.table, rank.col)

  }
  current.cols                               <- colnames(RNAseq.table)
  current.cols[RNAseq.features$rank.columns] <- header
  colnames(RNAseq.table)                     <- current.cols

  return(RNAseq.table)
}

Create.Module.groups <- function (annotation.db) {
  module.dict <- list()
  for(module in unique(annotation.db$Module)){
    module.dict[[module]] <- annotation.db[which(annotation.db$Module == module), 'Annotation']
  }
  return(list("module.dict" = module.dict,
              "all annotations in a module" = unique(unlist(module.dict, use.names = FALSE))
              )
         )


}

Filter.Low.Coverage <-  function (RNAseq.data, threshold = 4) {


  mat <- matrix(ncol=2,nrow=0)
  for(bin in RNAseq.data$features$bins) {
    data        <- RNAseq.data$table[which(RNAseq.data$table$Bin == bin),
                                     RNAseq.data$features$sample.columns]
    total       <- sum(data)
    LCG         <- length(which(data > threshold))
    n.samples   <- length(RNAseq.data$features$sample.columns)
    mat         <- rbind(mat, c((LCG/(nrow(data)* n.samples)), log2(total)))
  }


  bins.keep <- RNAseq.data$features$bins[which(mat[,1] >= 0.8)]

  # Put the trash out
  RNAseq.data$table <- RNAseq.data$table[which(RNAseq.data$table$Bin %in% bins.keep) ,]
  RNAseq.data$features$bins <- bins.keep
  RNAseq.data$features$annotation_presence_absence <- RNAseq.data$features$annotation_presence_absence[, bins.keep]
  return(RNAseq.data)

}
