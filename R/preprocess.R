
#' Pre-processing of RNAseq data
#' @description Wrapper of preprocess functions, including: Identifying Matrix
#' features, Normalization, Filtering for non informative lines
#' @param filepath RNAseq data file
#' @param normalize.method User defined method to normalize data
#' @param filter.low.coverage boolean expression to wether we should filter out
#' @param filter.genes.zero False or a number. The number indicates the number of zeros that can max be present. Genes with more zeros than the cut-off get pruned.
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
Pre_process_input <- function(file.path, database, normalize.method = T,
                              filter.method = "MAD",
                              filter.low.coverage = T,
                              filter.genes.zero = F,
                              normalization.features = NULL, taxon_file = NULL){

  cat("Reading table from file..\n")
  RNAseq.table     <- read.csv2(file.path, stringsAsFactors = F)

  RNAseq.table$Annotation[is.na(RNAseq.table$Annotation)] <- ""

  RNAseq.features  <- Get_matrix_features(RNAseq.table, database)
  RNAseq.features$annotation.db <- database

  # To be sure the corresponding collumns are converted to their correct type
  RNAseq.table$Bin <- as.character(RNAseq.table$Bin)

  for(sample.col in RNAseq.features$sample.columns){
    RNAseq.table[, sample.col] <- as.numeric(RNAseq.table[, sample.col])
  }
  if (normalize.method != FALSE){
    cat("Starting Normalization..\n\n")
    RNAseq.table <- Normalization(RNAseq.table,
                                  RNAseq.features,
                                  normalize.method,
                                  normalization.features)
  }else{
    cat("Skipping Normalization..\n")
  }
  if (filter.method    != FALSE){
    cat("Filtering genes..\n")
    RNAseq.table <- Filter(RNAseq.table,
                           RNAseq.features$sample.columns,
                           filter.method)
  }

  if (isTRUE(filter.low.coverage)) { #Defaults
    cat("Filtering low coverage genomes..\n")
    RNAseq.data <- Filter.Low.Coverage(list("table"    = RNAseq.table,
                                            "features" = RNAseq.features))

    #Assume the filter.low.coverage variable is a vector c(threshold, cutoff)
  } else if(! isFALSE(filter.low.coverage)){
    cat("Filtering low coverage genomes (alt)..\n")
    RNAseq.data <- Filter.Low.Coverage(list("table"    = RNAseq.table,
                                            "features" = RNAseq.features))
    RNAseq.data <- Filter.Low.Coverage(list("table"    = RNAseq.table,
                                            "features" = RNAseq.features),
                                       threshold = filter.low.coverage[1],
                                       cutoff = filter.low.coverage[2])

  }
  else {
    cat("Skipping Filter of low coverage genomes..\n")
    RNAseq.data <- list("table"    = RNAseq.table,
                        "features" = RNAseq.features)
  }

  if (! isFALSE(filter.genes.zero)){
    cat("Filtering genes with zero counts..\n")
    RNAseq.data <- Prune.Genes.Zeros(RNAseq.data, filter.genes.zero)
  }else{
    cat("Skipping filter of genes with zero counts..\n")
  }

  RNAseq.table     <- Create.Rank.Columns(RNAseq.table, RNAseq.features)


  #Add trait presence absence later
  #RNAseq.data$features$trait_presence_absence <- Get_trait_presence_absence(RNAseq.data)


  # all_kos <- c(RNAseq.data$features$annotation.db$`all annotations in a module`, unique(RNAseq.data$table$Annotation))
  # expansion <- Expand_module_database(RNAseq.data)
  # RNAseq.data$features$annotation.db          <- expansion$annotation.db
  # RNAseq.data$features$annotation.db$`all annotations in a module` <- all_kos

  cat("Calculating trait presence absence..\n")
   RNAseq.data$features$trait_presence_absence <- Get_trait_presence_absence(RNAseq.data)


  # Combine the table and the features in one object
  return(RNAseq.data)
}


#'@export
Combine_databases <- function(...){
  combined <- c(...)
  all_kos <- combined %>% unlist %>% unique
  return(list('module.dict' = combined,
              'all annotations in a module' = all_kos))
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
  bins <- as.character(sort(unique(RNAseq.table$Bin)))

  annotation_presence_absence <- Get_annotation_presence_absence(RNAseq.table, bins,annotation.db)
  return(list(
    "bins"                        = bins,
    "sample.columns"              = sample.columns,
    "rank.columns"                = ncol(RNAseq.table)+1:length(sample.columns),
    "annotation_presence_absence" = annotation_presence_absence
  ))
}

#' Calculate which KO is present for each Bin
#' @name Get_annotation_presence_absence
#' @description Calculate which KO is present for each Bin
#' @param RNAseq.table
#' @param bins a vector with all unique bins
#' @author JJM van Steenbrugge
Get_annotation_presence_absence <- function(RNAseq.table, bins, annotation.db){
  library(magrittr)



  database_kos <- unlist(annotation.db) %>% unique
  table_kos <- RNAseq.table$Annotation %>% unique

  all_kos <- c(database_kos, table_kos) %>% stringi::stri_remove_empty()

  .Lookup_bin <- function(bin){


    present_bins <- RNAseq.table[which(RNAseq.table$Bin == bin), "Annotation"]

    overlap <- all_kos %in% NULL
    try(
      overlap <- all_kos %in% present_bins
      )
    return(overlap)
  }

  annotation_presence_absence <- sapply(bins, .Lookup_bin)


  colnames(annotation_presence_absence) <- bins
  rownames(annotation_presence_absence) <- all_kos

  return( annotation_presence_absence )
}


Get_trait_presence_absence <- function(RNAseq.data){

  trait_pa <- sapply(RNAseq.data$features$annotation.db$module.dict, Get_trait_pa,
                     RNAseq.data = RNAseq.data, cutoff_fraction = 0.8) %>% rbind
  return(trait_pa)

}

Get_trait_pa <- function(RNAseq.data, module_kos, cutoff_fraction = 0.8) {

  return(sapply(as.character(RNAseq.data$features$bins),function(bin) {

    kos <- RNAseq.data$features$annotation_presence_absence[module_kos, bin]

    true_sum <- kos %>% sum
    if ((true_sum /length(module_kos) >= cutoff_fraction)) {
      return(T)
    }
    return(F)
  }))
}


#' @export
Add_BRITE_DB <- function(RNAseq.data, KEGG_BRITE){
  # Module database
  RNAseq.data$features$annotation.db$module.dict %>% c(., KEGG_BRITE)

  # Ko PA
  RNAseq.data$features$annotation_presence_absence <-  Get_annotation_presence_absence(RNAseq.data$table,
                                                                                       RNAseq.data$features$bins,
                                                                                       RNAseq.data$features$annotation.db)

  # Trait PA
  module_pa <- RNAseq.data$features$trait_presence_absence
  brite_pa <- sapply(KEGG_BRITE,  Get_trait_pa,
                     RNAseq.data = RNAseq.data, cutoff_fraction = 0.66) %>% rbind

  colnames(brite_pa) <- RNAseq.data$features$bins

  RNAseq.data$features$trait_presence_absence <-
      rbind(module_pa, brite_pa[, colnames(module_pa)])


  return(RNAseq.data)

}

#' Normalize RNAseq data
#' @param RNAseq.table
#' @param normalize.method variable that is either a function to be applied or
#' the default method for normalization of RNAseq data.
#' @author JJM van Steenbrugge
Normalization <- function(RNAseq.table, RNAseq.features,
                          normalize.method, normalization.features){
  if(typeof(normalize.method) == "closure"){
    return(normalize.method(RNAseq.table))
  } else if(normalize.method == "simple"){
    return(Normalize(RNAseq.table, RNAseq.features, normalization.features, simple = normalize.method))
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

Filter.Low.Coverage <-  function (RNAseq.data, threshold = 4, cutoff = 0.6) {


  mat <- matrix(ncol = 2, nrow = 0)
  for (bin in RNAseq.data$features$bins) {
    data <- RNAseq.data$table[
      which(RNAseq.data$table$Bin == bin),
      RNAseq.data$features$sample.columns
    ]
    total <- sum(data)
    LCG <- length(which(data >= threshold))
    n.samples <- length(RNAseq.data$features$sample.columns)
    mat <- rbind(mat, c((LCG / (nrow(data) * n.samples)), log2(total)))
  }


  bins.keep <- as.character(RNAseq.data$features$bins[which(mat[,1] >= cutoff)])

  # Put the trash out

  return(Prune_bins(RNAseq.data, bins.keep))

}


Prune.Genes.Zeros <- function(RNAseq.data, threshold = 0){
  rows.keep <- apply(RNAseq.data$table, 1,
                     function(row){
                       sample_cols <- row[RNAseq.data$features$sample.columns] %>% as.numeric
                       c <- length(which(sample_cols == 0))

                       return(c <= threshold)
                       }
                     )
  s <- sum(rows.keep)

  RNAseq.data$table <- RNAseq.data$table[rows.keep,]

  RNAseq.data$features$bins <- RNAseq.data$table$Bin %>% unique
  #
  # min_gene_cutoff <- 1
  # for(bin in RNAseq.data$features$bins){
  #   rows <- RNAseq.data$table[which(RNAseq.data$table$Bin==bin),] %>% nrow
  #   if (rows >= min_gene_cutoff){
  #     bins.keep <- c(bins.keep, bin)
  #   }
  # }

  ## Check if we pruned all genes from a genome
  return(RNAseq.data)

}

Prune_bins <- function(RNAseq.data, bins.keep){
  RNAseq.data$table <- RNAseq.data$table[which(RNAseq.data$table$Bin %in% bins.keep) ,]
  RNAseq.data$features$bins <- bins.keep
  RNAseq.data$features$annotation_presence_absence <- RNAseq.data$features$annotation_presence_absence[, bins.keep]
  RNAseq.data$features$trait_presence_absence      <- RNAseq.data$features$trait_presence_absence[, bins.keep]

  return(RNAseq.data)
}




#' Normalized RNAseq raw read counts
#'
#'
#' RNAseq raw read counts may by normalized based on various parameters including reads per sample,
#' reads mapped per genome, gene length etc. Here we implement the normalization of edgeR (citation)
#' which accounts for differences in both Sequencing Detph and RNA composition (see edgeR documentation
#' page 2.7.2 & 2.7.3). However, in metatranscritpomic studies, it may also be beneficial to
#' adjust for an additional source of compositional bias in which a single organisms may contribute
#' a high relative abundance of transcrtipts, resulting in an undersampling of other organisms.
#' Therefore, we provide an additional normalization step to normalize on a per genome/bin basis.
#'
#' @param RNAseq_Annotated_Matrix The original count matrix (See X for format details).
#' @param no_feature,ambiguous,not_aligned  A set of vectors equal to the number of samples,
#' containing the number of reads that had no feature,
#' where ambiguously mapped, or not aligned in their  (obtained from the mapping output).
#' @param gene_lengths A matrix with the length of each gene (genes must be in same order as input RNAseq_Annotated_Matrix)
#' @param method A string containing the method to use, either one of: ["default", "TMM", "RLE"].  In addition to the described default method, TMM and RLE from bioconductors edgeR
#' package are implemented as well
#' @export
#' @author BO Oyserman
#' @return The normalized read counts  of \code{Sample 1} ... \code{Sample N}.
#' @examples RNAseq_Normalize(RNAseq_Annotated_Matrix, no_feature,ambiguous, not_aligned)
#' @note \preformatted{To remove rows that have a 0 for its read counts:}
#' \code{RNAseq_Annotated_Matrix[apply(RNAseq_Annotated_Matrix[, SS:SE], 1, function(x) !any(x == 0)), ]}
#' \preformatted{Where SS and SE are the start and end columns of the samples (raw counts).}
Normalize <- function(RNAseq.table, RNAseq.features, normalization.features, simple){

  # normalized by total of non-rRNA reads per sample mapped
  sum_aligned         <- apply(RNAseq.table[, RNAseq.features$sample.columns],
                               2, sum)

  if (missing(simple)){
  total_nonRNA_reads  <- sum_aligned + normalization.features$no_feature
  total_nonRNA_reads  <-  total_nonRNA_reads + normalization.features$ambiguous
  total_nonRNA_reads  <-  total_nonRNA_reads  + normalization.features$not_aligned


  normalized_by_total <- total_nonRNA_reads / max(total_nonRNA_reads)
  } else{
  normalized_by_total <- normalization.features / max(normalization.features)
 }
  # An alternative would be to use the library size.
  # However, since the pool of rRNA is often both physically and bioinformatically removed,
  # the count of mRNA reads per sample is more intuitively relevant.
  # normalized_by_total <- library_size/max(library_size)


  RNAseq.table[,RNAseq.features$sample.columns] <- t(t(RNAseq.table[,RNAseq.features$sample.columns])
                                                     /normalized_by_total)

  # convert to log base 2
  # RNAseq_Annotated_Matrix[,SS:SE]<-log(RNAseq_Annotated_Matrix[,SS:SE],2)

  # replace -Inf with 0
  new_cols <- apply(RNAseq.table[, RNAseq.features$sample.columns], 2,
                        function(col) {
                          col[is.infinite(col)] <- 0
                          return(col)
                        })
  RNAseq.table[, RNAseq.features$sample.columns] <- new_cols

  return(Normalize_by_bin(RNAseq.table, RNAseq.features))



}



#' @author BO Oyserman
Normalize_by_bin <- function(RNAseq.table, RNAseq.features){



  # Step 1: Calculate the number of reads mapped to each bin in each sample (This may be a separate function)
  sum_reads_per_genome_matrix<-matrix(NA,
                                      nrow = length(RNAseq.features$bins),
                                      ncol = length(RNAseq.features$sample.columns))

  for (i in 1: length(RNAseq.features$bins) ){

    for (j in 2: (length(RNAseq.features$sample.columns) + 1) ) {

      bin_rows <- which(RNAseq.table[, "Bin"] == RNAseq.features$bins[i])
      sum_reads_per_genome_matrix[i, j - 1] <- sum(RNAseq.table[bin_rows, j])
    }
  }

  # Step 2: Calculate max per bin.
  # Divide each column (sample) per row in normalized_sum_reads_per_genome_matrix (each bin) by
  # the max count per bin

  normalized_sum_reads_per_genome_matrix<-sum_reads_per_genome_matrix/
    apply(sum_reads_per_genome_matrix,1,max)

  # Step 3: normalize reads by max mapped to a genome
  for (i in 1: length(RNAseq.features$bins)) {
    bin_rows <- which(RNAseq.table[, 'Bin'] == RNAseq.features$bins[i])

    RNAseq.table[bin_rows, RNAseq.features$sample.columns] <- t(t(RNAseq.table[bin_rows, RNAseq.features$sample.columns])/
                                                    normalized_sum_reads_per_genome_matrix[i, ])

  }

  return(RNAseq.table)
}

