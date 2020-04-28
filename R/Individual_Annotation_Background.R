#' Individual Random Background distributions
#' @name Individual Random Background distributions
#' @description A function that creates three different random background distributions
#' For each of them, two different genomes are sampled. Then for the first distribution,
#' from each of the two genomes, a random gene is selected. For the second distribution,
#' from each of the two genomes, a random gene that has an annotation is selected.
#' For the third distribution, an annotation is selected that is present in both genomes.
#' Then for each of the two genomes, a gene is selected with that annotation.
#' @author JJM van Steenbrugge
#' BO Oyserman
#' @param RNAseq.data
#' @param N Number of iterations
#' @param threads Number of cpu cores to be used
#' @param metrics Named list containing possibly multiple functions for distance
#' metrics
#' @export
Individual_Annotation_Background <- function(RNAseq.data,
                                             N       = 1000,
                                             threads = 3,
                                             metrics = list("PC"   = PC,
                                                            "NRED" = NRED)){
  require(doSNOW)
  cl <- makeSOCKcluster(threads)
  registerDoSNOW(cl)


  result <- foreach::foreach(i = 1: 3, .export = c('Random.Genes.bkgd',
                                                   'Random.Annotated.Genes.bkgd',
                                                   'Random.Identical.Annotated.Genes.bkgd'#,
                                                   # 'N',
                                                   # 'RNAseq.data',
                                                   # 'metrics'
                                                   )) %dopar% {
    RNAseq.data$annotation.only <- RNAseq.data$table[which(RNAseq.data$table$Annotation != ""),]
    if (i == 1){
      Random.Genes.bkgd(RNAseq.data, metrics, N)
    }else if (i == 2){
      Random.Annotated.Genes.bkgd.new(RNAseq.data, metrics, N)
    }else if (i == 3){
      Random.Identical.Annotated.Genes.bkgd.new(RNAseq.data, metrics, N)
    }
  }


  stopCluster(cl)
  # The names are essential
  names(result) <- c("Random Genes", "Random Annotated Genes", "Genes with the same annotation")
  return(result)
}

Generate_empty_output <- function(metrics, N){
  # Creating the output format
  output <- list()
  for(metric in names(metrics)){
    output[[metric]] <- rep(NA, N)
  }

  return(output)
}

Random.Annotated.Genes.bkgd.new <- function(RNAseq.data, metrics, N){

  output <- Generate_empty_output(metrics, N)
  rows.c <- nrow(RNAseq.data$annotation.only)

  for(i in 1:N){
    # Select two rows
    selected_rows <- sample.int(n=rows.c, size=2)

    for(metric.current in names(metrics)){
      row.A <- RNAseq.data$annotation.only[selected_rows[1], ]
      row.B <- RNAseq.data$annotation.only[selected_rows[2], ]

      distance <- metrics[[metric.current]](row.A, row.B, RNAseq.data$features)
      output[[metric.current]][i] <- distance
    }
  }

  return(output)
}


Get_non_unique_annotations <- function(RNAseq.data){
  counts <- RNAseq.data$annotation.only$Annotation %>% table
  non_uniq.counts <- which(counts>=2) %>% names

  return(non_uniq.counts)
}

Sample_random_genes_identical_annotation <- function(RNAseq.data, duplicate_annotations){
  random_gene <- sample(duplicate_annotations, 1)

  tab <- RNAseq.data$table[which(RNAseq.data$table$Annotation == random_gene),]
  tab.nrow <- nrow(tab)

  rows.selected <- sample.int(tab.nrow, 2)
  return(list("rowA" = tab[rows.selected[1], ],
               "rowB" = tab[rows.selected[2], ])
  )
}

Random.Identical.Annotated.Genes.bkgd.new <- function(RNAseq.data, metrics, N){
  output <- Generate_empty_output(metrics, N)
  duplicate_annotations <- Get_non_unique_annotations(RNAseq.data)

  for(i in 1:N){
    rows <- Sample_random_genes_identical_annotation(RNAseq.data, duplicate_annotations)
    for(metric.current in names(metrics)){
      row.A <- rows$rowA
      row.B <- rows$rowB

      distance <- metrics[[metric.current]](row.A, row.B, RNAseq.data$features)
      output[[metric.current]][i] <- distance
    }

  }

  return(output)
}

#' Background distribution of any two random genes in two random genomes
#' @author JJM van Steenbrugge
#' @param RNAseq.data
#' @param metrics
#' @param N
Random.Genes.bkgd <- function(RNAseq.data, metrics, N){

  # Creating the output format
  output <- list()
  for(metric in names(metrics)){
    output[[metric]] <- rep(NA, N)
  }


  for(i in 1:N){

    random.genomes     <- sample(RNAseq.data$features$bins, 2)

    positions.genome.A <- which(RNAseq.data$table$Bin == random.genomes[1])
    positions.genome.B <- which(RNAseq.data$table$Bin == random.genomes[2])

    position.A         <- sample(positions.genome.A, 1)
    position.B         <- sample(positions.genome.B, 1)

    for(metric.current in names(metrics)){
      row.A <- RNAseq.data$table[position.A, ]
      row.B <- RNAseq.data$table[position.B, ]

      # Call the distance metric function
      distance <- metrics[[metric.current]](row.A, row.B, RNAseq.data$features)
      output[[metric.current]][i] <- distance
    }
  }

  return(output)
}

#' Background distribution of any two random genes with an annotation in two
#' random genomes
#' @author JJM van Steenbrugge
#' @param RNAseq.data
#' @param metrics
#' @param random.genomes Allows reuse of this function in the module background function @see@seealso \code{\link{Random Background distributions of modules}}
#' @param N
#' @export
Random.Annotated.Genes.bkgd <- function(RNAseq.data, metrics, N, random.genomes){
  library(magrittr)
  out.terms = T
  RNAseq.data$annotation.only <- RNAseq.data$table[which(RNAseq.data$table$Annotation != ""),]
  #Pre select N pairs of two random genomes each
  if( missing(random.genomes) ) {
    random.genomes <- lapply(1:N, function(x) sample(RNAseq.data$features$bins, 2))
    out.terms = F
  }

  # Creating the output format
  output <- Generate_empty_output(metrics, N)

  #format used.terms
  used.terms <- list()
  unique_random.genomes <- random.genomes %>% unlist %>% unique
  for(genome in unique_random.genomes){
    used.terms[[as.character(genome)]] <- NA
  }

  for(i in 1:N){
    print(i)
    positions.genome.A <- which(RNAseq.data$annotation.only$Bin == random.genomes[[i]][1])

    positions.genome.B <- which(RNAseq.data$annotation.only$Bin == random.genomes[[i]][2])

    position.A         <- sample(positions.genome.A, 1)
    position.B         <- sample(positions.genome.B, 1)

    A.char <- as.character(random.genomes[[i]][1])
    B.char <- as.character(random.genomes[[i]][2])
    used.terms[[A.char]] <- c(used.terms[[A.char]],
                                              as.character(RNAseq.data$annotation.only[position.A,]$Annotation))
    used.terms[[B.char]] <- c(used.terms[[B.char]],
                                              as.character(RNAseq.data$annotation.only[position.B,]$Annotation))

   for(metric.current in names(metrics)){
     row.A <- RNAseq.data$annotation.only[position.A, ]
     row.B <- RNAseq.data$annotation.only[position.B, ]

      # Call the distance metric function
     distance <- metrics[[metric.current]](row.A, row.B, RNAseq.data$features)
     output[[metric.current]][i] <- distance
   }

  }
  if(out.terms){
    return(list("used terms" = used.terms,
                "scores"  = output))
  }else{
    return(output)
  }
}

#' Background distribution of two random random genes with the same
#' random annotation in two random genomes.
#' @name Random.Identical.Annotated.Genes.bkgd
#' @export
#' @author JJM van Steenbrugge
#' @param RNAseq.data
#' @param metrics
#' @param N
#' @param random.genomes Allows reuse of this function in the module background function @see@seealso \code{\link{Random Background distributions of modules}}
Random.Identical.Annotated.Genes.bkgd <- function(RNAseq.data, metrics, N, random.genomes){
  library(magrittr)
  #Pre select N pairs of two random genomes each
  if( missing(random.genomes)) {
    random.genomes <- lapply(1:N, function(x) sample(RNAseq.data$features$bins, 2))
    out.terms = F
  }

  # Creating the output format
  output <- list()
  for( metric in names(metrics)) {
    output[[metric]] <- rep(NA, N)
  }

  used.terms <- list()
  unique_random.genomes <- random.genomes %>% unlist %>% unique
  for(genome in unique_random.genomes){
    used.terms[[as.character(genome)]] <- NA
  }

  for( i in 1:N) {
    positions.genome.A <- RNAseq.data$annotation.only[which(RNAseq.data$annotation.only$Bin == random.genomes[[i]][1]), ]
    positions.genome.B <- RNAseq.data$annotation.only[which(RNAseq.data$annotation.only$Bin == random.genomes[[i]][2]), ]

    pool.A <- as.character(positions.genome.A$Annotation)
    pool.B <- as.character(positions.genome.B$Annotation)

    annotations.overlap <- pool.A[pool.A %in% pool.B]
    random.annotation   <- sample(annotations.overlap, 1)

    A.char <- as.character(random.genomes[[i]][1])
    B.char <- as.character(random.genomes[[i]][2])
    used.terms[[A.char]] <- c(used.terms[[A.char]],random.annotation)
    used.terms[[B.char]] <- c(used.terms[[B.char]],random.annotation)

    # Take the first occurence
    position.genome.A <-
      positions.genome.A[which(positions.genome.A$Annotation == random.annotation), ] [1, ]

    position.genome.B <-
      positions.genome.B[which(positions.genome.B$Annotation == random.annotation), ] [1, ]

    # Tmp fix
    position.A         <- position.genome.A
    position.B         <- position.genome.B


    for(metric.current in names(metrics)){

      row.A <- position.A
      row.B <- position.B

      # Call the distance metric function
      distance <- metrics[[metric.current]](row.A, row.B, RNAseq.data$features)

      output[[metric.current]][i] <- distance
    }

  }

  return(list("used terms" = used.terms,
              "scores"  = output))

}


