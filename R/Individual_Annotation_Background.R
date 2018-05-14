#' Individual Random Background distributions
#' @name Individual Random Background distributions
#' @description A function that creates three different random background distributions
#' For each of them, two different genomes are sampled. Then for the first distribution,
#' from each of the two genomes, a random gene is selected. For the second distribution,
#' from each of the two genomes, a random gene that has an annotation is selected.
#' For the third distribution, an annotation is selected that is present in both genomes.
#' Then for each of the two genomes, a gene is selected with that annotation.
#' @author JJM van Steenbrugge
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
                                                   'Random.Identical.Annotated.Genes.bkgd',
                                                   'N',
                                                   'RNAseq.data',
                                                   'metrics')) %dopar% {
    RNAseq.data$annotation.only <- RNAseq.data$table[which(RNAseq.data$table$Annotation != ""),]
    if (i == 1){
      Random.Genes.bkgd(RNAseq.data, metrics, N)
    }else if (i == 2){
      Random.Annotated.Genes.bkgd(RNAseq.data, metrics, N)
    }else if (i == 3){
      Random.Identical.Annotated.Genes.bkgd(RNAseq.data, metrics, N)
    }
  }

  stopCluster(cl)
  names(result) <- c("Random Genes", "Random Annotated Genes", "Genes with the same annotation")
  return(result)
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
#' @param N
Random.Annotated.Genes.bkgd <- function(RNAseq.data, metrics, N){

  # Creating the output format
  output <- list()
  for(metric in names(metrics)){
    output[[metric]] <- rep(NA, N)
  }

    for(i in 1:N){

    random.genomes     <- sample(RNAseq.data$features$bins, 2)

    positions.genome.A <- which(RNAseq.data$annotation.only$Bin == random.genomes[1])

    positions.genome.B <- which(RNAseq.data$annotation.only$Bin == random.genomes[2])

    position.A         <- sample(positions.genome.A, 1)
    position.B         <- sample(positions.genome.B, 1)


    for(metric.current in names(metrics)){
      row.A <- RNAseq.data$annotation.only[position.A, ]
      row.B <- RNAseq.data$annotation.only[position.B, ]

      # Call the distance metric function
      distance <- metrics[[metric.current]](row.A, row.B, RNAseq.data$features)
      output[[metric.current]][i] <- distance
    }

  }

  return(output)
}

#' Background distribution of two random random genes with the same
#' random annotation in two random genomes
#' @author JJM van Steenbrugge
#' @param RNAseq.data
#' @param metrics
#' @export
#' @param N
Random.Identical.Annotated.Genes.bkgd <- function(RNAseq.data, metrics, N, random.genomes){

  #Pre select N pairs of two random genomes each
  if(missing(random.genomes)){
    random.genomes <- rep(list(sample(RNAseq.data$features$bins, 2)), N)
  }

  # Creating the output format
  output <- list()
  for(metric in names(metrics)){
    output[[metric]] <- rep(NA, N)
  }



  for(i in 1:N){

    #random.genomes     <- sample(RNAseq.data$features$bins, 2)

    positions.genome.A <- RNAseq.data$annotation.only[which(RNAseq.data$annotation.only$Bin == random.genomes[[i]][1]), ]

    positions.genome.B <- RNAseq.data$annotation.only[which(RNAseq.data$annotation.only$Bin == random.genomes[[i]][2]), ]

    pool.A <- as.character(positions.genome.A$Annotation)
    pool.B <- as.character(positions.genome.B$Annotation)

    annotations.overlap <- pool.A[pool.A %in% pool.B]
    random.annotation   <- sample(annotations.overlap, 1)

    # Take the first occurence
    position.genome.A <-
      positions.genome.A[which(positions.genome.A$Annotation == random.annotation), ] [1,]

    position.genome.B <-
      positions.genome.B[which(positions.genome.B$Annotation == random.annotation), ] [ 1, ]

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

  return(output)
}


