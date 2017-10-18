#' Individual Random Background distributions
Individual_Annotation_Background <- function(RNAseq.data,
                                             N       = 1000,
                                             cores   = 2,
                                             metrics = list("PC"   = PC,
                                                            "NRED" = NRED)){
  require(doSNOW)
  cl <- makeSOCKcluster(3)
  registerDoSNOW(cl)


  result <- foreach::foreach(i = 1: 3, .export = c('Random.Genes.bkgd',
                                                   'Random.Annotated.Genes.bkgd',
                                                   'Random.Identical.Annotated.Genes.bkgd',
                                                   'N',
                                                   'RNAseq.data',
                                                   'metrics')) %dopar% {
    if(i == 1){
      Random.Genes.bkgd(RNAseq.data, metrics, N)
    }else if(i == 2){
      Random.Annotated.Genes.bkgd(RNAseq.data, metrics, N)
    }else if (i == 3){
      Random.Identical.Annotated.Genes.bkgd(RNAseq.data, metrics, N)
    }
  }

  stopCluster(cl)
  names(result) <- c("Random Genes", "Random Annotated Genes", "Genes with the same annotation")
  return(result)
}



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
      print(metric.current)
      row.A <- RNAseq.data$table[position.A, ]
      row.B <- RNAseq.data$table[position.B, ]

      # Call the distance metric function
      distance <- metrics[[metric.current]](row.A, row.B, RNAseq.data$features)
      print(distance)
      output[[metric.current]][i] <- distance
    }
  }

  return(output)
}

Random.Annotated.Genes.bkgd <- function(RNAseq.data, metrics, N){

  # Creating the output format
  output <- list()
  for(metric in names(metrics)){
    output[[metric]] <- rep(NA, N)
  }


    for(i in 1:N){

    random.genomes     <- sample(RNAseq.data$features$bins, 2)

    positions.genome.A <- which(RNAseq.data$table$Bin == random.genomes[1] &
                                RNAseq.data$table$Annotation != "")

    positions.genome.B <- which(RNAseq.data$table$Bin == random.genomes[2] &
                                RNAseq.data$table$Annotation != "")

    position.A         <- sample(positions.genome.A, 1)
    position.B         <- sample(positions.genome.B, 1)


    for(metric.current in names(metrics)){
      print(metric.current)
      row.A <- RNAseq.data$table[position.A, ]
      row.B <- RNAseq.data$table[position.B, ]

      # Call the distance metric function
      distance <- metrics[[metric.current]](row.A, row.B, RNAseq.data$features)
      print(distance)
      output[[metric.current]][i] <- distance
    }

  }

  return(output)
}

Random.Identical.Annotated.Genes.bkgd <- function(RNAseq.data, metrics, N){

  # Creating the output format
  output <- list()
  for(metric in names(metrics)){
    output[[metric]] <- rep(NA, N)
  }


  for(i in 1:N){

    random.genomes     <- sample(RNAseq.data$features$bins, 2)

    positions.genome.A <- RNAseq.data$table[which(RNAseq.data$table$Bin == random.genomes[1] &
                                  RNAseq.data$table$Annotation != ""), ]

    positions.genome.B <- RNAseq.data$table[which(RNAseq.data$table$Bin == random.genomes[2] &
                                  RNAseq.data$table$Annotation != ""), ]

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
      print(metric.current)
      row.A <- position.A
      row.B <- position.B

      # Call the distance metric function
      distance <- metrics[[metric.current]](row.A, row.B, RNAseq.data$features)
      print(distance)
      output[[metric.current]][i] <- distance
    }

  }

  return(output)
}



# Example metrics
PC <- function(rowA, rowB, RNAseq.features){
  return(cor(as.numeric(rowA[RNAseq.features$sample.columns]),
             as.numeric(rowB[RNAseq.features$sample.columns])
             )
         )
}

NRED <- function(rowA, rowB, RNAseq.features) {
  r.A <- as.numeric(rowA[RNAseq.features$rank.columns])
  r.B <- as.numeric(rowB[RNAseq.features$rank.columns])
  return(
    sum((r.A - r.B) * (r.A - r.B))
  )
}
