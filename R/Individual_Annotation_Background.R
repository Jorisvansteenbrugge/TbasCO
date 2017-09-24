#' Individual Random Background distributions
Individual_Annotation_Background <- function(RNAseq.data,
                                             N       = 1000,
                                             cores   = 2,
                                             metrics = list("PC"   = PC,
                                                            "NRED" = NRED)){
  require(doSNOW)
  cl <- makeSOCKcluster(3)
  registerDoSNOW(cl)


  result <- foreach::foreach(i = 1: 3, .export = c('bg1',
                                                   'bg2',
                                                   'bg3',
                                                   'N',
                                                   'RNAseq.data',
                                                   'metrics')) %dopar% {
    if(i == 1){
      bg1(RNAseq.data, metrics)
    }else if(i == 2){
      bg2()
    }else{
      bg3()
    }
                             }

  stopCluster(cl)

  return(result)
}



bg1 <- function(RNAseq.data, metrics, N){

  for(i in 1:N){

    random.genomes <- sample(RNAseq.data$features$bins, 2)

    positions.genome.A <- which(RNAseq.data$table$Bin == random.genomes[1])
    positions.genome.B <- which(RNAseq.data$table$Bin == random.genomes[2])

    #sample one gene from both
    position.A <- sample(positions.genome.A, 1)
    position.B <- sample(positions.genome.B, 1)


  }

}

bg2 <- function(){
  return(rep("2", 10))
}

bg3 <- function(){
  return(rep("3", 10))
}



# Example metrics
PC <- function(rowA, rowB, RNAseq.features){
  return(cor(rowA[RNAseq.features$sample.columns],
             rowB[RNAseq.features$sample.columns]
             )
         )
}

NRED <- function(rowA, rowB, RNAseq.features) {
  sum((rowA[RNAseq.features$rank.columns] - rowB[RNAseq.features$rank.columns])
      * (rowA[RNAseq.features$rank.columns] - rowB[RNAseq.features$rank.columns]))
}
