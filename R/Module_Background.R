#' Generate a background distibution of composite Z scores
#' @name Random Background distributions of modules
#' @author JJM van Steenbrugge
#' @param RNAseq.data
#' @param N Number of genes to include in the random module
#' @param Z Number of iterations
#' @param range NOT SUPPORTED YET
#' @param metrics Named list containing possibly multiple functions for distance
#' @param threads Number of cpu cores to be used
#' @export
Module_Background <- function(RNAseq.data,
                              Z_scores,
                              N,
                              Z,
                              metrics,
                              threads = 2){
  require(doSNOW)
  cl <- makeSOCKcluster(threads)
  registerDoSNOW(cl)



  result <- foreach::foreach(i = 1: N) %dopar% {
    #pick two genomes
    RNAseq.data$annotation.only <- RNAseq.data$table[which(RNAseq.data$table$Annotation != ""),]
    random.genomes <- sample(RNAseq.data$features$bins, 2)
    random.genomes <- rep(list(random.genomes), Z) #Workaround to re-use function below

    distances <- Random.Identical.Annotated.Genes.bkgd(RNAseq.data, metrics, Z, random.genomes)
    # Get Z scores

    # Jaccard distance

    # This is quick and dirty
    # composite.distance <- mean(distances$PC)
    #   mean((-Random_Zscore_Pearson_Distances) +
    #        Random_Zscore_Euclidean_Distances,
    #      na.rm = TRUE)[1]


  }
}
