#' Create an all versus all comparison for each annotation term
#' @name Pairwise comparison
#' @description
#' @param filepath RNAseq data file
#' @param normalize.method Variable
#' @export
#' @author JJM van Steenbrugge
Calc_Pairwise_Annotation_Distance <- function(RNAseq.data, annotation.db,
                                              distance.metrics, bkgd.individual.Zscores){
  if(missing(annotation.db)){
    annotation.db <- RNAseq.data$features$annotation.db
  }

  nbins <- length(RNAseq.data$features$bins)
  #dit voor elke KO term
  ANNOTATION <- 'K03641'

  distances <- list()


  #combined[x,y] <- 1
  for(metric in names(distance.metrics)){
    distances[[metric]] <- matrix(data = NA, nrow = nbins, ncol = nbins)

    for (x in 1:(nbins-1)) {
      genome.A <- RNAseq.data$table[which(RNAseq.data$table$Bin == RNAseq.data$features$bins[x]),]
      genome.A.genes <- genome.A[which(genome.A$Annotation == ANNOTATION),]
      #In case you have multiple, sample 1
      genome.A.gene <- genome.A.genes[sample(nrow(genome.A.genes),size=1),]

      for (y in (x+1):nbins) {
        genome.B <- RNAseq.data$table[which(RNAseq.data$table$Bin == RNAseq.data$features$bins[x]),]
        genome.B.genes <- genome.B[which(genome.B$Annotation == ANNOTATION),]
        #In case you have multiple, sample 1
        genome.B.gene <- genome.B.genes[sample(nrow(genome.B.genes),size=1),]

        distances[[metric]][x,y] <- distance.metrics[[metric]](genome.A.gene,genome.B.gene, RNAseq.data$features)

      }
    }
  }

  distances.Z <- .Convert_zscores(distances, distance.metrics, bkgd.individual.Zscores)

  #Quick and dirty
  composite.distance <- (-distances.Z$PC) + distances.Z$NRED
  # combine them yo
}
