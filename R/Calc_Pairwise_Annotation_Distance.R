#' Create an all versus all comparison for each annotation term
#' @name Pairwise comparison of Annotations
#' @description For each annotation term, a comparison is done for each genome
#' versus all other genomes, 1 versus 1.
#' @param filepath RNAseq data file
#' @param distance.metrics Variable
#' @param annotation.db
#' @param bkgd.individual.Zscores
#' @export
#' @author JJM van Steenbrugge
Calc_Pairwise_Annotation_Distance <- function(RNAseq.data, annotation.db,
                                              distance.metrics, bkgd.individual.Zscores,
                                              show.progress = F){
  if(missing(annotation.db)){
    annotation.db <- RNAseq.data$features$annotation.db
  }


  #This is a nested function to make the foreach loop below more straightforward
  .Pairwise_Distance <- function(annotation_term){
    distances <- list()


    #combined[x,y] <- 1
    for(metric in names(distance.metrics)){
      distances[[metric]] <- matrix(data = NA, nrow = nbins, ncol = nbins)

      for (x in 1: (nbins - 1)) {

        genome.A       <- RNAseq.data$table[which(RNAseq.data$table$Bin == RNAseq.data$features$bins[x]),]
        genome.A.genes <- genome.A[which(genome.A$Annotation == annotation_term),]

        #In case you have multiple, sample 1
        genome.A.gene  <- genome.A.genes[sample(nrow(genome.A.genes),size=1),]

        for (y in (x + 1): nbins) {

          genome.B       <- RNAseq.data$table[which(RNAseq.data$table$Bin == RNAseq.data$features$bins[x]),]
          genome.B.genes <- genome.B[which(genome.B$Annotation == annotation_term),]

          #In case you have multiple, sample 1
          genome.B.gene  <- genome.B.genes[sample(nrow(genome.B.genes),size=1),]

          distances[[metric]][x, y] <- distance.metrics[[metric]](genome.A.gene,genome.B.gene, RNAseq.data$features)
        }
      }
    }

    distances.Z <- .Convert_zscores(distances, distance.metrics, bkgd.individual.Zscores)

    #Quick and dirty
    return( (-distances.Z$PC) + distances.Z$NRED)
  }

  nbins <- length(RNAseq.data$features$bins)



  require(doSNOW)
  cl <- makeSOCKcluster(4)
  registerDoSNOW(cl)




  #Runs in parallel
  composite.distances <- foreach::foreach(i = 1: length(annotation.db$`all annotations in a module`),
                                          .export = c('.Pairwise_Distance',
                                                      'distance.metrics',
                                                      'nbins',
                                                      'RNAseq.data',
                                                      'annotation.db',
                                                      '.Convert_zscores',
                                                      'bkgd.individual.Zscores'),
                                          .verbose = show.progress) %do%{

      return(.Pairwise_Distance(annotation.db$`all annotations in a module`[i]))
  }

  stopCluster(cl)
  names(composite.distances) <- annotation.db$`all annotations in a module`
  return(composite.distances)
}
