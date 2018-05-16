#' Identify Trait Attributes
#' @name Identify Trait Attributes
#' @description For each trait (or module) Attributes are calculated
#' @param filepath RNAseq data file
#' @param distance.metrics Variable
#' @param annotation.db
#' @param bkgd.individual.Zscores
#' @export
#' @author JJM van Steenbrugge
Identify_Trait_Attributes <- function(RNAseq_Annotated_Matrix,pairwise.distances,
                                      annotation.db, threads = 4){
  require(doSNOW)
  cl <- makeSOCKcluster(threads)
  registerDoSNOW(cl)

  # Allows for easy subset testing
  if(missing(annotation.db)){
    annotation.db <- RNAseq.data$features$annotation.db
  }

  trait.names <- names(annotation.db$module.dict)

  foreach::foreach(i = 1: length(trait.names),
                          .export = c(),
                          .verbose = show.progress) %dopar%{

      module.terms             <- annotation.db$module.dict[[i]]

      jaccard.module.distance  <- .Calc_Jaccard_Module(RNAseq.data,
                                                      module.terms)

      avg.zscore.module        <- .Calc_Avg_Zscore_Module(RNAseq.data,
                                                         pairwise.distances,
                                                         module.terms)
      pairwise.module.distance <- (1-jaccard.module.distance)
      jpe.distance             <- avg.zscore.module * pairwise.module.distance
  }



  stopCluster(cl)

}


.Calc_Jaccard_Module <- function(RNAseq.data, module.terms){

  nbins <- length(RNAseq.data$features$bins)

  jaccard.module.distance <- matrix(data = NA,
                                    nrow = nbins, ncol = nbins,
                                    dimnames = list(RNAseq.data$features$bins,
                                                    RNAseq.data$features$bins))

  for (x in 1: (nbins - 1)) {
    bin.a <- RNAseq.data$features$bins[x]
    for (y in (x + 1): nbins) {
      bin.b <- RNAseq.data$features$bins[y]
      jaccard.module.distance[x,y]  <- .Calc_Jaccard(RNAseq.data,
                                                     c(bin.a, bin.b),
                                                     module.terms)
    }
  }
  return(jaccard.module.distance)
}

.Calc_Avg_Zscore_Module <- function(RNAseq.data, pairwise.distances,
                                    module.terms){

  nbins <- length(RNAseq.data$features$bins)

  avg.zscore.module <- matrix(data = NA,
                              nrow = nbins, ncol = nbins,
                              dimnames = list(RNAseq.data$features$bins,
                                              RNAseq.data$features$bins))

  for (x in 1: (nbins - 1)) {
    bin.a <- RNAseq.data$features$bins[x]
    for (y in (x + 1): nbins) {
      bin.b <- RNAseq.data$features$bins[y]

      avg.zscore.module[x,y] <- mean(unlist(pairwise.distances[module.terms]),
                                     na.rm = TRUE)

    }
  }

  return(avg.zscore.module)
}


