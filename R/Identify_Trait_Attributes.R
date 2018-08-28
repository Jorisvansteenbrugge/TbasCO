#' Identify Trait Attributes
#' @name Identify Trait Attributes
#' @description For each trait (or module) Attributes are calculated
#' @param RNAseq.data Collection of multple components, include RNA seq data, annotations, etc.
#' See \code{\link{Pre_process_input}} for the full list.
#' @param pairwise.distances Named list containing possibly multiple functions for distance
#' @param annotation.db List containing a dictionary like structure with trait as names
#' and annotation as values.
#' @param bkgd.individual.Zscores
#' @export
#' @author JJM van Steenbrugge
#' BO Oyserman
Identify_Trait_Attributes <- function(RNAseq.data,pairwise.distances,
                                      annotation.db, threads = 4){
  require(doSNOW)

  cl <- snow::makeSOCKcluster(threads)
  registerDoSNOW(cl)

  # Allows for easy subset testing
  if(missing(annotation.db)){
    annotation.db <- RNAseq.data$features$annotation.db
  }




  trait.names <- names(annotation.db$module.dict)

  output.list <- foreach::foreach(i = 1: length(trait.names),
                          .export = c('.Calc_Jaccard_Module',
                                      '.Calc_Jaccard',
                                      '.Calc_Avg_Zscore_Module'),
                          .verbose = F) %dopar%{
      # In case a bin has no presence for the trait

      module.terms             <- annotation.db$module.dict[[i]]

      jaccard.module.distance  <- .Calc_Jaccard_Module(RNAseq.data,
                                                      module.terms)

      avg.zscore.module        <- .Calc_Avg_Zscore_Module(RNAseq.data,
                                                         pairwise.distances,
                                                         module.terms)
      pairwise.module.distance <- (1-jaccard.module.distance)
      module.distances         <- avg.zscore.module * pairwise.module.distance
      module.distances.table   <- as.data.frame(subset(reshape2::melt(module.distances),
                                         value != 0))
      colnames(module.distances.table)      <- c("Genome_1","Genome_2","Distance")
      module.distances.table.reduced        <- module.distances.table[which(
                                                module.distances.table[, 3] < 0), ]
      module.distances.table.reduced[, 3]   <- (-module.distances.table.reduced[, 3])

      graph                    <- igraph::graph_from_data_frame(module.distances.table.reduced,
                                                                 directed = FALSE)
      clusters                 <- igraph::cluster_louvain(graph,
                                                          weights = igraph::E(graph)$Distance)

      return( list("avg.zscore.module"        = avg.zscore.module,
                   "module.distances"         = module.distances,
                   "pairwise.module.distance" = pairwise.module.distance,
                   "JPE_distance_Table"       = module.distances.table,
                   "clusters"                 = clusters)
      )
                          }

  stopCluster(cl)
  names(output.list) <- trait.names
  return(output.list)
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

#' @export
.Calc_Avg_Zscore_Module <- function(RNAseq.data, pairwise.distances,
                                    module.terms){

  nbins <- length(RNAseq.data$features$bins)

  avg.zscore.module <- matrix(data = NA,
                              nrow = nbins, ncol = nbins,
                              dimnames = list(RNAseq.data$features$bins,
                                              RNAseq.data$features$bins))

  for (x in 1: (nbins - 1)) {
    for (y in (x + 1): nbins) {

        d <- mean(sapply(module.terms, function(term) pairwise.distances[[term]][x,y] ),
                                     na.rm = T)
        avg.zscore.module[x,y] <- d
    }
  }

  return(avg.zscore.module)
}


