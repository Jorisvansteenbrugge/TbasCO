#' Prune Trait Attributes
#' @name Prune Trait Attributes
#' @description For each trait attribute, statistical significance is calculated
#' @param filepath RNAseq data file
#' @param distance.metrics Variable
#' @param annotation.db
#' @param bkgd.individual.Zscores
#' @export
#' @author JJM van Steenbrugge
Prune_Trait_Attributes <- function(trait.attributes, bkgd.traits, features,
                                   p.threshold = 0.05, completion.threshold = 0.75,
                                   pairwise.distances){

  annotation.db <- features$annotation.db

  .Filter_Completion <- function(trait.terms, bins, features, completion.threshold){
    n <- length(trait.terms)
    bins_remove <- c()

    for (bin in bins) {
      pa.current <- features$annotation_presence_absence[trait.terms,bin]


      completion <- sum(pa.current) / n

      if (! completion >= completion.threshold) {
        bins_remove <- c(bins_remove, bin)
      }

    }
    return(bins_remove)

  }

  .Calc_P <- function(trait, p.threshold, annotation.db){
    n.genes <- length(annotation.db$module.dict[[trait]])
    trait.attributes.current <- trait.attributes[[trait]]

    cluster.attributes <- list()

    clusters <- levels(as.factor(trait.attributes.current$clusters$membership))

    bins <- trait.attributes.current$cluster$names

    for(cluster in clusters){
      bins.cluster <- bins[which(trait.attributes.current$clusters$membership == cluster)]

      # number of bins <= 2 would end up with only 1 avg zscore which is not enough for a t.test
      if ( length(bins.cluster) <= 2 ) {
        next()
      }

     # Filter for completion
      bins_remove <- .Filter_Completion(annotation.db$module.dict[[trait]],
                                        bins.cluster, features,
                                        completion.threshold )

      if ( length(bins_remove) == length(bins.cluster) ) {
        next()
      }else if (length(bins_remove) == 0) {
        next()
      } else {
        bins.cluster <- bins.cluster[-which(bins.cluster %in% bins_remove)]
      }



      bin.zscores  <- trait.attributes.current$avg.zscore.module[bins.cluster,bins.cluster]

      p.val <- NA
      tryCatch({
        p.val <- t.test(bin.zscores,bkgd.traits[[as.character(n.genes)]], alternative = 'less')$p.value
        p.val <- p.adjust(p.val, method = 'BH', length(clusters))

        if( p.val <= p.threshold){
          cluster.attributes[[cluster]] <- list("genomes"= unique(bins.cluster),
                                                "p-val"  = p.val)
        }
      },
      error = function(cond) {

      })
    }



    return(cluster.attributes)
  }

  trait.names <- names(annotation.db$module.dict)

  trait.attributes.pruned <- list()
  for(i in 1: length(trait.names)) {
    trait.attribute <- .Calc_P(trait.names[i], p.threshold, annotation.db)
    # Backup check if there are not sig attributes
    if (length(trait.attribute) == 0){
      pos.sig <- Identify_Significance_Trait(trait.names[i],
                                             RNAseq.data,
                                             pairwise.distances,
                                             bkgd.traits)
      if (length(pos.sig) > 0) {
        trait.attribute <- list('1' = pos.sig)
      }
    }
    trait.attributes.pruned[[i]] <-  trait.attribute
  }

  names(trait.attributes.pruned) <- trait.names
  return(trait.attributes.pruned)
}

