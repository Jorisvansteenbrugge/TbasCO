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
    pa.current <- features$annotation_presence_absence[trait.terms,bins]

    if(n == 1){
      completion <- sum(pa.current) / n
    }else{
      complete <- apply(pa.current,1,function(x){if(sum(x)>=1){return(1)}else{return(0)}})
      completion <- sum(complete) / n
    }
    return(completion >= completion.threshold)

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
      if(! .Filter_Completion(annotation.db$module.dict[[trait]], bins, features, completion.threshold )){
        next()
      }



      bin.zscores  <- trait.attributes.current$avg.zscore.module[bins.cluster,bins.cluster]

      p.val <- t.test(bin.zscores,bkgd.traits[[as.character(n.genes)]], alternative = 'less')$p.value
      p.val <- p.adjust(p.val, method = 'BH', length(clusters))

      if( p.val <= p.threshold){
        cluster.attributes[[cluster]] <- list("genomes"= bins.cluster,
                                              "p-val"  = p.val)
      }
    }



    return(cluster.attributes)
  }

  trait.names <- names(annotation.db$module.dict)

  trait.attributes.pruned <- list()
  for(i in 1: length(trait.names)) {
    trait.attribute <- .Calc_P(trait.names[i], p.threshold, annotation.db)
    # Backup check if there are not sig attributes
    if (length(trait.attribute) == 0){
    # try(
      pos.sig <- Identify_Significance_Trait(trait.names[i],
                                             RNAseq.data,
                                             pairwise.distances,
                                             bkgd.traits)
      if (length(pos.sig) > 0) {
        trait.attribute <- list('1' = pos.sig)
      }

    # ,silent = T)
    #
    }

    trait.attributes.pruned[[i]] <-  trait.attribute
  }

  names(trait.attributes.pruned) <- trait.names
  return(trait.attributes.pruned)
}

