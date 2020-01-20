#' Prune Trait Attributes
#' @name Prune Trait Attributes
#' @description For each trait attribute, statistical significance is calculated
#' @param filepath RNAseq data file
#' @param distance.metrics Variable
#' @param annotation.db
#' @param bkgd.individual.Zscores
#' @export
#' @author JJM van Steenbrugge
Prune_Trait_Attributes <- function(trait.attributes, bkgd.traits, RNAseq.data,
                                   p.threshold = 0.05, pairwise.distances, bkgd.individual.Zscores,
                                   annotation.db, trait_presence_absence){

  features <- RNAseq.data$features

  if (missing(annotation.db)){
    annotation.db <- features$annotation.db
  }
  if (missing(trait_presence_absence)) {
    trait_presence_absence <- features$trait_presence_absence
  }


  .Filter_Completion <- function(features, trait, bins){

    print(trait)
    print(bins)
    bin.completions <- as.logical(trait_presence_absence[as.character(bins)], trait)
    return(features$bins[!bin.completions])

  }

  .Calc_P_pairwise <- function(bins, trait){
    kos <- RNAseq.data$features$annotation.db$module.dict[[trait]]
    distances <- kos %>% sapply(function(ko) {
      dists <- pairwise.distances[[ko]] [bins, bins]
      return(dists[ which(!is.na(dists)) ])
    }) %>% unlist %>% as.numeric

    zscores <-  bkgd.individual.Zscores$zscores$`Genes with the same annotation`

    bkgd <- (-zscores$PC) + zscores$NRED


    pval <- 30
    try(pval <- t.test(distances, bkgd,
               alternative = 'less')$p.value
        ,silent = T)

    return(
      pval
    )


  }

  .Calc_P <- function(trait, p.threshold, annotation.db){
    n.genes <- length(annotation.db$module.dict[[trait]])
    trait.attributes.current <- trait.attributes[[trait]]

    cluster.attributes <- list()

    clusters <- levels(as.factor(trait.attributes.current$clusters$membership))

    bins <- trait.attributes.current$cluster$names

    for(cluster in clusters){
      bins.cluster <- bins[which(trait.attributes.current$clusters$membership == cluster)]

      # Filter for completion
      bins_remove <- .Filter_Completion(features, trait, bins.cluster)

      if ( length(bins_remove) == length(bins.cluster) ) {
        next()
      } else if( length(bins_remove) > 0 ) {
        bins.cluster <- bins.cluster[!bins.cluster %in% bins_remove]
      }



      # number of bins <= 2 would end up with only 1 avg zscore which is not enough for a t.test
      if ( length(bins.cluster) < 2 ) {
        next()
      } else if ( length(bins.cluster) == 2){
        p.val <- .Calc_P_pairwise(bins.cluster, trait)

      }
      else {
        bin.zscores  <- trait.attributes.current$avg.zscore.module[bins.cluster,bins.cluster]

        p.val <- 1
        try(
          {p.val <- t.test(bin.zscores,bkgd.traits[[as.character(n.genes)]], alternative = 'less')$p.value},
          silent = T
        )

        p.val <- p.adjust(p.val, method = 'BH', length(clusters))

      }


      tryCatch(
        {
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
    print( trait.names[i] )
    trait.attribute <- .Calc_P( trait.names[i], p.threshold, annotation.db )
   # Backup check if there are not sig attributes
   # if (length(trait.attribute) == 0){
   #   pos.sig <- Identify_Significance_Trait(trait.names[i],
   #                                          RNAseq.data,
   #                                          pairwise.distances,
   #                                          bkgd.traits)
   #   if (length(pos.sig) > 0) {
   #     trait.attribute <- list('1' = pos.sig)
   #   }
   # }
    trait.attributes.pruned[[i]] <-  trait.attribute
  }

  names(trait.attributes.pruned) <- trait.names
  return(trait.attributes.pruned)
}
