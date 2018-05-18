#' Prune Trait Attributes
#' @name Prune Trait Attributes
#' @description For each trait attribute, statistical significance is calculated
#' @param filepath RNAseq data file
#' @param distance.metrics Variable
#' @param annotation.db
#' @param bkgd.individual.Zscores
#' @export
#' @author JJM van Steenbrugge
Prune_Trait_Attributes <- function(trait.attributes, bkgd.modules, annotation.db,
                                   p.threshold = 0.05, threads = 4){

  .Calc_P <- function(trait, p.threshold){
    n.genes <- length(annotation.db$module.dict[[trait]])
    trait.attributes.current <- trait.attributes[[trait]]

    cluster.attributes <- list()


    clusters <- levels(as.factor(trait.attributes.current$clusters$membership))

    bins <- trait.attributes.current$cluster$names

    for(cluster in clusters){
      bins.cluster <- bins[which(trait.attributes.current$clusters$membership == cluster)]
      bin.zscores  <- trait.attributes.current$avg.zscore.module[bins.cluster,bins.cluster]

      p.val = tryCatch({
      p.val <- t.test(bin.zscores,bkgd.modules[[as.character(n.genes)]], alternative = 'less')$p.value
      }, warning = function(w){
        return(w)
      }, error = function(e){
        # Expecting this to happen if there is no bg_distance_module calculated
        #  for the module size, or if there are not enough bins in a cluster
      })

      if(! is.null(p.val) && p.val < p.threshold){
        cluster.attributes[[cluster]] <- list("genomes"= bins.cluster,
                                              "p-val"  = p.val)
      }

    }
    return(cluster.attributes)
  }

  trait.names <- names(annotation.db$module.dict)


  trait.attributes.pruned <- foreach::foreach(i = 1: length(trait.names)) %do%{
    .Calc_P(trait.names[i], p.threshold)
  }


  names(trait.attributes.pruned) <- trait.names
  return(trait.attributes.pruned)
}
