#' Prune Trait Attributes
#' @name Prune Trait Attributes
#' @description For each trait attribute, statistical significance is calculated
#' @param filepath RNAseq data file
#' @param distance.metrics Variable
#' @param annotation.db
#' @param bkgd.individual.Zscores
#' @export
#' @author JJM van Steenbrugge
Prune_Trait_Attributes <- function(trait.attributes, bkgd.modules, features,
                                   p.threshold = 0.05, completion.threshold = 0.75,
                                   threads = 4){

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
    return(completion >= 0.75)

  }

  .Calc_P <- function(trait, p.threshold, annotation.db){
    n.genes <- length(annotation.db$module.dict[[trait]])
    trait.attributes.current <- trait.attributes[[trait]]

    cluster.attributes <- list()

    clusters <- levels(as.factor(trait.attributes.current$clusters$membership))

    bins <- trait.attributes.current$cluster$names

    for(cluster in clusters){
      bins.cluster <- bins[which(trait.attributes.current$clusters$membership == cluster)]


      #Filter for completion
      if(! .Filter_Completion(annotation.db$module.dict[[trait]], bins, features, 0.75 )){
        next()
      }


      #

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
    .Calc_P(trait.names[i], p.threshold, annotation.db)
  }

  # for(i in 1:length(trait.names)){
  #   .Calc_P(trait.names[i], p.threshold, annotation.db)
  # }

  names(trait.attributes.pruned) <- trait.names
  return(trait.attributes.pruned)
}
