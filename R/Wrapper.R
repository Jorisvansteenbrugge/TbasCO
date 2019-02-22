
#' GetTraitAttributes
#' @description Wrapper of trait identification functions, including statistical validation
#' @param N Number of iterations for the random background distributions
#' @param metrics
#' @param p number of CPU threads to use simultaneously
#' @export
#' @examples
#' GetTraitAttributes()
#' GetTraitAttributes(RNAseq.data, N = 1000, p = 4)
#' @return trait.attributes.pruned
#' @author JJM van Steenbrugge
GetTraitAttributes <- function(RNAseq.data = RNAseq.data, N = 1000, metrics = distance.metrics,
                               p = 4, show_progres = TRUE){

  if (show_progres) print("Calcluating random background for genes")
  bkgd.individual         <- Individual_Annotation_Background(RNAseq.data,
                                                              N       = N,
                                                              metrics = metrics,
                                                              threads = p)

  bkgd.individual.Zscores <- Calc_Z_scores(bkgd.individual, distance.metrics)

  if (show_progres) print("Calculating random background for traits")
  bkgd.traits             <- Random_Trait_Background(RNAseq.data,
                                                     bkgd.individual.Zscores,
                                                     N = N,
                                                     metrics = metrics,
                                                     threads = p)

  if (show_progres) print("Calculating pairwise distances between annotations")
  pairwise.distances      <- Calc_Pairwise_Annotation_Distance(RNAseq.data,
                                                               RNAseq.data$features$annotation.db,
                                                               metrics,
                                                               bkgd.individual.Zscores,
                                                               show.progress = F,
                                                               threads = p)
  if (show_progres) print("Identifying Trait Attributes")
  trait.attributes        <- Identify_Trait_Attributes(RNAseq.data = RNAseq.data,
                                                       pairwise.distances = pairwise.distances,
                                                       threads = p)

  if (show_progres) print("Pruning Trait Attributes")
  trait.attributes.pruned <- Prune_Trait_Attributes(trait.attributes, bkgd.traits,
                                                    RNAseq.data,
                                                    p.threshold = 0.05,
                                                    pairwise.distances = metrics,bkgd.individual.Zscores = bkgd.individual.Zscores)


  return(list("Traitattributes" = trait.attributes,
              "Traitattributes sig" = trait.attributes.pruned,
              'Random background Individual' = bkgd.individual,
              'zscores individual' = bkgd.individual.Zscores,
              "Random backgrounds Traits" = bkgd.traits,
              "Pairwise distances" = pairwise.distances)
         )
}
