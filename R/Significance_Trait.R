#' Identify_Significance_Trait
#' @name Identify_Significance_Trait
#' @description Gives a verdict wether a trait is significant as a whole.
#' @param trait character string of the trait.
#' @param RNAseq.data Collection of multple components, include RNA seq data,
#' annotations, etc. See \code{\link{Pre_process_input}} for the full list.
#' @param distance.metrics Named list containing distance functions as values.
#' @param bkgd.traits Random Background distributions of Traits.
#' @param bkgd.individual.Zscores
#' @export
#' @author JJM van Steenbrugge
Identify_Significance_Trait <- function (trait, RNAseq.data,
                                        pairwise.distances, bkgd.traits) {


    distances <- c()

    for (KO in RNAseq.data$features$annotation.db$module.dict[[trait]] ) {
      genomes <- RNAseq.data$table[which(RNAseq.data$table$Annotation == KO),
                                   'Bin']
      distances <- c(distances, pairwise.distances[[KO]][genomes, genomes])

    }

    distances <- as.numeric(distances[-which(is.na(distances))])




  p.val <- NA

  n.ko <- as.character(length(
    RNAseq.data$features$annotation.db$module.dict[[trait]]))

  try(p.val <- t.test(distances, bkgd.traits[[n.ko]])$p.value,
      silent = T)
  try(
    if(! is.na(p.val)  && p.val <= 0.05) {
      return(list("genomes" = genomes,
                  "p-val"   = p.val))
    } else {
      return(list())
    }
    ,silent =T)



}
