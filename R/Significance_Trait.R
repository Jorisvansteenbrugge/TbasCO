#' Identify_Significance_Trait
#' @name Identify_Significance_Trait
#' @description Gives a verdict wether a trait is significant as a whole.
#' @param trait character string of the trait.
#' @param RNAseq.data Collection of multple components, include RNA seq data, annotations, etc. See \code{\link{Pre_process_input}} for the full list.
#' @param distance.metrics Named list containing distance functions as values.
#' @param bkgd.traits Random Background distributions of Traits.
#' @param bkgd.individual.Zscores
#' @export
#' @author JJM van Steenbrugge
Identify_Significance_Trait <- function (trait, RNAseq.data, distance.metrics,
                                        bkgd.individual.Zscores, bkgd.traits) {


    distances.list <- list()

    for (KO in RNAseq.data$features$annotation.db$module.dict[[trait]] ) {
      data.rows <- RNAseq.data$table[which(RNAseq.data$table$Annotation == KO),]

      genomes.present <- unique(data.rows$Bin)

      if (length(genomes.present) <= 1) {
        next()
      }

      dist.matrix <- matrix(data = NA,
                            ncol = length(genomes.present),
                            nrow = length(genomes.present))
      for(x in 1: (length(genomes.present) - 1)) {

        gen.x <- data.rows[which(data.rows$Bin == genomes.present[x]),][1,]

        for (y in (x+1):length(genomes.present)) {
          gen.y <- data.rows[which(data.rows$Bin == genomes.present[y]),][1,]

          distances.metrics <- list()
          for (metric in names(distance.metrics)) {
            dist   <- distance.metrics[[metric]](gen.x, gen.y, RNAseq.data$features)
            dist.Z <- (dist - bkgd.individual.Zscores$mu$`Random Annotated Genes`[[metric]] ) /
               bkgd.individual.Zscores$sd$`Random Annotated Genes`[[metric]]


            distances.metrics[[metric]] <- dist.Z

          }




          # composite
          comp.dist <- (-distances.metrics$PC) + distances.metrics$NRED
          # jaccard
          jacc.dist <- .Calc_Jaccard(RNAseq.data, random.genomes = c(genomes.present[x],
                                                                     genomes.present[y]),
                        used.terms = RNAseq.data$features$annotation.db$module.dict[[trait]])

          norm.dist <- (comp.dist * (1 - jacc.dist))

          dist.matrix[x,y] <- norm.dist
        }
      }

      distances.list[[KO]] <- dist.matrix

    } # KOs

    distances <- unlist(distances.list)
    distances <- as.numeric(distances[-which(is.na(distances))])
    distances.all[[metric]] <- distances



  p.val <- NA

  n.ko <- as.character(length(
    RNAseq.data$features$annotation.db$module.dict[[module]]))

  try(p.val <- t.test(distances, bkgd.traits[[n.ko]])$p.value,
      silent = T)

  try(
    if(p.val <= 0.05) {
      return(TRUE)
    } else {
      return(FALSE)
    }
    ,silent =T)



}
