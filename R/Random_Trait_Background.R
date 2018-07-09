#' Generate a background distibution of composite Z scores
#' @name Random Background distributions of Traits
#' @description Create a random background distribution of distances between genes in a random trait. When multiple distance metrics
#' are provided, a composite Z score is calculated. The background distribution is determined by first random sampling two genomes.
#' Between the two genomes the overlap in annotations is determined. Based on the overlap a trait
#' is artificially created containing N annotations. For each annotation in the trait, the distances are
#' calculated between Genome A and Genome B, as described in the previous section
#' @param RNAseq.data Collection of multple components, include RNA seq data, annotations, etc. See \code{\link{Pre_process_input}} for the full list.
#' @param N Number of genes to include in the random module
#' @param Z Number of iterations
#' @param metrics Named list containing possibly multiple functions for distance
#' @param threads Number of cpu cores to be used
#' @export
#' @author JJM van Steenbrugge
#' BO Oyserman
Random_Trait_Background <- function(RNAseq.data,
                                    bkgd.individual.Zscores,
                                    N,
                                    Z,
                                    metrics){


  .Procedure <- function(RNAseq.data,
                         bkgd.individual.Zscores,
                         N,
                         Z,
                         metrics){

    #Runs in parallel
    result <- foreach::foreach(i = 1:N, .export = c("Random.Annotated.Genes.bkgd",
                                                     ".Convert_zscores",
                                                     ".Calc_Jaccard")) %do% {
      #pick two genomes
      RNAseq.data$annotation.only <- RNAseq.data$table[which(RNAseq.data$table$Annotation != ""),]
      random.genomes              <- sample(RNAseq.data$features$bins, 2)
      random.genomes.combis       <- rep(list(random.genomes), Z) #Workaround to re-use function below

      distances   <- Random.Annotated.Genes.bkgd(RNAseq.data, metrics, Z, random.genomes.combis)
      print(distances)
      distances.Z <- .Convert_zscores(distances$scores, metrics, bkgd.individual.Zscores)


      # This is quick and dirty composition
      composite.distance          <- mean( (-distances.Z$PC) + distances.Z$NRED,
                                           na.rm = T)[1]

      # Jaccard distance
      jaccard.distance            <- .Calc_Jaccard(RNAseq.data, random.genomes,
                                                   distances$`used terms`)

      #normalize composite scores with jaccard distance
      return( composite.distance * (1 - jaccard.distance) )
    }
  }



  require(doSNOW)

  bkgd.modules <- list()

  for(Z.size in Z){
    bkgd.modules[[as.character(Z.size)]] <- unlist(.Procedure(RNAseq.data,
                                                              bkgd.individual.Zscores,
                                                              N,
                                                              Z.size,
                                                              metrics))
  }




  return(bkgd.modules)
}

#' @export
.Convert_zscores <- function(distances, metrics, bkgd.individual.Zscores){
  Zscores <- bkgd.individual.Zscores

  distances.Z <- list()

  for(metric in names(metrics)){
    distances.Z[[metric]] <- ( distances[[metric]] - Zscores$mu$`Random Annotated Genes`[[metric]] ) /
                                  Zscores$sd$`Random Annotated Genes`[[metric]]
  }

  return(distances.Z)
}

#' @export
.Calc_Jaccard <- function(RNAseq.data, random.genomes, used.terms){
  presence.absence <- RNAseq.data$features$annotation_presence_absence
  PA.genome.A <- presence.absence[used.terms, as.character(random.genomes[1])]
  PA.genome.B <- presence.absence[used.terms, as.character(random.genomes[2])]

  Jaccard_Distance <- 1 - (sum(which(PA.genome.A == 1) %in% which(PA.genome.A == 1))  /
                             (sum(which(PA.genome.A == 0) %in% which(PA.genome.A == 1))  +
                                sum(which(PA.genome.A == 1)%in%which(PA.genome.A == 0))   +
                                sum(which(PA.genome.A == 1)%in%which(PA.genome.A == 1))))
  return(Jaccard_Distance)
}
