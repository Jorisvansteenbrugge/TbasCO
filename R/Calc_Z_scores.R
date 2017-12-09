
#' Calculate Z scores of individual background distributions
#' @description Function to convert different random background distribution values of possible different distance metrics into Z_scores
#' @param bkgd.individual Named list containing (three) different background distributions.
#' @param metrics Named list containing possibly multiple functions for distance.
#' @export
#' @example Calc_Z_scores(bkgd.individual, metrics)
#' @authors JJM van Steenbrugge
#' BO Oyserman
Calc_Z_scores <- function(bkgd.individual, metrics){
  mu.values <- list()
  sd.values <- list()

  # This can probably be a lot shorter
  for(distribution in names(bkgd.individual)){

    mu.values[[distribution]] <- list()
    sd.values[[distribution]] <- list()

    current.bkgd <- bkgd.individual[[distribution]]

    for(metric in names(metrics)){
      mu    <- mean(current.bkgd[[metric]])
      stdev <- sd(current.bkgd[[metric]])

      mu.values[[distribution]][[metric]] <- mu
      sd.values[[distribution]][[metric]] <- stdev
    }

  }


  #applying them
  distributions <- c("Random Annotated Genes",
                     "Genes with the same annotation")
  Z_scores.output <- list()

  for(distribution in distributions){
    Z_scores.output[[distribution]] <- list()

    current.bkgd <- bkgd.individual[[distribution]]
    for(metric in names(metrics)){

      current.idv.bkgd <- bkgd.individual[[distribution]][[metric]]
      current.mu <- mu.values[[distribution]][[metric]]
      current.sd <- sd.values[[distribution]][[metric]]

      Z_scores.output[[distribution]][[metric]] <- (current.idv.bkgd - current.mu) /
        current.sd
    }
  }
  return(Z_scores.output)
}
