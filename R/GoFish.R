#' Identify whether an annotation term has a certain expression pattern in one
#' or more genomes
#' @name Go Fish
#' @description Identify whether a gene has a certain expression pattern in one or more genomes.
#' @param trait character string of the Annotation term
#' @param RNAseq.data Collection of multple components, include RNA seq data, annotations, etc. See \code{\link{Pre_process_input}} for the full list.
#' @param pattern character string of the kind of pattern to look for. The string
#' should be one of the following: {'increase', 'decrease', 'peak', 'vale}
#' @export
#' @example Go_Fish('K00927', RNAseq.data, type = 'peak')
#' @author JJM van Steenbrugge
Go_Fish <- function(trait, RNAseq.data, pattern){

  if(missing(pattern)){
    points <- .Draw_pattern()
    lines(points$x, points$y)
    text(3,0.5,'Using this pattern to match')
  }



  annotations <- RNAseq.data$features$annotation.db$module.dict[[trait]]

  annotation.rna <- RNAseq.data$table[which(
  RNAseq.data$table$Annotation %in% annotations),]

  tfs <- apply(annotation.rna, 1, .process_row)
  return(tfs)

}


.Draw_pattern <- function(nsamples = 6){
  plot(c(1,nsamples), c(0,1), type='n',
       xlab = 'Time points',
       ylab = 'Normalized Rank')
  cat("Click on six points in the graph, then press finish in the plot window")

  return(locator(nsamples))

}

.Draw_pattern()
