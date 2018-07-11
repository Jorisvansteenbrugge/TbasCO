#' Identify whether an annotation term has a certain expression pattern in one
#' or more genomes
#' @name Go Fish
#' @description Identify whether a gene has a certain expression pattern in one or more genomes.
#' @param RNAseq.data Collection of multple components, include RNA seq data, annotations, etc. See \code{\link{Pre_process_input}} for the full list.
#' @export
#' @example Go_Fish('K00927', RNAseq.data, type = 'peak')
#' @author JJM van Steenbrugge
Go_Fish <- function(RNAseq.data){

  points <- .Draw_pattern()
  lines(points$x, points$y)
  text(3,0.5,'Using this pattern to match')

  for (i in 1:length(names(RNAseq.data$features$annotation.db$module.dict))) {
    trait <- RNAseq.data$features$annotation.db$module.dict[[i]]
    genes <- RNAseq.data$table[which(RNAseq.data$table$Annotation %in% trait),]


    verdicts <- apply(genes, 1, function(gene){
      gene <- gene[RNAseq.data$features$rank.columns]
      return(.Test_pattern(points$y, gene))
    })

    if (TRUE %in% verdicts) {
      print(names(RNAseq.data$features$annotation.db$module.dict)[i])
    }
  }

}

.Test_pattern <- function(pattern, gene_pattern, margin = 0.15){
  verdict <- T
  for( i in 1: length(pattern) ) {
    lower <- pattern[i] - margin
    upper <- pattern[i] + margin

    if(gene_pattern[i] < lower || gene_pattern[i] > upper) {
      verdict <- F
    }

  }
  return(verdict)
}

.Draw_pattern <- function(nsamples = 6){
  plot(c(1,nsamples), c(0,1), type='n',
       xlab = 'Time points',
       ylab = 'Normalized Rank')
  cat("Pick a rank value (Y axis) for each timme point by clicking in the graph")

  return(locator(nsamples))

}

