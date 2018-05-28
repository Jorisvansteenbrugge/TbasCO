#' Identify whether an annotation term has a certain expression pattern in one
#' or more genomes
#' @name Go Fish
#' @description Identify whether a gene has a certain expression pattern in one or more genomes.
#' @param annotation character string of the Annotation term
#' @param RNAseq.data Collection of multple components, include RNA seq data, annotations, etc. See \code{\link{Pre_process_input}} for the full list.
#' @param pattern character string of the kind of pattern to look for. The string
#' should be one of the following: {'increase', 'decrease', 'peak', 'vale}
#' @export
#' @author JJM van Steenbrugge
Go_Fish <- function(annotation, RNAseq.data, pattern){

  .process_row <- function(row){
    sample.n <- length(RNAseq.data$features$sample.columns)
    row <- row[RNAseq.data$features$sample.columns]
    if      ( pattern == 'increase' ) {
      res <- sapply(2:sample.n,  function(x){
        if(row[x] > row[(x-1)]){
          return(T)
        }else {
          return(F)
        }
      })
    }
    else if ( pattern == 'decrease' ) {
      res <- sapply(2:sample.n,  function(x){
        if(row[x] < row[(x-1)]){
          return(T)
        }else {
          return(F)
        }
      })
    }
    else if ( pattern == 'peak'     ) {
      res <- sapply(2: (sample.n - 1),  function(x){
        if ( row[(x-1)] < row[x] &&  row[(x+1)] < row[x]) {
          return(T)
        }
        else{
          return(F)
        }
      })
    }
    else if ( pattern == 'vale'     ) {
      res <- sapply(2: (sample.n - 1),  function(x){
        if ( row[(x-1)] > row[x] &&  row[(x+1)] > row[x]) {
          return(T)
        }
        else{
          return(F)
        }
      })
    }

    if      ( T %in% res        ) {
      return(T)
    }
    else                          {
      return(F)
    }
  }



  annotation.rna <- RNAseq.data$table[which(
    RNAseq.data$table$Annotation == annotation),]


  tfs <- apply(annotation.rna, 1, .process_row)
  return(tfs)

}


