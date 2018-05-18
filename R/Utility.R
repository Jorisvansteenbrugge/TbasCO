#' @export
#' @author JJM van Steenbrugge
Traitattributes_To_Sbsmatrix <- function(trait.attributes.pruned, bins) {
  sbs.matrix   <- matrix(nrow = length(bins),
                         ncol = 1)
  sbs.matrix   <- sbs.matrix[, -1]
  traits.names <- names(trait.attributes.pruned)
  cols <- c()

  # Traits
  for(trait.name in traits.names) {
    trait <- trait.attributes.pruned[[trait.name]]
    if(length(trait) == 0){
      next()
    }
    # Attributes
    for( i in 1: length(trait)) {
      attribute.bins <- trait[[i]]$genomes
      occurences     <- bins %in% attribute.bins
      sbs.matrix <- cbind(sbs.matrix, occurences)
      cols <- c(cols, paste(trait.name,'.',i, sep = ''))
      #include names
    }
    print(trait.name)
  }
  colnames(sbs.matrix) <- cols
  sbs.matrix <- as.data.frame(sbs.matrix)
  cols <- sapply(sbs.matrix, is.logical)
  sbs.matrix[,cols] <- lapply(sbs.matrix[,cols], as.numeric)
  return(as.matrix(sbs.matrix))
}


