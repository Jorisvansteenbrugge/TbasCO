#' @export
#' @author JJM van Steenbrugge
Traitattributes_To_Sbsmatrix <- function(trait.attributes, bins) {
  sbs.matrix   <- matrix(nrow = length(bins),
                         ncol = 1)
  sbs.matrix   <- sbs.matrix[, -1]
  traits.names <- names(trait.attributes)
  cols <- c()

  # Traits
  for(trait.name in traits.names) {
    trait <- trait.attributes[[trait.name]]
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


#' @author Thomas Schlesinger
calcmfrow<- function(x){
  temp <- sqrt(x)

  if((temp %% floor(temp))==0){
    return(c(temp,temp))
  } else if((temp %% floor(temp))<0.5){
    return(c(floor(temp),ceiling(temp)))
  } else {
    return(c(ceiling(temp),ceiling(temp)))
  }
}

#' Plot Trait Attribute
#' @name Plot Trait Attribute
#' @description Plot the expression profile for all genomes expressing that attribute.
#' For each annotation a separate plot is created where the thick black line represents
#' the mean expression of all genomes with that annotation. The dashed blue and red lines
#' represent the mean + 0.5*sd and the mean - 0.5 * sd, respectively
#' @param trait.attribute The name of the trait attribute (as character string)
#' @param trait.attributes The full list of all trait attributes
#' @param RNAseq data
#' @example Plot_Trait_Attribute('M00027.2', trait.attributes.pruned, RNAseq.data)
#' @export
#' @author JJM van Steenbrugge
Plot_Trait_Attribute <- function(trait.attribute, trait.attributes,
                                 RNAseq.data){

  trait.attribute.s <- unlist(strsplit(x = trait.attribute,split = '[.]'))

  trait.annotations <- RNAseq.data$features$annotation.db$module.dict[[trait.attribute.s[1]]]
  attribute.genomes <- trait.attributes[[trait.attribute.s[1]]][[trait.attribute.s[2]]]$genomes

  par(mfrow=calcmfrow(length(trait.annotations)))

  for(annotation in trait.annotations){
    annotation.expression <- RNAseq.data$table[which(RNAseq.data$table$Annotation == annotation &
                                                       RNAseq.data$table$Bin %in% attribute.genomes),]
    if(nrow(annotation.expression) == 1){
      expression <- annotation.expression[RNAseq.data$features$sample.columns]
      plot(as.character(expression), type = 'l', ylab = 'Expression value', xlab = "Points in time",
          col="black", lwd="3", ylim = c(0, max(expression)), main = annotation)
      next()
    }
    mean.cols <- apply(annotation.expression[RNAseq.data$features$sample.columns],
                       2, mean)
    sd.cols   <- apply(annotation.expression[RNAseq.data$features$sample.columns],
                       2, sd) / 2

    mean.psd <- mean.cols + sd.cols
    mean.msd <- mean.cols - sd.cols

    max.val <- max(mean.psd)
    if(! is.nan(max.val)){
     plot(mean.cols,type = 'l', ylab = 'Expression value', xlab = "Points in time",
          col="black", lwd="3", ylim = c(0, max.val), main = annotation)
      points(mean.psd, type='l', col="blue", lty=2)
      points(mean.msd, type='l', col="red", lty=2)
    }else{
      plot(c(0,20,60,100,100,100), type='n', main = annotation, ylab="Expression value",
           xlab = "Points in time")
    }
  }
}
#' Plot_Background_Individual_Genes
#' @name Plot Background Individual Genes
#' @description Plotting of each background distribution.
#' @param bkgd.individual.Zscores
#' @export
#' @author BO Oyserman
Plot_Background_Individual_Genes <- function(bkgd.individual.Zscores){
  require(hexbin)
  require(RColorBrewer)
  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))

  random.genes.hexb <- hexbin(bkgd.individual.Zscores$zscores$`Random Genes`$PC,
                              bkgd.individual.Zscores$zscores$`Random Genes`$NRED)

  random.annotated.genes.hexb <- hexbin(bkgd.individual.Zscores$zscores$`Random Annotated Genes`$PC,
                                        bkgd.individual.Zscores$zscores$`Random Annotated Genes`$NRED)

  random.identical.annotated.genes.hexb <- hexbin(bkgd.individual.Zscores$zscores$`Genes with the same annotation`$PC,
                                                  bkgd.individual.Zscores$zscores$`Genes with the same annotation`$NRED)

  plot(random.genes.hexb,
       colramp=rf,mincnt=1, maxcnt=max(random.identical.annotated.genes.hexb@count),
       xlab="PC",ylab="NRED", main="Random Genes")

  plot(random.annotated.genes.hexb,
       colramp=rf,mincnt=1, maxcnt=max(random.identical.annotated.genes.hexb@count),
       xlab="PC",ylab="NRED", main="Random Annotated Genes")

  plot(random.identical.annotated.genes.hexb,
       colramp=rf,mincnt=1, maxcnt=max(random.identical.annotated.genes.hexb@count),
       xlab="PC",ylab="NRED", main="Genes with the same annotation")

}



