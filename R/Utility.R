#' Plot Trait Attribute
#' @name Plot Trait Attribute
#' @description Convert a list with trait attributes to a matrix where each
#' collumn is a trait attribute and each row a genome. Presence or absence is
#' indicated with either a zero (0) or one (1).
#' @param trait.attributes The full list of all trait attributes
#' @param genomes A vector containing the names of all genomes present in the
#' sample.
#' @export
#' @author JJM van Steenbrugge
Traitattributes_To_Sbsmatrix <- function(trait.attributes, genomes) {
  sbs.matrix   <- matrix(nrow = length(genomes),
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
      attribute.genomes <- trait[[i]]$genomes
      occurences     <- genomes %in% attribute.genomes
      sbs.matrix <- cbind(sbs.matrix, occurences)
      cols <- c(cols, paste(trait.name,'.',i, sep = ''))
      #include names
    }
  }
  colnames(sbs.matrix) <- cols
  sbs.matrix <- as.data.frame(sbs.matrix)
  cols <- sapply(sbs.matrix, is.logical)
  sbs.matrix[,cols] <- lapply(sbs.matrix[,cols], as.numeric)

  rownames(sbs.matrix) <- genomes
  return(as.matrix(sbs.matrix))
}

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
#' @description Plot the expression profile for all genomes expressing that
#' attribute. For each annotation a separate plot is created where the thick
#' black line represents the mean expression of all genomes with that annotation.
#' The dashed blue and red lines represent the mean + 0.5*sd and the mean -
#' 0.5 * sd, respectively.
#' @param trait.attribute The name of the trait attribute (as character string)
#' @param trait.attributes The full list of all trait attributes
#' @param RNAseq.data Collection of multple components, include RNA seq data,
#' annotations, etc. See \code{\link{Pre_process_input}} for the full list.
#' @example Plot_Trait_Attribute('M00027.2', trait.attributes.pruned, RNAseq.data)
#' @export
#' @author JJM van Steenbrugge
Plot_Trait_Attribute_Expression <- function(trait.attribute, trait.attributes,
                                 RNAseq.data){

  trait.attribute.s <- unlist(strsplit(x = trait.attribute,split = '[.]'))

  trait.annotations <- RNAseq.data$features$annotation.db$module.dict[[trait.attribute.s[1]]]
  attribute.genomes <- trait.attributes[[trait.attribute.s[1]]][[trait.attribute.s[2]]]$genomes

  par(mfrow=calcmfrow(length(trait.annotations)))

  for(annotation in trait.annotations){
    annotation.expression <- RNAseq.data$table[which(RNAseq.data$table$Annotation == annotation &
                                                       RNAseq.data$table$Bin %in% attribute.genomes),]
    if(nrow(annotation.expression) == 1){
      expression <- annotation.expression[RNAseq.data$features$rank.columns]
      plot(as.character(expression), type = 'l', ylab = 'Expression value', xlab = "Points in time",
          col="black", lwd="3", ylim = c(0, max(expression)), main = annotation)
      next()
    }
    mean.cols <- apply(annotation.expression[RNAseq.data$features$rank.columns],
                       2, mean)
    sd.cols   <- apply(annotation.expression[RNAseq.data$features$rank.columns],
                       2, sd) / 2

    mean.psd <- mean.cols + sd.cols
    mean.msd <- mean.cols - sd.cols

    max.val <- max(mean.psd)
    if(! is.nan(max.val)){
     plot(mean.cols,type = 'l', ylab = 'Expression value', xlab = "Points in time",
          col="black", lwd="3", ylim = c(0, 1), main = annotation)
      points(mean.psd, type='l', col="blue", lty=2)
      points(mean.msd, type='l', col="red", lty=2)
    }else{
      plot(c(0,20,60,100,100,100), type='n', main = annotation, ylab="Expression value",
           xlab = "Points in time")
    }
  }
}

#' Network Trait Genomes
#' @name Network Trait Genomes
#' @description Draw a network of trait(-attributes) with its associated genomes
#' in Cytoscape (\url{http://www.cytoscape.org/}) automatically via Cytoscape's
#' REST API..
#' @param trait.names a vector containing the names of traits, as provided in the
#' RNAseq.data$features$annotation.db$module.dict object
#' @export
#' @example Network_Trait_Genomes(c("M00002", "M00007"), trait.attributes.pruned)
#' @author JJM van Steenbrugge
Network_Trait_Genomes    <- function(trait.names, trait.attributes.pruned,
                                     genome.phylogeny){
  if(!"RCy3" %in% installed.packages()){
    source("https://bioconductor.org/biocLite.R")
    biocLite("RCy3")
  }
  library(RCy3)
  library(randomcoloR)


  nodes <- matrix(nrow=0, ncol=3)
  colnames(nodes) <- c("id","group","phylogeny")
  edges <- matrix(nrow=0, ncol=2)
  colnames(edges) <- c("source","target")

  added_nodes <- c()


  for(trait in trait.names) {
    attributes <- trait.attributes.pruned[[trait]]

    for(attribute in names(attributes)) {

      name <- paste(trait, attribute,
                    sep = '.')
      nodes <- rbind(nodes, c(name, 'trait','-'))

      for(genome in attributes[[attribute]]$genomes) {

        if(! genome %in% added_nodes ) {
          if(! missing(genome.phylogeny)) {
            nodes <- rbind(nodes, c(genome, 'genome', genome.phylogeny[[genome]]))
          } else {
            nodes <- rbind(nodes, c(genome, 'genome', '.'))
          }
          added_nodes <- c(added_nodes, genome)
        }

        edges <- rbind(edges, c(genome, name))
      }
    }
  }

  createNetworkFromDataFrames(data.frame(nodes,stringsAsFactors = F),
                              data.frame(edges,stringsAsFactors = F),
                              title= paste(trait.names, collapse = ','),
                              collection="Network Traits and Genomes")

  setVisualStyle('Nested Network Style')

  values <- c ('trait',  'genome')
  shapes <- c ('rectangle', 'circle')
  setNodeShapeMapping ('group', values, shapes)

  # Genome node colors ----
  setNodeColorBypass(node.names = nodes[which(nodes[,2] == 'genome'), 1],
                     new.colors = '#b30000'

                     )
  # Custom family colors
  if(! missing(genome.phylogeny)) {
    families <- unique(sapply(genome.phylogeny, function(x) {return(x)} ))

    for (i in 1:length(families)) {
      setNodeColorBypass(node.names = nodes[which(nodes[,3] == families[i]),1],
                         new.colors = randomcoloR::randomColor(count = 1,
                                                               luminosity = 'dark'))
    }
  }

  # Label text colors
  setNodeLabelColorBypass(node.names = nodes[which(nodes[,2] == 'genome'), 1],
                          new.colors = '#ffffff')

  setNodeSizeBypass(node.names = nodes[which(nodes[,2] == 'genome'), 1],
                    new.sizes = 40)

  # Trait node colors ----
  setNodeColorBypass(node.names = nodes[which(nodes[,2] == 'trait'), 1],
                     new.colors = '#737373')


}

#' Plot_Background_Individual_Genes
#' @name Plot Background Individual Genes
#' @description Plotting of each background distribution.
#' @param bkgd.individual.Zscores A list containing (composite) Zscores based on
#' a random background distribution.
#' See \code{\link{Calc_Z_scores}} for more information.
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

  cnt.max <- max(c(random.identical.annotated.genes.hexb@count,
                   random.annotated.genes.hexb@count,
                   random.genes.hexb@count))
  plot(random.genes.hexb,
       colramp=rf,mincnt=1, maxcnt=cnt.max,
       xlab="PC",ylab="NRED", main="Random Genes")

  plot(random.annotated.genes.hexb,
       colramp=rf,mincnt=1, maxcnt=cnt.max,
       xlab="PC",ylab="NRED", main="Random Annotated Genes")

  plot(random.identical.annotated.genes.hexb,
       colramp=rf,mincnt=1, maxcnt=cnt.max,
       xlab="PC",ylab="NRED", main="Genes with the same annotation")
}

#' Plot_Background_Modules
#' @name Plot Background Modules
#' @description Plotting of all background distributions of different sizes
#' together. Each distribution is represented by a different colour. This plot
#' hints wether or not the background distribution is indeed random.
#' @param bkgd.traits A collection of random background distributions of traits.
#' See \code{\link{Random_Trait_Background}}.
#' @export
#' @author JJM van Steenbrugge
Plot_Background_Modules          <- function(bkgd.traits){
  if(length(bkgd.traits) < 1 ){
    return("There are no random background distributions calculated")
  }
  colours           <- rainbow(length(bkgd.traits))
  bkgd.names        <- names(bkgd.traits)
  legend            <- sapply(bkgd.names,
                              function(x) paste('N = ', x, sep=""))

  plot(density(bkgd.traits[[1]],
               na.rm=TRUE),
       ylim     = c(0, 1),
       xlim     = c(-4, 4),
       ylab     = "Density",
       xlab     = "Composite Zscore",
       cex.main = 0.75,
       col      = colours[1],
       main     = "Random background distributions of module sizes 2 - 20")

  for(i in 2: length(bkgd.traits)){
    points(density(bkgd.traits[[i]],
                   na.rm=TRUE),
           type = "l",
           col  = colours[i])
  }

  legend("topleft",
         legend = legend,
         col    = colours,
         pch    = 1,
         cex    = 0.7)
}

#' Calculate Association Rules
#' @name Calculate Association Rules
#' @description Calculate Association Rules using the apriori algorithm (as
#' implemented in the arules package). This function lets the user decide on the
#' number of rules to produce and estimates parameters to produce those results.
#' @param sbs.trait.attributes Matrix where each collumn is a trait attribute
#' and each row a genome. Presence or absence is indicated with either a zero
#' (0) or one (1).
#' @param lhs A vector containing one or multiple terms for the left-hand-side
#' (antecedent) of the association rules. For more information see
#' \url{https://en.wikipedia.org/wiki/Association_rule_learning}.
#' @param rhs A vector containing one or multiple terms for the right-hand-side
#' (consequent) of the association rules. For more information see
#' \url{https://en.wikipedia.org/wiki/Association_rule_learning}.
#' @param N The number of association rules to create (e.g. N = 100 will result
#' in 100 association rules).
#' @export
#' @author JJM van Steenbrugge
Association_Rules <- function(sbs.trait.attributes,
                              lhs, rhs, N){
  require(arules)

  .Reach_N <- function(N, lhs = c(), rhs= c()){
    n <- 0
    support = 1.0
    confidence = 1.0
    while (n < N) {
      apri <- apriori(sbs.trait.attributes,
                      parameter = list(support    = support,
                                       confidence = confidence,
                                       minlen     = 2),
                      appearance = list(lhs = lhs,
                                        rhs = rhs)
                      )

      rules <- as(apri, 'data.frame')
      n <- nrow(rules)
      print(n)
      support <- support - 0.1
      if(support <= 0){
        support <- 1.0
        confidence <- confidence - 0.1
      }
    }
    return(rules)
  }




  if (missing(lhs) && missing(rhs)) {
    rules <- .Reach_N(N)
  } else if (!missing(lhs) && !missing(rhs)) {
    rules1 <- .Reach_N(N/2, lhs = lhs)
    rules2 <- .Reach_N(N/2, rhs = rhs)

    rules1 <- rules1[order(rules1[,1]), ]
    rules2 <- rules2[order(rules2[,1]), ]

    rules  <- rbind(rules1[1:(N/2),], rules2[1:(N/2),])

  } else if (!missing(lhs)) {
    rules <- .Reach_N (N,lhs = lhs)
  } else if (!missing(rhs)) {
    rules <- .Reach_N (N,rhs = rhs)
  }

  return(rules)
}

#' Draw_Expression
#' @name Draw Expression
#' @description Creates a drawing of a trait (possibly with multiple different
#' trait-attributes) where all its composing genes are plotting individually
#' with their corresponding expression values.
#' @param trait character string with the name of a trait.
#' @param RNAseq.data Collection of multple components, include RNA seq data,
#' annotations, etc.
#' @seealso \code{\link{Pre_process_input}} for the full list of parameters,
#' \url{https://github.com/Jorisvansteenbrugge/TcT/wiki/The-RNAseq.data-object}
#' for more information.
#' @param trait.attributes.pruned The full list of all pruned trait attributes.
#' @export
#' @author JJM van Steenbrugge
Draw_Expression <- function(trait, RNAseq.data, trait.attributes.pruned) {

  y.max        <- 110
  x.max        <- 35


  attributes   <- trait.attributes.pruned[[trait]]
  n.attributes <- length(attributes)
  genes        <- RNAseq.data$features$annotation.db$module.dict[[trait]]
  n.genes      <- length(genes)

  plot(c(0,x.max), c(0,y.max), type='n', xlab = '',ylab='', axes = F)

  y.coords <- seq(100, 1    , by = -15) # With this setting max of 7 genes
  x.coords <- seq(0  , x.max, by =  11)

  # For each gene
  for(i in 1:n.genes) {
    gene <- genes[i]

    # Gene name
    graphics::text(x = ( x.coords[1] + 5 ),
                   y = ( y.coords[i] - 5 ),
                   labels = gene)

    # for each attribute
    for(y in 1:n.attributes) {
#      Draw Labels
       graphics::text(x = ( x.coords[y+1] + 5), y = y.max, labels = paste(trait,
                                                                        y, sep = '.'))

       rect(xleft  = x.coords[y+1]     , ybottom = y.coords[i] - 10,
             xright = x.coords[y+1] + 10, ytop    = y.coords[i])

        # Draw expression line
        genomes     <- attributes[[y]]$genomes
        genomes.rna <- RNAseq.data$table[which(RNAseq.data$table$Bin %in% genomes),]
        genomes.rna.gene <- genomes.rna[which(genomes.rna$Annotation == gene),
                                        RNAseq.data$features$rank.columns]
        genomes.rna.gene.mean <- apply(genomes.rna.gene, 2, mean)


        # Expression ranks
        y.bot <- y.coords[i] - 10
        rank.pos <- (genomes.rna.gene.mean * 10) + y.bot
        x.pos    <- seq(1,10, length.out = 6)+ x.coords[y+1]
        graphics::lines(x.pos,
                        rank.pos
        )

        # Generate box colours
        values <- seq(0.1,1,by = 0.1)
        cols <- colorRamp(c('white','red'))
        colours <- rgb(cols(values)/255)

        # Draw all the boxes!
        box.coords <- seq(0,10,length.out = (length(genomes.rna.gene.mean) +1))
        for(z in 2:length(box.coords)) {
          rank <- genomes.rna.gene.mean[z-1]
          col = 'red'
          # I am not proud of this
           if(rank <= 1 && rank > 0.9){
             col <- colours[1]
             }
          else if(rank <= 0.9 && rank > 0.8){
            col <- colours[2]
          } else if(rank <= 0.8 && rank > 0.7){
            col <- colours[3]
          }else if(rank <= 0.7 && rank > 0.6){
            col <- colours[4]
          }else if(rank <= 0.6 && rank > 0.5){
            col <- colours[5]
          } else if(rank <= 0.5 && rank > 0.4){
            col <- colours[6]
          }else if(rank <= 0.4 && rank > 0.3){
            col <- colours[7]
          }else if(rank <= 0.3 && rank > 0.2){
            col <- colours[8]
          }else if(rank <= 0.2 && rank > 0.1){
            col <- colours[9]
          }else if(rank <= 0.1 && rank >= 0.0){
            col <- colours[10]
          }
          rect(xleft  = x.coords[y+1] + box.coords[(z-1)] , ybottom = y.coords[i] - 10,
               xright = x.coords[y+1] + box.coords[z], ytop    = y.coords[i],
               col=col)

          # Expression ranks
          y.bot <- y.coords[i] - 10
          rank.pos <- (genomes.rna.gene.mean * 10) + y.bot
          x.pos    <- seq(1,10, length.out = 6)+ x.coords[y+1]
          graphics::lines(x.pos,
                          rank.pos, col = 'black'
          )

        }




    }
  }
}
