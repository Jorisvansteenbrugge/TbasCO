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
  sbs.matrix <- matrix(
    nrow = length(genomes),
    ncol = 1
  )
  sbs.matrix <- sbs.matrix[, -1]
  traits.names <- names(trait.attributes)
  cols <- c()

  # Traits
  for (trait.name in traits.names) {
    trait <- trait.attributes[[trait.name]]
    if (length(trait) == 0) {
      next()
    }
    # Attributes
    for (i in 1:length(trait)) {
      attribute.genomes <- trait[[i]]$genomes
      occurences <- genomes %in% attribute.genomes
      sbs.matrix <- cbind(sbs.matrix, occurences)
      cols <- c(cols, paste(trait.name, ".", i, sep = ""))
      # include names
    }
  }
  colnames(sbs.matrix) <- cols
  sbs.matrix <- as.data.frame(sbs.matrix)
  cols <- sapply(sbs.matrix, is.logical)
  sbs.matrix[, cols] <- lapply(sbs.matrix[, cols], as.numeric)

  rownames(sbs.matrix) <- genomes
  return(as.matrix(sbs.matrix))
}

calcmfrow <- function(x) {
  temp <- sqrt(x)

  if ((temp %% floor(temp)) == 0) {
    return(c(temp, temp))
  } else if ((temp %% floor(temp)) < 0.5) {
    return(c(floor(temp), ceiling(temp)))
  } else {
    return(c(ceiling(temp), ceiling(temp)))
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
Plot_Trait_Attribute_Expression <- function(trait.attribute,
                                            trait.attributes.pruned,
                                            RNAseq.data) {
  trait.attribute.s <- unlist(strsplit(x = trait.attribute, split = "[.]"))

  trait.annotations <- RNAseq.data$features$annotation.db$module.dict[[trait.attribute.s[1]]]
  attribute.genomes <- trait.attributes[[trait.attribute.s[1]]][[trait.attribute.s[2]]]$genomes

  par(mfrow = calcmfrow(length(trait.annotations)))

  for (annotation in trait.annotations) {
    annotation.expression <- RNAseq.data$table[which(RNAseq.data$table$Annotation == annotation &
      RNAseq.data$table$Bin %in% attribute.genomes), ]
    if (nrow(annotation.expression) == 1) {
      expression <- annotation.expression[RNAseq.data$features$rank.columns]
      plot(as.character(expression),
        type = "l", ylab = "Expression value", xlab = "Points in time",
        col = "black", lwd = "3", ylim = c(0, max(expression)), main = annotation
      )
      next()
    }
    mean.cols <- apply(
      annotation.expression[RNAseq.data$features$rank.columns],
      2, mean
    )
    sd.cols <- apply(
      annotation.expression[RNAseq.data$features$rank.columns],
      2, sd
    ) / 2

    mean.psd <- mean.cols + sd.cols
    mean.msd <- mean.cols - sd.cols

    max.val <- max(mean.psd)
    if (!is.nan(max.val)) {
      plot(mean.cols,
        type = "l", ylab = "Expression value", xlab = "Points in time",
        col = "black", lwd = "3", ylim = c(0, 1), main = annotation
      )
      points(mean.psd, type = "l", col = "blue", lty = 2)
      points(mean.msd, type = "l", col = "red", lty = 2)
    } else {
      plot(c(0, 20, 60, 100, 100, 100),
        type = "n", main = annotation, ylab = "Expression value",
        xlab = "Points in time"
      )
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
#' @examples Network_Trait_Genomes(c("M00002", "M00007"), trait.attributes.pruned)
#' \dontrun{
#' shared_ten_genomes <- c('M00335', 'M00120', 'M00144', 'M00149', 'M00260', 'M00082', 'M00360', 'M00022', 'M00023')
#' niche_traits <- c("M00150", 'M00064', 'M00123', 'M00186', 'M00244', 'M00223', 'M00323', 'M00328', 'M00453', 'M00515', 'M00523', 'M00579', 'M00772')
#' }
#' @author JJM van Steenbrugge
Network_Trait_Genomes <- function(trait.names, trait.attributes.pruned,
                                  genome.phylogeny) {
  if (!"RCy3" %in% installed.packages()) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("RCy3")
  }
  library(RCy3)
  library(randomcoloR)


  nodes <- matrix(nrow = 0, ncol = 3)
  colnames(nodes) <- c("id", "group", "phylogeny")
  edges <- matrix(nrow = 0, ncol = 2)
  colnames(edges) <- c("source", "target")

  added_nodes <- c()


  for (trait in trait.names) {
    attributes <- trait.attributes.pruned[[trait]]

    for (attribute in names(attributes)) {
      name <- paste(trait, attribute,
        sep = "."
      )
      nodes <- rbind(nodes, c(name, "trait", "-"))

      for (genome in attributes[[attribute]]$genomes) {
        if (!genome %in% added_nodes) {
           if (!missing(genome.phylogeny)) {
             nodes <- rbind(nodes, c(genome, "genome", genome.phylogeny[[genome]]))
           } else {
            nodes <- rbind(nodes, c(genome, "genome", "."))
          }
          added_nodes <- c(added_nodes, genome)
        }

        edges <- rbind(edges, c(genome, name))
      }
    }
  }

  createNetworkFromDataFrames(data.frame(nodes, stringsAsFactors = F),
    data.frame(edges, stringsAsFactors = F),
    title = paste(trait.names, collapse = ","),
    collection = "Network Traits and Genomes"
  )

  setVisualStyle("Nested Network Style")

  values <- c("trait", "genome")
  shapes <- c("rectangle", "circle")
  setNodeShapeMapping("group", values, shapes)

  # Genome node colors ----
  setNodeColorBypass(
    node.names = nodes[which(nodes[, 2] == "genome"), 1],
    new.colors = "#b30000"
  )
  # Custom family colors
  if (!missing(genome.phylogeny)) {
    families <- unique(sapply(genome.phylogeny, function(x) {
      return(x)
    }))

    for (i in 1:length(families)) {
      setNodeColorBypass(
        node.names = nodes[which(nodes[, 3] == families[i]), 1],
        new.colors = randomcoloR::randomColor(
          count = 1,
          luminosity = "dark"
        )
      )
    }
  }



  # Label text colors
  setNodeLabelColorBypass(
    node.names = nodes[which(nodes[, 2] == "genome"), 1],
    new.colors = "#ffffff"
  )

  setNodeSizeBypass(
    node.names = nodes[which(nodes[, 2] == "genome"), 1],
    new.sizes = 40
  )

  # Trait node colors ----
  setNodeColorBypass(
    node.names = nodes[which(nodes[, 2] == "trait"), 1],
    new.colors = "#737373"
  )
}

#' Plot_Background_Individual_Genes
#' @name Plot Background Individual Genes
#' @description Plotting of each background distribution.
#' @param bkgd.individual.Zscores A list containing (composite) Zscores based on
#' a random background distribution.
#' See \code{\link{Calc_Z_scores}} for more information.
#' @export
#' @author BO Oyserman
Plot_Background_Individual_Genes <- function(bkgd.individual.Zscores) {
  require(hexbin)
  require(RColorBrewer)
  rf <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

  all_scores_x <- c(
    bkgd.individual.Zscores$zscores$`Random Genes`$PC,
    bkgd.individual.Zscores$zscores$`Random Annotated Genes`$PC,
    bkgd.individual.Zscores$zscores$`Genes with the same annotation`$PC
  )
  
   all_scores_y <- c(
    bkgd.individual.Zscores$zscores$`Random Genes`$NRED,
    bkgd.individual.Zscores$zscores$`Random Annotated Genes`$NRED,
    bkgd.individual.Zscores$zscores$`Genes with the same annotation`$NRED
  )

  random.genes.hexb <- hexbin(bkgd.individual.Zscores$zscores$`Random Genes`$PC,
    bkgd.individual.Zscores$zscores$`Random Genes`$NRED,
    ybnds = c(min(all_scores_y), max(all_scores_y),
    xbnds = c(min(all_scores_x), max(all_scores_x)
  )

  random.annotated.genes.hexb <- hexbin(bkgd.individual.Zscores$zscores$`Random Annotated Genes`$PC,
    bkgd.individual.Zscores$zscores$`Random Annotated Genes`$NRED,
    ybnds = c(min(all_scores_y), max(all_scores_y),
    xbnds = c(min(all_scores_x), max(all_scores_x)
  )

  random.identical.annotated.genes.hexb <- hexbin(bkgd.individual.Zscores$zscores$`Genes with the same annotation`$PC,
    bkgd.individual.Zscores$zscores$`Genes with the same annotation`$NRED,
    ybnds = c(min(all_scores_y), max(all_scores_y),
    xbnds = c(min(all_scores_x), max(all_scores_x)
  )

  cnt.max <- max(c(
    random.identical.annotated.genes.hexb@count,
    random.annotated.genes.hexb@count
  ))

  plot(random.genes.hexb,
  colramp = rf, mincnt = 1, maxcnt = cnt.max,
  xlab = "PC", ylab = "NRED", main = "Random Genes"
  )

  plot(random.annotated.genes.hexb,
    colramp = rf, mincnt = 1, maxcnt = cnt.max,
    xlab = "PC", ylab = "NRED", main = "Random Annotated Genes"
  )

  plot(random.identical.annotated.genes.hexb,
    colramp = rf, mincnt = 1, maxcnt = cnt.max,
    xlab = "PC", ylab = "NRED", main = "Genes with the same annotation"
  )
}

Plot_Metric_Comparison <- function(bkgd.individual) {
  t_test_KO_random_pearson <- t.test(bkgd.individual$`Random Annotated Genes`$PC,
    bkgd.individual$`Random Genes`$PC,
    alternative = "less"
  ) # x > y (NULL)

  t_test_KO_random_euclidean <- t.test(bkgd.individual$`Random Annotated Genes`$NRED,
    bkgd.individual$`Random Genes`$NRED,
    alternative = "greater"
  ) # x > y (NULL)

  par(
    mfrow = c(2, 2),
    mar = c(3, 3, 3, 1)
  )
  # plot 1
  plot(density(bkgd.individual$`Random Genes`$PC,
    adjust = 2,
    na.rm = TRUE
  ),
  ylim = c(0, 1),
  xlab = "",
  ylab = "",
  main = ""
  )
  points(density(bkgd.individual$`Random Annotated Genes`$PC, adjust = 2),
    typ = "l",
    col = "blue"
  )
  mtext(paste("p-value = ", signif(t_test_KO_random_pearson$p.value, 2)),
    side = 3,
    col = "blue",
    padj = 2,
    cex = 0.75
  )
  title(
    ylab = "Density",
    line = 2,
    cex.lab = 1
  )
  title(
    xlab = "PC",
    line = 2,
    cex.lab = 1
  )

  # plot 2
  plot(density(bkgd.individual$`Random Annotated Genes`$PC,
    adjust = 2
  ),
  ylim = c(0, 1),
  xlab = "",
  ylab = "",
  main = " "
  )
  points(density(bkgd.individual$`Genes with the same annotation`$PC,
    adjust = 2
  ),
  typ = "l",
  col = "red"
  )
  mtext(paste("p-value = ", signif(t_test_KO_random_pearson$p.value, 2)),
    side = 3,
    col = "red",
    padj = 2,
    cex = 0.75
  )
  title(
    ylab = "Density",
    line = 2,
    cex.lab = 1
  )
  title(
    xlab = "PC",
    line = 2,
    cex.lab = 1
  )

  # plot 3
  plot(density(bkgd.individual$`Random Genes`$NRED,
    adjust = 2
  ),
  typ = "l",
  ylim = c(0, 1.25),
  xlab = "",
  ylab = "",
  main = ""
  )
  points(density(bkgd.individual$`Random Annotated Genes`$NRED, adjust = 2),
    typ = "l",
    col = "blue"
  )
  title(
    ylab = "Density",
    line = 2,
    cex.lab = 1
  )
  title(xlab = "NRED", line = 2, cex.lab = 1)
  mtext(paste("p-value = ", signif(t_test_KO_random_euclidean$p.value, 2)),
    side = 3, col = "blue", padj = 2, cex = .75
  )

  # plot 4
  plot(density(bkgd.individual$`Random Annotated Genes`$NRED,
    adjust = 2
  ),
  typ = "l",
  ylim = c(0, 1.25),
  xlab = "",
  ylab = "",
  main = ""
  )
  points(density(bkgd.individual$`Random Annotated Genes`$NRED,
    adjust = 2
  ),
  typ = "l",
  col = "red"
  )
  title(
    ylab = "Density",
    line = 2,
    cex.lab = 1
  )
  title(
    xlab = "NRED",
    line = 2,
    cex.lab = 1
  )
  title(" \n\nComparison of random & functional \n pairwise comparisons",
    outer = TRUE
  )
  mtext(paste(
    "p-value = ",
    signif(t_test_KO_random_euclidean$p.value, 2)
  ),
  side = 3,
  col = "red",
  padj = 2,
  cex = .75
  )
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
Plot_Background_Modules <- function(bkgd.traits) {
  if (length(bkgd.traits) < 1) {
    return("There are no random background distributions calculated")
  }
  colours <- rainbow(length(bkgd.traits))
  bkgd.names <- names(bkgd.traits)
  legend <- sapply(
    bkgd.names,
    function(x) paste("N = ", x, sep = "")
  )

  plot(density(bkgd.traits[[1]],
    na.rm = TRUE
  ),
  ylim = c(0, 1),
  xlim = c(-4, 4),
  ylab = "Density",
  xlab = "Composite Zscore",
  cex.main = 0.75,
  col = colours[1],
  main = "Random background distributions of modules"
  )

  for (i in 2:length(bkgd.traits)) {
    points(density(bkgd.traits[[i]],
      na.rm = TRUE
    ),
    type = "l",
    col = colours[i]
    )
  }

  legend("topleft",
    legend = legend,
    col = colours,
    pch = 1,
    cex = 0.7
  )
}

#' Calculate Association Rules
#' @name Calculate Association Rules
#' @description Calculate Association Rules using the apriori algorithm (as
#' implemented in the arules package). This function lets the user decide on the
#' number of rules to produce and estimates parameters to produce those results.
#' @param sbs.trait.attributes Matrix where each collumn is a trait attribute
#' and each row a genome. Presence or absence is indicated with either a one
#' (1) or zero (0).
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
                              lhs, rhs, N) {
  require(arules)

  .Reach_N <- function(N, lhs = c(), rhs = c()) {
    n <- 0
    support <- 1.0
    confidence <- 1.0
    while (n < N) {
      apri <- apriori(sbs.trait.attributes,
        parameter = list(
          support = support,
          confidence = confidence,
          minlen = 2
        ),
        appearance = list(
          lhs = lhs,
          rhs = rhs
        )
      )

      rules <- as(apri, "data.frame")
      n <- nrow(rules)
      print(n)
      support <- support - 0.1
      if (support <= 0) {
        support <- 1.0
        confidence <- confidence - 0.1
      }
    }
    return(rules)
  }



  if (missing(lhs) && missing(rhs)) {
    rules <- .Reach_N(N)
  } else if (!missing(lhs) && !missing(rhs)) {
    rules1 <- .Reach_N(N / 2, lhs = lhs)
    rules2 <- .Reach_N(N / 2, rhs = rhs)

    rules1 <- rules1[order(rules1[, 1]), ]
    rules2 <- rules2[order(rules2[, 1]), ]

    rules <- rbind(rules1[1:(N / 2), ], rules2[1:(N / 2), ])
  } else if (!missing(lhs)) {
    rules <- .Reach_N(N, lhs = lhs)
  } else if (!missing(rhs)) {
    rules <- .Reach_N(N, rhs = rhs)
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
  y.max <- 140
  x.max <- 60


  attributes <- trait.attributes.pruned[[trait]]
  n.attributes <- length(attributes)
  genes <- RNAseq.data$features$annotation.db$module.dict[[trait]]
  n.genes <- length(genes)

  plot(c(0, as.numeric(x.max)), c(0, y.max), type = "n", xlab = "", ylab = "", axes = F)

  y.coords <- seq(y.max, 1, by = -15) # With this setting max of 7 genes
  x.coords <- seq(0, x.max, by = 11)

  # For each gene
  for (i in 1:n.genes) {
    gene <- genes[i]

    # Gene name
    graphics::text(
      x = (x.coords[1] + 5),
      y = (y.coords[i] - 5),
      labels = gene
    )

    # for each attribute
    for (y in 1:n.attributes) {
      #      Draw Labels
      graphics::text(x = (x.coords[y + 1] + 5), y = y.max, labels = paste(trait,
        y,
        sep = "."
      ))

      rect(
        xleft = x.coords[y + 1], ybottom = y.coords[i] - 10,
        xright = x.coords[y + 1] + 10, ytop = y.coords[i]
      )

      # Draw expression line
      genomes <- attributes[[y]]$genomes
      genomes.rna <- RNAseq.data$table[which(RNAseq.data$table$Bin %in% genomes), ]
      genomes.rna.gene <- genomes.rna[
        which(genomes.rna$Annotation == gene),
        RNAseq.data$features$rank.columns
      ]
      genomes.rna.gene.mean <- apply(genomes.rna.gene, 2, mean)

      # Generate box colours
      values <- seq(0.1, 1, by = 0.1)
      cols <- colorRamp(c("white", "red"))
      colours <- rgb(cols(values) / 255)

      # Draw all the boxes!
      box.coords <- seq(0, 10, length.out = (length(genomes.rna.gene.mean) + 1))
      for (z in 2:length(box.coords)) {
        rank <- genomes.rna.gene.mean[z - 1]
        col <- "red"
        if (is.nan(rank) || is.na(rank)) {
          rank <- 1
        }
        # I am not proud of this
        if (rank <= 1 && rank > 0.9) {
          col <- colours[1]
        }
        else if (rank <= 0.9 && rank > 0.8) {
          col <- colours[2]
        } else if (rank <= 0.8 && rank > 0.7) {
          col <- colours[3]
        } else if (rank <= 0.7 && rank > 0.6) {
          col <- colours[4]
        } else if (rank <= 0.6 && rank > 0.5) {
          col <- colours[5]
        } else if (rank <= 0.5 && rank > 0.4) {
          col <- colours[6]
        } else if (rank <= 0.4 && rank > 0.3) {
          col <- colours[7]
        } else if (rank <= 0.3 && rank > 0.2) {
          col <- colours[8]
        } else if (rank <= 0.2 && rank > 0.1) {
          col <- colours[9]
        } else if (rank <= 0.1 && rank >= 0.0) {
          col <- colours[10]
        }
        rect(
          xleft = x.coords[y + 1] + box.coords[(z - 1)], ybottom = y.coords[i] - 10,
          xright = x.coords[y + 1] + box.coords[z], ytop = y.coords[i],
          col = col
        )
        for (line.idx in 1:nrow(genomes.rna.gene)) {
          y.bot <- y.coords[i] - 10
          rank.pos <- (genomes.rna.gene[line.idx, ] * 10) + y.bot
          x.pos <- seq(1, 10, length.out = 6) + x.coords[y + 1]
          graphics::lines(
            x.pos,
            rank.pos
          )
        }
      }
    }
  }
}

#' Get module categories
#' @description Returns a list of the kegg categories with their corresponding modules
#' @export
Get_module_categories <- function(pythonfile){
  reticulate::source_python(pythonfile)
  parsed <- get_module_categories(kegg_categories)

  cat.genes <- parsed[[1]]
  names(cat.genes) <- parsed[[2]]

  return(cat.genes)
}

#' Get Metric Dist
#' @param metric either {'PC','NRED'}
#' @param cat.genes list("Carbon")
getMetricDist <- function(metric, cat.genes, distance.metrics, RNAseq.data,
                          bkgd.individual, bkgd.individual.Zscores) {
  mean.distances <- matrix(nrow = 0, ncol = 4)

  for (category in names(cat.genes)) {
    distances.list <- list()
    cat(paste(category, "\n"))

    for (KO in cat.genes[[category]]) {
      data.rows <- RNAseq.data$table[which(RNAseq.data$table$Annotation == KO), ]

      genomes.present <- unique(data.rows$Bin)
      if (length(genomes.present) <= 1) {
        next()
      }
      dist.matrix <- matrix(data = NA, ncol = length(genomes.present), nrow = length(genomes.present))
      for (x in 1:(length(genomes.present) - 1)) {
        gen.x <- data.rows[which(data.rows$Bin == genomes.present[x]), ][1, ]

        for (y in (x + 1):length(genomes.present)) {
          gen.y <- data.rows[which(data.rows$Bin == genomes.present[y]), ][1, ]
          dist.matrix[x, y] <- distance.metrics[[metric]](gen.x, gen.y, RNAseq.data$features)
        }
      }

      distances.list[[KO]] <- dist.matrix
    } # ko

    distances <- unlist(distances.list)

    p.val <- NA
    p.c <- "not significant"

    try(p.val <- t.test(distances, bkgd.individual$`Random Annotated Genes`[[metric]],
      na.action = na.omit
    )$p.value,
    silent = T
    )

    try(if (p.val <= 0.05) {
      p.c <- "significant"
    } else {
      p.c. <- "not significant"
    }
    ,
    silent = T
    )


    cex <- median(distances, na.rm = T)
    if (is.null(cex)) {
      cex <- NA
    }


    # actual distances
    mean.dist <- mean(distances, na.rm = T)
    mean.dist.Z <- (mean.dist - bkgd.individual.Zscores$mu$`Random Annotated Genes`[[metric]]) /
      bkgd.individual.Zscores$sd$`Random Annotated Genes`[[metric]]



    mean.distances <- rbind(
      mean.distances,
      c(mean.dist.Z, p.c, category, cex)
    )
  }

  return(mean.distances)
}

#' Get Metric Dist module
#' @export
getMetricDistModule <- function(metric, cat.modules) {
  mean.distances <- matrix(nrow = 0, ncol = 5)

  for (category in names(cat.modules)) {
    print(category)
    for (module in cat.modules[[category]]) {
      distances.list <- list()

      for (KO in RNAseq.data$features$annotation.db$module.dict[[module]]) {
        data.rows <- RNAseq.data$table[which(RNAseq.data$table$Annotation == KO), ]

        genomes.present <- unique(data.rows$Bin)
        if (length(genomes.present) <= 1) {
          next()
        }
        dist.matrix <- matrix(data = NA, ncol = length(genomes.present), nrow = length(genomes.present))
        for (x in 1:(length(genomes.present) - 1)) {
          gen.x <- data.rows[which(data.rows$Bin == genomes.present[x]), ][1, ]

          for (y in (x + 1):length(genomes.present)) {
            gen.y <- data.rows[which(data.rows$Bin == genomes.present[y]), ][1, ]
            dist.matrix[x, y] <- distance.metrics[[metric]](gen.x, gen.y, RNAseq.data$features)
          }
        }

        distances.list[[KO]] <- dist.matrix
      } # KOs

      distances <- unlist(distances.list)
      distances <- as.numeric(distances[-which(is.na(distances))])

      p.val <- NA
      p.c <- "not significant"
      n.ko <- as.character(length(
        RNAseq.data$features$annotation.db$module.dict[[module]]
      ))

      try(p.val <- t.test(distances, bkgd.traits[[n.ko]])$p.value,
        silent = T
      )

      try(if (p.val <= 0.05) {
        p.c <- "significant"
      } else {
        p.c. <- "not significant"
      }
      ,
      silent = T
      )
      # actual distances
      mean.dist <- median(distances, na.rm = T)
      mean.dist.Z <- (mean.dist - bkgd.individual.Zscores$mu$`Random Annotated Genes`[[metric]]) /
        bkgd.individual.Zscores$sd$`Random Annotated Genes`[[metric]]
      mean.distances <- rbind(
        mean.distances,
        c(mean.dist.Z, p.val, p.c, category, module)
      )
    } # module
  }

  return(mean.distances)
}

#' @export
Plot_Pathway_genes <- function(metric_name, distance.metrics, RNAseq.data,
                               bkgd.individual, bkgd.individual.Zscores) {
  # Catogerized Modules
  cat.modules <- read.csv("/home/joris/categorized_modules2.csv", sep = ";", header = F)

  # Hacking categories because the order is known ----
  categories <- c(
    rep("P_Energy metabolism", 4),
    rep("P_Carbohydrate and lipid metabolism", 10),
    rep("P_Nucleotide and AA metabolism", 12),
    rep("P_Secondary metabolism", 2),
    rep("S_Energy metabolism", 2),
    rep("S_Genetic information processing", 10),
    rep("S_Environmental information processing", 9),
    rep("F_Metabolism", 2),
    rep("F_Environemntal information processing", 2),
    rep("F_Cellular processes", 1),
    rep("Si_Gene set", 4)
  )
  #
  cat.genes <- list()

  # Do this better in the future
  con <- file("/home/joris/categorized_modules2.csv", "r")
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
      break
    }
    line <- unlist(strsplit(line, ";"))

    cat.genes[[line[1]]] <- unique(line[2:length(line)])
  }

  close(con)


  # Get Zscores ----
  metric.Z <- getMetricDist(
    metric_name, cat.genes, distance.metrics, RNAseq.data,
    bkgd.individual, bkgd.individual.Zscores
  )

  library(ggplot2)
  ids <- 1:nrow(metric.Z)
  metric.Z <- cbind(metric.Z, ids)
  metric.Z <- cbind(metric.Z, categories)

  colnames(metric.Z) <- c("value", "sig", "collection", "cex", "ids", "categories")


  metric.Z <- as.data.frame(metric.Z, stringsAsFactors = F)
  metric.Z$value <- as.numeric(metric.Z$value)
  metric.Z$ids <- as.numeric(metric.Z$ids)
  metric.Z$ids <- factor(metric.Z$ids, levels = metric.Z$ids)

  metric.Z$cex <- exp(as.numeric(metric.Z$cex) / mean(as.numeric(metric.Z$cex),
    na.rm =
    ))

  metric.plot <- ggplot(
    metric.Z, aes(
      x = ids, y = value,
      colour = categories,
      shape = factor(sig)
    )
    # ,size = cex
  ) +
    geom_point(size = 3) +
    xlab("") +
    ylab(paste(metric_name, "Z score", sep = " ")) +
    facet_grid(. ~ categories, scales = "free_x", space = "free_x") +
    geom_hline(yintercept = 0) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )


  # Combine plot
  # plot(gridExtra::arrangeGrob(pearson.plot, nred.plot))

  return(list(
    "data" = metric.Z,
    "plot" = metric.plot
  ))
}

Plot_Pathway_modules <- function() {


  sig.pathways <- nred.Z[which(nred.Z$sig == "significant"), "collection"]
  cat.modules.sig <- cat.modules[sig.pathways]


  nred.modules.Z <- getMetricDistModule("NRED", cat.modules.sig)
  nred.modules.Z <- cbind(nred.modules.Z, (1:nrow(nred.modules.Z)))
  colnames(nred.modules.Z) <- c("value", "pval", "sig", "pathway", "module", "ids")
  nred.modules.Z.df <- as.data.frame(nred.modules.Z, stringsAsFactors = F)
  nred.modules.Z.df$value <- as.numeric(nred.modules.Z.df$value)
  nred.modules.Z.df$ids <- as.numeric(nred.modules.Z.df$ids)
  nred.modules.Z.df$ids <- factor(nred.modules.Z.df$ids, levels = nred.modules.Z.df$ids)

  nred.modules.Z.sig <- nred.modules.Z.df[which(nred.modules.Z.df$sig == "significant"), ]

  ggplot(nred.modules.Z.sig, aes(x = ids, y = value, colour = pathway, shape = factor(sig))) +
    geom_point(size = 2) +
    xlab("") +
    ylab("NRED Z score") +
    facet_grid(. ~ pathway, scales = "free_x", space = "free_x") +
    geom_hline(yintercept = 0) +
    guides(colour = FALSE)


  # for each pathway the % of sig
  percentage_sig <- list()
  for (pathway in unique(nred.modules.Z.df$pathway)) {
    pathway_rows <- nred.modules.Z.df[which(
      nred.modules.Z.df$pathway == pathway
    ), ]

    total <- nrow(pathway_rows)
    sig <- nrow(pathway_rows[which(pathway_rows$sig == "significant"), ])
    percentage_sig[[pathway]] <- ((sig * 100) / total)
    cat(paste(pathway, ((sig * 100) / total), sep = ";"), "\n")
  }


  pc.modules.Z <- getMetricDistModule("PC", cat.modules.sig)
  pc.modules.Z <- cbind(pc.modules.Z, (1:nrow(pc.modules.Z)))
  colnames(pc.modules.Z) <- c("value", "pval", "sig", "pathway", "module", "ids")
  pc.modules.Z.df <- as.data.frame(pc.modules.Z, stringsAsFactors = F)
  pc.modules.Z.df$value <- as.numeric(pc.modules.Z.df$value)
  pc.modules.Z.df$ids <- as.numeric(pc.modules.Z.df$ids)
  pc.modules.Z.df$ids <- factor(pc.modules.Z.df$ids, levels = pc.modules.Z.df$ids)

  pc.modules.Z.sig <- pc.modules.Z.df[which(pc.modules.Z.df$sig == "significant"), ]


  sigboth <- pc.modules.Z.sig[which(pc.modules.Z.sig$module %in%
    nred.modules.Z.sig$module), ]
  sigboth.pruned <- sigboth[which(sigboth$module %in% prune_lalala(sigboth)), ]
}


Plot_Trait_Expression <- function(trait, subset_genomes) {
  dev.off()
  .getRank <- function(genome, rows) {
    genome.rows <- rows[
      which(rows$Bin == genome),
      RNAseq.data$features$rank.columns
    ]

    x <- genome.rows[sample.int(nrow(genome.rows), 1), ]

    return(x)
  }


  annotations <- RNAseq.data$features$annotation.db$module.dict[[trait]]
  par(mfrow = calcmfrow(length(annotations)))



  for (KO in annotations) {
    print(KO)
    rows <- RNAseq.data$table[which(RNAseq.data$table$Annotation == KO), ]

    if (missing(subset_genomes)) {
      genomes <- unique(rows$Bin)
    } else {
      genomes <- unique(rows$Bin)
      genomes <- subset_genomes[which(subset_genomes %in% genomes)]
    }

    print(genomes)
    if (length(genomes) == 0) {
      next()
    }

    plot(1:length(RNAseq.data$features$rank.columns), .getRank(genomes[1], rows),
      ylim = c(0, 1), main = KO, type = "l"
    )
    for (genome in genomes) {
      lines(1:length(RNAseq.data$features$rank.columns), .getRank(genome, rows))
    }
  }
}

Plot_Venn <- function(trait.attributes.pruned) {
  tap <- trait.attributes.pruned
  b16 <- c()
  b39 <- c()
  total <- c()

  for (trait.idx in 1:length(tap)) {
    t <- tap[[trait.idx]]
    #####
    if (length(t) == 0) {
      next
    }
    #####

    for (ta.idx in 1:length(t)) {
      ta <- t[[ta.idx]]
      ta.name <- paste(names(tap)[trait.idx],
        ta.idx,
        sep = "."
      )
      total <- c(total, ta.name)
      genomes <- ta$genomes

      c <- F

      if ("16" %in% genomes) {
        b16 <- c(b16, ta.name)
        c <- T
      }

      if ("39" %in% genomes) {
        b39 <- c(b39, ta.name)
        c <- T
      }

      # if (! c) {
      #   total <- c(total, ta.name)
      # }
    }
  }

  o12 <- total[which(total %in% b39)]
  o23 <- b39  [which(b39 %in% b16)]
  o13 <- total[which(total %in% b16)]
  o123 <- o12  [which(o12 %in% b16)]

  VennDiagram::draw.triple.venn(length(total),
    length(b39),
    length(b16),
    length(o12),
    length(o23),
    length(o13),
    length(o123),
    category = c("Total", "b39", "b16")
  )

  .getGenomeTraits <- function(genome) {
    tas <- c()

    for (trait.idx in 1:length(tap)) {
      t <- tap[[trait.idx]]
      #####
      if (length(t) == 0) {
        next
      }
      #####

      for (ta.idx in 1:length(t)) {
        ta <- t[[ta.idx]]
        ta.name <- paste(names(tap)[trait.idx],
          ta.idx,
          sep = "."
        )

        genomes <- ta$genomes

        c <- F

        if (genome %in% genomes) {
          tas <- c(tas, ta.name)
          c <- T
        }
      }
    }
    return(tas)
  }
  library(RCy3)
  nodes <- as.matrix(RNAseq.data$features$bins, ncol = 1)
  colnames(nodes) <- "id"

  edges <- matrix(ncol = 3, nrow = 0)
  for (x in 1:(length(RNAseq.data$features$bins) - 1)) {
    x.tas <- .getGenomeTraits(RNAseq.data$features$bins[x])
    for (y in (x + 1):length(RNAseq.data$features$bins)) {
      y.tas <- .getGenomeTraits(RNAseq.data$features$bins[y])
      overlap <- length(x.tas[which(x.tas %in% y.tas)])
      edges <- rbind(edges, c(
        RNAseq.data$features$bins[x],
        RNAseq.data$features$bins[y],
        overlap
      ))
    }
  }
  colnames(edges) <- c("source", "target", "weight")
  weights <- as.numeric(edges[, "weight"])
  weights <- (weights / min(weights))^2

  edges[, 3] <- weights

  createNetworkFromDataFrames(data.frame(nodes, stringsAsFactors = F),
    data.frame(edges, stringsAsFactors = F),
    title = "title",
    collection = "Network Traits and Genomes"
  )
}

Shared_traits <- function( RNAseq.data ) {

  interesting_bins <- c('16','39','22','29', '32')
  other_bins <- RNAseq.data$features$bins[which(!(RNAseq.data$features$bins %in% interesting_bins))]
  rows_shared <- RNAseq.data$features$trait_presence_absence %>% apply(1, function(row){
    if (sum(row[interesting_bins]) == length(interesting_bins) && sum(row[other_bins]) == 0) {
      return(T)
    } else{F}
  })

  shared_traits <- rows_shared[which(rows_shared)] %>% names


  distances <- pairwise.distances[RNAseq.data$features$annotation.db$module.dict$M00150] %>% sapply(function(x) {
    x[c("16","39"),c("16","39")]
  }) %>% .[which(!is.na(.))]

  t.test(distances, bkgd.traits$`4`, alternative = 'less')


}

Plot_Shared_Attributes <- function(trait.attributes.pruned, RNAseq.data) {
  library(magrittr)
  trait.names <- list()

  accum_bins <- c('16',"39")

  for (trait.name in names(trait.attributes.pruned)) {
    trait <- trait.attributes.pruned[[trait.name]]

    for (attribute.name in names(trait)) {
      attribute <- trait[[attribute.name]]
      genomes <- attribute$genomes
      if (sum(accum_bins %in% genomes) == 2) {
        trait.names[[trait.name]] <- attribute.name
      }
    }
  }

  other_bins <- RNAseq.data$features$bins[-which(RNAseq.data$features$bins %in% accum_bins)]

  bin_occurences <- sapply(other_bins, function(bin) {
    occurences <- sapply(1:length(trait.names), function(i) {
      ta <- trait.attributes.pruned[[names(trait.names)[i]]][[trait.names[[i]]]]
      if (bin %in% ta$genomes) {
        return(1)
      }
      else {
        return(0)
      }
    })

    return(sum(occurences))
  }) %>%
    sort(decreasing = T) %T>%
    plot(
      type = "l", xaxt = "n",
      xlab = "Genomic bins",
      ylab = "Number of attributes shared with CAA",
      ylim = c(0,30)
    )
  axis(1, at = 1:length(bin_occurences), labels = names(bin_occurences))
}

Plot_Redundancy_Traits <- function(RNAseq.data) {
  library(ggplot2)
  library(magrittr)

  ta.pa <- apply(RNAseq.data$features$trait_presence_absence, 1, function(row) {
    return(sum(row))
  })

  ta.pa <- ta.pa[which(ta.pa != 0)]
  sort(ta.pa) %>% barplot(xaxt='n', ylim = c(0,20))

}

#' @export
Calc_TnA_redundancy <- function(RNAseq.data ) {
  point.matrix <- matrix(ncol=4, nrow=0)

  t.pa <- RNAseq.data$features$trait_presence_absence

  combinations <- combn(RNAseq.data$features$bins, 2, simplify = F)

  for (pair in combinations){
    try({
      A <- pair[1] %>% as.character
      B <- pair[2] %>% as.character

      traits.A <- t.pa[, A] %>% which(. == T) %>% names
      traits.B <- t.pa[, B] %>% which(. == T) %>% names
      overlap.traits <- intersect(traits.A, traits.B) %>% length

      attributes.A <- which(sbs.trait.attributes[A,]  == 1) %>% names
      print(B)
      attributes.B <- which(sbs.trait.attributes[B,]  == 1) %>% names
      overlap.attributes <- intersect(attributes.A, attributes.B) %>% length

      point.matrix <- rbind(point.matrix,
                            c(A, B,
                              overlap.traits,
                              overlap.attributes))
    })

  }
  colnames(point.matrix) <- c("A", "B", 'traits', 'attributes')

  point.df <- data.frame(point.matrix, stringsAsFactors = F)
  point.df$traits %<>% as.numeric
  point.df$attributes %<>% as.numeric
  all_bins <- unique(c(point.df[, 1], point.df[, 2]))
  sum_traits <- NULL
  sum_attributes <- NULL

  for (i in all_bins) {
    bin_rows <- which(point.df[, 1] == i | point.df[, 2] == i)
    sum_traits <- c(sum_traits, sum(point.df[bin_rows, 3]))
    sum_attributes <- c(sum_attributes, sum(point.df[bin_rows, 4]))
  }

  return(cbind(all_bins, sum_traits, sum_attributes))
}

Plot_traits_vs_attributes <- function(model_bin) {
  point.matrix <- matrix(ncol=3, nrow=0)
  point.matrix_39 <- matrix(ncol=3, nrow=0)

  t.pa <- RNAseq.data$features$trait_presence_absence

  combinations <- combn(RNAseq.data$features$bins, 2, simplify = F)

  if(!missing(model_bin)){
    combinations_39 <- lapply(RNAseq.data$features$bins, function(bin){
      if(!bin==model_bin){
      return(c(bin, model_bin))
      }
    })
  }


  for (pair in combinations){
    try({
      A <- pair[1] %>% as.character
      B <- pair[2] %>% as.character

      traits.A <- t.pa[, A] %>% which(. == T) %>% names
      traits.B <- t.pa[, B] %>% which(. == T) %>% names
      overlap.traits <- intersect(traits.A, traits.B) %>% length

      attributes.A <- which(sbs.trait.attributes[A,]  == 1) %>% names
      print(B)
      attributes.B <- which(sbs.trait.attributes[B,]  == 1) %>% names
      overlap.attributes <- intersect(attributes.A, attributes.B) %>% length

      vs <- paste(pair, collapse = 'vs')

      point.matrix <- rbind(point.matrix,
                            c(vs,
                              overlap.traits,
                              overlap.attributes))
    })

  }
  colnames(point.matrix) <- c("pairs", 'traits', 'attributes')

  point.df <- data.frame(point.matrix, stringsAsFactors = F)
  point.df$traits %<>% as.numeric
  point.df$attributes %<>% as.numeric
  ylim_plot1 <- c(0,max(point.df$attributes)+5)
  xlim_plot1 <- c(0, max(point.df$traits)+5)

  plot(x = as.numeric(point.matrix[,2]),
       y = as.numeric(point.matrix[,3]),
       xlab = '# Overlap Traits',
       ylab = '# Overlap Attributes',
       main = "Pairwise Genome Comparisons",
       ylim = ylim_plot1,
       xlim = xlim_plot1)

  model <- lm(attributes~traits, data = point.df)
  abline(model$coefficients)

  cinterval <- confint(model, level=.99)
  abline(cinterval[,1])
  abline(cinterval[,2])


  if(!missing(model_bin)){
  for (pair in combinations_39){
    try({
      A <- pair[1] %>% as.character
      B <- pair[2] %>% as.character

      traits.A <- t.pa[, A] %>% which(. == T) %>% names
      traits.B <- t.pa[, B] %>% which(. == T) %>% names
      overlap.traits <- intersect(traits.A, traits.B) %>% length

      attributes.A <- which(sbs.trait.attributes[A,]  == 1) %>% names
      print(B)
      attributes.B <- which(sbs.trait.attributes[B,]  == 1) %>% names
      overlap.attributes <- intersect(attributes.A, attributes.B) %>% length

      vs <- paste(pair, collapse = 'vs')

      point.matrix_39 <- rbind(point.matrix_39,
                            c(vs,
                              overlap.traits,
                              overlap.attributes))
    })

  }
  colnames(point.matrix_39) <- c("pairs", 'traits', 'attributes')

  point.df <- data.frame(point.matrix_39, stringsAsFactors = F)
  point.df$traits %<>% as.numeric
  point.df$attributes %<>% as.numeric


  points(x = as.numeric(point.matrix_39[,2]),
       y = as.numeric(point.matrix_39[,3]),
       xlab = '# Overlap Traits',
       ylab = '# Overlap Attributes',
       main = "Pairwise Genome Comparisons",
       ylim = ylim_plot1,
       xlim = xlim_plot1,
       pch=19,
       cex=1.5,
       col="red")

  text(x=10,y=ylim_plot1[2]-5, paste("bin_",model_bin))
#  model_39 <- lm(attributes~traits, data = point.df)
#  abline(model_39$coefficients, col="red")

#  cinterval_39 <- confint(model_39, level=.99)
#  abline(cinterval_39[,1], col="red")
#  abline(cinterval_39[,2], col="red")
  #polygon

  }
  # identify(x= as.numeric(point.matrix[,2]), y = as.numeric(point.matrix[,3]), labels = point.matrix[,1])
}

#' @export
Plot_traits_vs_attributes_highlight <- function(model_bin, string_to_colors) {
  point.matrix <- matrix(ncol=3, nrow=0)
  point.matrix_39 <- matrix(ncol=3, nrow=0)

  t.pa <- RNAseq.data$features$trait_presence_absence

  combinations <- combn(RNAseq.data$features$bins, 2, simplify = F)

  if(!missing(model_bin)){
    combinations_39 <- lapply(RNAseq.data$features$bins, function(bin){
      if(!bin==model_bin){
        return(c(bin, model_bin))
      }
    })
  }


  for (pair in combinations){
    try({
      A <- pair[1] %>% as.character
      B <- pair[2] %>% as.character

      traits.A <- t.pa[, A] %>% which(. == T) %>% names
      traits.B <- t.pa[, B] %>% which(. == T) %>% names
      overlap.traits <- intersect(traits.A, traits.B) %>% length

      attributes.A <- which(sbs.trait.attributes[A,]  == 1) %>% names
      print(B)
      attributes.B <- which(sbs.trait.attributes[B,]  == 1) %>% names
      overlap.attributes <- intersect(attributes.A, attributes.B) %>% length

      vs <- paste(pair, collapse = 'vs')

      point.matrix <- rbind(point.matrix,
                            c(vs,
                              overlap.traits,
                              overlap.attributes))
    })

  }
  colnames(point.matrix) <- c("pairs", 'traits', 'attributes')

  point.df <- data.frame(point.matrix, stringsAsFactors = F)
  point.df$traits %<>% as.numeric
  point.df$attributes %<>% as.numeric
  ylim_plot1 <- c(0,max(point.df$attributes)+5)
  xlim_plot1 <- c(0, max(point.df$traits)+5)

  plot(x = as.numeric(point.matrix[,2]),
       y = as.numeric(point.matrix[,3]),
       xlab = '# Overlap Traits',
       ylab = '# Overlap Attributes',
       main = "Pairwise Genome Comparisons",
       ylim = ylim_plot1,
       xlim = xlim_plot1)

  model <- lm(attributes~traits, data = point.df)
  abline(model$coefficients)

  cinterval <- confint(model, level=.99)
  abline(cinterval[,1])
  abline(cinterval[,2])


  if(!missing(model_bin)){
    for (pair in combinations_39){
      try({
        A <- pair[1] %>% as.character
        B <- pair[2] %>% as.character

        traits.A <- t.pa[, A] %>% which(. == T) %>% names
        traits.B <- t.pa[, B] %>% which(. == T) %>% names
        overlap.traits <- intersect(traits.A, traits.B) %>% length

        attributes.A <- which(sbs.trait.attributes[A,]  == 1) %>% names
        print(B)
        attributes.B <- which(sbs.trait.attributes[B,]  == 1) %>% names
        overlap.attributes <- intersect(attributes.A, attributes.B) %>% length

        vs <- paste(pair, collapse = 'vs')

        point.matrix_39 <- rbind(point.matrix_39,
                                 c(vs,
                                   overlap.traits,
                                   overlap.attributes))
      })

    }
    colnames(point.matrix_39) <- c("pairs", 'traits', 'attributes')

    point.df <- data.frame(point.matrix_39, stringsAsFactors = F)
    point.df$traits %<>% as.numeric
    point.df$attributes %<>% as.numeric


    points(x = as.numeric(point.matrix_39[,2]),
           y = as.numeric(point.matrix_39[,3]),
           xlab = '# Overlap Traits',
           ylab = '# Overlap Attributes',
           main = "Pairwise Genome Comparisons",
           ylim = ylim_plot1,
           xlim = xlim_plot1,
           pch=19,
           cex=1.5,
           col=string_to_colors)

    text(x=10,y=ylim_plot1[2]-5, paste("bin_",model_bin))
    #  model_39 <- lm(attributes~traits, data = point.df)
    #  abline(model_39$coefficients, col="red")

    #  cinterval_39 <- confint(model_39, level=.99)
    #  abline(cinterval_39[,1], col="red")
    #  abline(cinterval_39[,2], col="red")
    #polygon

  }
  # identify(x= as.numeric(point.matrix[,2]), y = as.numeric(point.matrix[,3]), labels = point.matrix[,1])
}

Get_KEGG_Sub_Modules <- function(){
  library(reticulate)
  library(magrittr)

  source_python('/home/joris/TcT/python/parse_module_definition.py')

  # module_names <- RNAseq.data$features$annotation.db$module.dict %>% names

  start <- Sys.time()
  sub_modules <- lapply(module_names, function(module) {
    print(paste0("Parsing: ", module))

  if (length(grep('M[0-9]{5}', module)) != 0 ) {
    return(parse_module(module) )
  }else{
   return(list(RNAseq.data$features$annotation.db$module.dict[[module]]))
  }

  })
  end <- Sys.time()

  names(sub_modules) <- module_names

}

#' Identify whether an annotation term has a certain expression pattern in one
#' or more genomes
#' @name Go Fish
#' @description Identify whether a gene has a certain expression pattern in one or more genomes.
#' @param RNAseq.data Collection of multple components, include RNA seq data, annotations, etc. See \code{\link{Pre_process_input}} for the full list.
#' @export
#' @example Go_Fish('K00927', RNAseq.data, type = 'peak')
#' @author JJM van Steenbrugge
Go_Fish <- function(RNAseq.data){

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
    cat("Pick a rank value (Y axis) for each time point by clicking in the graph (start at 1)")

    positions <- list('x' = c(),
                      'y' = c()
    )

    for(i in 1:nsamples){
      pos <- locator(1)
      positions$x <- c(positions$x, pos$x)
      positions$y <- c(positions$y, pos$y)

      plot(positions,
           xlab = 'Time points',
           ylab = 'Normalized Rank',
           xlim = c(1, nsamples),
           ylim = c(0,1),
           type = 'b')

    }

    return(positions)

  }

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

#' Model Module
#' @name Model_Module
#' @description Creates a bar plot of pairwise comparisons between the traits of a model genome.
#' If a distance is less than the 99% quantile of the background distances for a module of similar size,
#' The bar is presented in dark grey.
#' @param trait character string with the name of a trait.
#' @param RNAseq.data Collection of multple components, include RNA seq data,
#' annotations, etc.
#' @seealso \code{\link{Pre_process_input}} for the full list of parameters,
#' \url{https://github.com/Jorisvansteenbrugge/TcT/wiki/The-RNAseq.data-object}
#' for more information.
#' @param RNAseq.data
#' @param trait.attributes
#' @param Model_Bin
#' @param Module_Pool
#' @param Bin_Order
#' @param Yrange It is optional to provide a Yrange
#' @examples Module_Model_List <- Model_Module(RNAseq.data, trait.attributes, '39', Module_Names, c(-2.5,2.5), bkgd.traits)
#'\dontrun{
#' Model_Bin <- 39
#' Bin_Order <- c(48,32,31,29,22,11,39,16,53,45,42,28,20,25,19,8,36,26,17)
#' Yrange <- c(-2.5,2.5)
#' Module_Names_ATP_synthesis <- c("M00144","M00143","M00146","M00147","M00149","M00150","M00148","M00162","M00151","M00152","M00154","M00155","M00153","M00417","M00416","M00156","M00157","M00159")
#'}
#' @export
#' @author BO Oyserman
#'
Model_Module <- function(RNAseq.data, trait.attributes, Model_Bin, Module_Names,bkgd.traits) {


  annotation.db <- RNAseq.data$features$annotation.db

  # Step one, take the module list and match it to the module dictionary psoitions.
# This is done by parsing out the module.dict
# There are various parts of this function
# 1) making a matrix of pairwise distances to the model
# 2) calculate the order on the X and Y axis
#     X - axis: the order of the genomes based on similarity to the model
#     Y - axis: the order of the modules based on their similarity to the model


  Module_Positions         <- match(Module_Names,
                                   names(annotation.db$module.dict))
  Module_Pool              <- NULL
  Model_Comparison_Matrix  <- NULL
  Module_Name_vector       <- NULL
  Model_Sig_Matrix         <- NULL
  Fish_Backgrounds_trimmed <- NULL

  # Using non-disjunctive form
  # Module_lengths           <- lapply(annotation.db$module.dict[Module_Names],
  #                                    length) %>% as.numeric

  Module_lengths           <- d_module_size_range_all[which(names(annotation.db$module.dict)%in%Module_Names)]

  Module_background_distribution_index <- match(as.numeric(Module_lengths),
                                                as.numeric(names(bkgd.traits)))

  Fish_Backgrounds <- c()
  for (i in Module_background_distribution_index){

    Fish_Background  <- quantile(unlist( bkgd.traits[i]),
                                 probs = seq(0, 0.1, .01))[6]
    Fish_Backgrounds <- c(Fish_Backgrounds,
                          Fish_Background)

  }

  # For each module name, get the number of genes in the module. This needs to be updated to use the disjunctive normalized forms.

  for (i in 1: length(Module_Names)) {
    Module_Pool[i] <- list(trait.attributes[[ Module_Positions[i] ]] [2])
  }

  # Reorder the bins based on the input variable, hashed out for updated version in which order is defined based on similarity
  # If a predetermined bin order is  supplied ...

  # Bin_Order_Index <- match(Bin_Order, rownames(Module_Pool[[1]] [[1]]))

  # Make Model Comparison Matrix, filtering modules that were not present in the model organism

  for (i in 1: length(Module_Pool) ) {

    # Create full matrix from each module by filling in the bottom of the traingle
    fullmatrix                                      <- Module_Pool[[i]] [[1]]
    fullmatrix[lower.tri(fullmatrix, diag = FALSE)] <- fullmatrix[upper.tri(fullmatrix,
                                                                            diag = FALSE)]
    # Parse out the vector of the model organisms and all pairwise comparisons
    Model_Comparison_vector <- fullmatrix[,
                                           which(dimnames(fullmatrix) [[1]] == Model_Bin
                                                  )]
    # print (Model_Comparison_vector)

    # If the vector is empty, go to the next module. Else, cbind the vector to a growing matrix containing all Module comparisons for that Model genome
    if(sum(is.na(Model_Comparison_vector)) == length(Model_Comparison_vector)){
      next
    } else {
      Model_Comparison_Matrix  <- cbind(Model_Comparison_Matrix,
                                        Model_Comparison_vector)
      Module_Name_vector       <- c(Module_Name_vector,
                                    Module_Names[i])
      # Only keep relavent backgrounds
      Fish_Backgrounds_trimmed <- c(Fish_Backgrounds_trimmed,
                                    Fish_Backgrounds[i])
    }

  }


  # Clean up the values, converting NaN to NA
  Model_Comparison_Matrix[ which(Model_Comparison_Matrix == "NaN") ] < -NA

  # Identify the bin orders their similarity with the model, and create an index
  Sum_Similarity  <- apply(Model_Comparison_Matrix, 1, sum, na.rm = TRUE)
  Bin_Order       <- Model_Comparison_Matrix %>% apply(1, sum, na.rm = T) %>%
                          sort %>% names

  Bin_Order_Index <- match(Bin_Order,
                           rownames(Model_Comparison_Matrix)) %>% rev

  # Next make a matrix of whether a module is significant
  for (i in 1: ncol(Model_Comparison_Matrix) ) {

    sig_bins         <- (Model_Comparison_Matrix[,i] <= Fish_Backgrounds_trimmed[i])
    Model_Sig_Matrix <- cbind(Model_Sig_Matrix,
                              sig_bins)
  }

  # Name the columns by module name
  colnames(Model_Comparison_Matrix) <- Module_Name_vector
  names(Fish_Backgrounds_trimmed)   <- Module_Name_vector
  colnames(Model_Sig_Matrix)        <- Module_Name_vector

  # Identify the module orders based on their similarity with the model, and create an index.
  # Similarity is based on the number of significantly similar

  Sum_Similarity_M   <- apply(Model_Sig_Matrix, 2,
                              sum, na.rm = TRUE)
  Module_Order       <- names(sort(apply(Model_Sig_Matrix, 2,
                                         sum, na.rm = TRUE)))
  Module_Order_Index <- match(Module_Order,
                              colnames(Model_Comparison_Matrix))


  # Return the various things calculated

  Model_Module_List <- list("Model_Comparison_Matrix"  = Model_Comparison_Matrix,
                            "Fish_Backgrounds_trimmed" = Fish_Backgrounds_trimmed,
                            "Model_Sig_Matrix"         = Model_Sig_Matrix,
                            "Bin_Order_Index"          = Bin_Order_Index,
                            "Module_Order_Index"       = Module_Order_Index,
                            "Module_Order"             = Module_Order
                            )
  return(Model_Module_List)




}


#
#' Plot_Model_Module
#' @name Plot_Model_Module
#' @description Creates a bar plot of pairwise comparisons between the traits of a model genome.
#' If a distance is less than the 99% quantile of the background distances for a module of similar size,
#' The bar is presented in dark grey.
#' @seealso \code{\link{Pre_process_input}} for the full list of parameters,
#' \url{https://github.com/Jorisvansteenbrugge/TcT/wiki/The-RNAseq.data-object}
#' for more information.
#' @param Model_Module_List
#' @param Model_Bin
#' @param Module_Names
#' @export
#' @author BO Oyserman

Plot_Model_Module <- function(Model_Module_List, Model_Bin, Module_Names, margins, sortbygenome) {

num_genes <-  length( Model_Module_List[[4]])
par(mfrow = margins,
    mar   = c(2.1, 1, 2.1, 1))

# How it sorted? If by a particular genome then do the following
if (length(sortbygenome)==1) {
  j=0
  # if it is the first plot in a row, plot the labels

  for (i in order(Model_Module_List$Model_Comparison_Matrix[sortbygenome,])) {
    j = j +1
    if (j %in% (which(1:(margins[1]%*%margins[2]) %in% as.numeric(margins[2]%*%seq(1,margins[1]))==TRUE)-(margins[2]-1))) {
      barplot( Model_Module_List[[1]][ Model_Module_List[[4]], 1 ],
               col    = NA,
               border = NA,
               axes   = FALSE,
               xlim   = c(-2,1),
               ylim   = c(0,19),
               yaxt   = 'n',
               xaxt   = 'n')
      text(x = rep(-1,19),
           y = seq(from = 1, to = (num_genes-1), by = (num_genes-2) / (num_genes-1)),
           labels = rownames(Model_Module_List[[3]])[Model_Module_List[[4]]], cex = 1)

    } else {

      sig_colors                       <- rep("white", length(Model_Module_List[[4]]))
      sig_colors[Model_Module_List[[3]][,i]] <- "gray0"

      barplot(Model_Module_List[[1]][Model_Module_List[[4]], i],
              xlim = c(-2,1), horiz = TRUE,
              main = colnames(Model_Module_List[[1]])[i],
              col  = sig_colors[Model_Module_List[[4]]],
              yaxt = 'n'
      )
      abline(v = 0, lwd = 1)
    }
  }
  # How it sorted? If NOT by a particular genome then sort by most significant modules
} else {
  j=0

  for (i in rev(Model_Module_List[[5]])) {

    j = j +1
    if (j %in% (which(1:(margins[1]%*%margins[2]) %in% as.numeric(margins[2]%*%seq(1,margins[1]))==TRUE)-(margins[2]-1))) {
      barplot( Model_Module_List[[1]][ Model_Module_List[[4]], 1 ],
               col    = NA,
               border = NA,
               axes   = FALSE,
               xlim   = c(-2,1),
               ylim   = c(0,19),
               yaxt   = 'n',
               xaxt   = 'n')
      text(x = rep(-1,19),
           y = seq(from = 1, to = (num_genes-1), by = (num_genes-2) / (num_genes-1)),
           labels = rownames(Model_Module_List[[3]])[Model_Module_List[[4]]], cex = 1)

    } else {

      sig_colors                        <- rep("white", length(Model_Module_List[[4]]))
      sig_colors[Model_Module_List[[3]][,i]] <- "gray0"

      barplot(Model_Module_List[[1]][Model_Module_List[[4]], i],
              xlim = c(-2,1), horiz = TRUE,
              main = colnames(Model_Module_List[[1]])[i],
              cex.main = 0.75,
              col  = sig_colors[Model_Module_List[[4]]],
              yaxt = 'n'
      )
      abline(v = 0, lwd = 1)
    }
  }
}

}


#
#' Plot_Model_Module_Simple
#' @name Plot_Model_Module
#' @description Creates a bar plot of pairwise comparisons between the traits of a model genome.
#' If a distance is less than the 99% quantile of the background distances for a module of similar size,
#' The bar is presented in dark grey.
#' @seealso \code{\link{Pre_process_input}} for the full list of parameters,
#' \url{https://github.com/Jorisvansteenbrugge/TcT/wiki/The-RNAseq.data-object}
#' for more information.
#' @param Model_Module_List
#' @param Model_Bin
#' @param Module_Names
#' @export
#' @author BO Oyserman

Plot_Model_Module <- function(Model_Module_List, Model_Bin, Module_Names, margins, sortbygenome) {

  num_genes <-  length( Model_Module_List[[4]])
  par(mfrow = margins,
      mar   = c(2.1, 1, 2.1, 1))

  # How it sorted? If by a particular genome then do the following
  if (length(sortbygenome)==1) {

    for (i in order(Model_Module_List$Model_Comparison_Matrix[sortbygenome,])) {
        sig_colors                       <- rep("white", length(Model_Module_List[[4]]))
        sig_colors[Model_Module_List[[3]][,i]] <- "gray0"

        barplot(Model_Module_List[[1]][Model_Module_List[[4]], i],
                xlim = c(-2,1), horiz = TRUE,
                main = colnames(Model_Module_List[[1]])[i],
                col  = sig_colors[Model_Module_List[[4]]],
                yaxt = 'n'
        )
        abline(v = 0, lwd = 1)
      }
    # How it sorted? If NOT by a particular genome then sort by most significant modules
  } else {

    for (i in rev(Model_Module_List[[5]])) {
        sig_colors                        <- rep("white", length(Model_Module_List[[4]]))
        sig_colors[Model_Module_List[[3]][,i]] <- "gray0"

        barplot(Model_Module_List[[1]][Model_Module_List[[4]], i],
                xlim = c(-2,1), horiz = TRUE,
                main = colnames(Model_Module_List[[1]])[i],
                cex.main = 0.75,
                col  = sig_colors[Model_Module_List[[4]]],
                yaxt = 'n'
        )
        abline(v = 0, lwd = 1)
      }
    }
  }



##' Automatically convert a vector of strings into a color for easy plotting
##'
##' @title Convert between strings to colors
##' @param string a vector of strings representing groups.
##' @param colors a vector of colors, one for each unique element in \code{string}.
##' @export
##' @return a vector of colors, one for each element in \code{string}
##' @aliases string.to.colors stringtocolor stringToColors string.to.color
##' @author Dustin Fife
##' @seealso \code{\link{number.to.colors}}
##' @examples
##' groups = sample(LETTERS[1:5], size=100, replace=TRUE)
##' plot(rnorm(100), rnorm(100), col=string.to.colors(groups))
##' plot(rnorm(100), rnorm(100), col=string.to.colors(groups),
##'    pch=as.numeric(string.to.colors(groups, colors=c(16:20))))
##' @note This function can also be used to specify pch values, cex values, or any other plotting values
##' the user may wish to differ across groups. See examples.
string.to.colors = function(string, colors=NULL){
  if (is.factor(string)){
    string = as.character(string)
  }
  if (!is.null(colors)){
    if (length(colors)!=length(unique(string))){
      break("The number of colors must be equal to the number of unique elements.")
    }
    else {
      conv = cbind(unique(string), colors)
    }
  } else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN=function(x){conv[which(conv[,1]==x),2]}))
}


##' Automatically convert a vector of numbers into a color for easy plotting
##'
##' @title Convert from numbers to colors
##' @param value a vector of numbers.
##' @param colors a vector of two or more colors representing the inflection points of the gradients, passed to \code{\link{colorRampPalette}}.
##' @param num The number of unique intervals for colors. Chose larger numbers for finer gradients (higher resolution).
##' @export
##' @return a vector of colors.
##' @aliases number.to.color numbers.to.colors integers.to.colors integer.to.colors numberToColors numberToColor
##' @author Dustin Fife
##' @seealso \code{\link{string.to.colors}} \code{\link{colorRampPalette}} \code{\link{gradient.legend}}
##' @examples
##' #### plot three variables, represent the third with a color
##' d = mvrnorm(100, mu=c(0,0,0), Sigma=matrix(c(1, .6, .6, .6, 1, .6, .6, .6, 1), ncol=3))
##' plot(d[,1:2], col=number.to.colors(d[,3], c("black", "pink")), pch=16)
number.to.colors = function(value, colors=c("red", "blue"), num=100){
  cols = colorRampPalette(colors)(num)
  cols = 	cols[findInterval(value, vec=seq(from=min(value), to=max(value), length.out=num))]
  cols
}



##' Create a gradiented legend
##'
##' @param y the variable used to create the gradient, typically in \code{\link{number.to.colors}}
##' @param yrange The range of y values. If y is supplied, it will pulls these from the actual y values.
##' @param cols The color gradients to be used that are passed to \code{\link{colorRampPalette}}.
##' @param location The location of the subplot, expressed in fractions of the entire plot (left x, right x,
##' bottom y, top y).
##' @param n the number of values used for the gradient. Higher numbers make a higher resolution
##' @param ... other arguments passed to image
##' @aliases gradient
##' @seealso \code{\link{number.to.colors}} \code{\link{colorRampPalette}}
##' @author Dustin Fife
##' @export
##' @examples
##' y = rnorm(100); x = .6*y + rnorm(100,0,sqrt(1-.6^2))
##' randnum = runif(100)
##' plot(x,y, col=number.to.colors(randnum), pch=16)
##' gradient.legend(randnum, xlab="", ylab="")
gradient.legend = function(y=NULL, yrange=NULL, cols = c("red", "blue"),location=c(.075,.3,.575,.975), n=100,...){

  #### if they don't supply a y, make sure they give range
  if (is.null(y) & is.null(yrange)){
    stop("You must supply either a y value or a y range.")
  }
  if (is.null(yrange)){
    yrange = range(y, na.rm=T)
  }

  #### set location
  par(fig=location, new=T, las=1, ps=9)
  z=matrix(1:n,nrow=1)
  x=1
  y=seq(yrange[1],yrange[2],len=n)
  my.colors = colorRampPalette(cols)
  image(x,y,z,col=my.colors(n),axes=FALSE,...)
  axis(2)
  par(fig=c(0,1,0,1))
}
