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
#' @example Network_Trait_Genomes(c("M00002", "M00007"), trait.attributes.pruned)
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

  all_scores <- c(
    bkgd.individual.Zscores$zscores$`Random Annotated Genes`$PC,
    bkgd.individual.Zscores$zscores$`Random Annotated Genes`$NRED,
    bkgd.individual.Zscores$zscores$`Genes with the same annotation`$PC,
    bkgd.individual.Zscores$zscores$`Genes with the same annotation`$NRED
  )



  random.annotated.genes.hexb <- hexbin(bkgd.individual.Zscores$zscores$`Random Annotated Genes`$PC,
    bkgd.individual.Zscores$zscores$`Random Annotated Genes`$NRED,
    ybnds = c(min(all_scores), max(all_scores))
  )

  random.identical.annotated.genes.hexb <- hexbin(bkgd.individual.Zscores$zscores$`Genes with the same annotation`$PC,
    bkgd.individual.Zscores$zscores$`Genes with the same annotation`$NRED,
    ybnds = c(min(all_scores), max(all_scores))
  )

  cnt.max <- max(c(
    random.identical.annotated.genes.hexb@count,
    random.annotated.genes.hexb@count
  ))


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
  cat.modules <- list()
  con <- file("/home/joris/kegg_pathway_module.csv", "r")
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
      break
    }
    line <- unlist(strsplit(line, ";"))

    cat.modules[[line[1]]] <- unique(line[2:length(line)])
  }

  close(con)

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
      ylab = "Number of attributes shared with CAA"
    )
  axis(1, at = 1:length(bin_occurences), labels = names(bin_occurences))
}

Plot_Redundancy_Traits <- function(RNAseq.data) {
  library(ggplot2)
  library(magrittr)

  ta.pa <- apply(RNAseq.data$features$trait_presence_absence, 1, function(row) {
    return(sum(row))
  })

  ta.matrix <- data.frame(cbind(as.numeric(ta.pa), names(ta.pa)))
  colnames(ta.matrix) <- c("count", "trait")

  ta.matrix$trait <- factor(ta.matrix$trait, levels = ta.matrix$trait[order(-ta.matrix$count)])



  barplot.matrix <- matrix(ncol = 3, nrow = 0)
  colnames(barplot.matrix) <- c("Trait", "Attribute", "Count")

  # Highest number of attributes in 1 trait
  max.ta <- sapply(trait.attributes.pruned, length) %>% max
  colour.scheme <- grey.colors( n = ( max.ta) )

  for(trait.name in names(trait.attributes.pruned)) {

    trait <- trait.attributes.pruned[[trait.name]]
    attribute.n <- trait %>% length


    if (attribute.n >= 1 ) {

      for(ta.idx in 1:attribute.n){
        ta <- trait[[ta.idx]]

        count <- length(ta$genomes) %>% as.numeric
        ta.name <- sprintf("Attribute%s", as.character(ta.idx))
        barplot.matrix <- rbind(barplot.matrix,
                                c(trait.name, ta.name, count))
      }

      max.genomes     <- ta.pa[trait.name]
      current.genomes <- barplot.matrix[which(barplot.matrix[,"Trait"] == trait.name), 'Count'] %>%
        as.numeric %>% sum

    }
 }

  barplot.df <- as.data.frame(barplot.matrix)

  barplot.df$Count <- as.numeric(barplot.df$Count)


  barplot.df$Trait <- factor(barplot.df$Trait, levels = levels(ta.matrix$trait) )

  ggplot(data = barplot.df, aes(x = Trait, y = Count, fill = Attribute, colour = 'black')) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values = colour.scheme)


  #
  # ggplot(data = ta.matrix, aes(x = trait, y = count)) +
  #   geom_bar(stat = "identity")
}

Plot_traits_vs_attributes <- function() {
  point.matrix <- matrix(ncol=3, nrow=0)
  t.pa <- RNAseq.data$features$trait_presence_absence

  for (pair in combn(RNAseq.data$features$bins, 2, simplify = F)){
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

  plot(x = as.numeric(point.matrix[,2]),
       y = as.numeric(point.matrix[,3]),
       xlab = '# Overlap Traits',
       ylab = '# Overlap Attributes',
       main = "Pairwise Genome Comparisons")


  model <- lm(attributes~traits, data = point.df)
  abline(model$coefficients)

  cinterval <- confint(model, level=.99)
  abline(cinterval[,1])
  abline(cinterval[,2])
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
