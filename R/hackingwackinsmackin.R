prune_lalala <- function(sigboth) {

  .filter_Genome <- function(genome) {
    g_rows <- rows[which(rows$Bin == genome),]

    u.anno <- unique(as.character(g_rows$Annotation))
    c <- sum(trait.anno %in% u.anno  )

    fraction <- c/length(trait.anno)

    if (fraction >= 0.75) {
      return (T)
    } else {
      return (F)
    }

  }

  out <- c()
  for (trait in sigboth$module) {
    rows <- RNAseq.data$table[which(RNAseq.data$table$Annotation %in%
                                      RNAseq.data$features$annotation.db$module.dict[[trait]]),]
    trait.anno <- RNAseq.data$features$annotation.db$module.dict[[trait]]
    genomes <- unique(rows$Bin)

    return.vals <- sapply(genomes, .filter_Genome)

    if (sum(return.vals >= 1)) {
      out <- c(out, trait)
    }

  }

  return(out)
}


Expr_per_KO <- function () {

  sizes <- sapply(RNAseq.data$features$annotation.db$`all annotations in a module`,
                  function(KO) {
                    rows <- RNAseq.data$table[which(RNAseq.data$table$Annotation == KO),]
                    if (nrow(rows) >= 1){
                      return(sum(rows[,RNAseq.data$features$rank.columns]))
                    } else {
                      return(1000000)
                    }

                  }
           )

}


