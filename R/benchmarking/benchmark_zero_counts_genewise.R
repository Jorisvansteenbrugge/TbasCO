library(TbasCO)
library(magrittr)

options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)


database <- Combine_databases(kegg_brite_20191208, kegg_module_20190723)

RNAseq.data <- Pre_process_input(args[1],
                                 database = database,
                                 normalize.method    = F,
                                 filter.method       = 'stdev',
                                 filter.low.coverage = F,
                                 filter.genes.zero = args[2])
pre_ret = list("RNAseq.data" = RNAseq.data)
save(pre_ret, file = args[4])

PC <- function(rowA, rowB, RNAseq.features){
  return(cor(as.numeric(rowA[RNAseq.features$sample.columns]),
             as.numeric(rowB[RNAseq.features$sample.columns])
  )
  )
}

# Calculates the Normalized Rank Euclidean Distance
NRED <- function(rowA, rowB, RNAseq.features) {
  r.A <- as.numeric(rowA[ RNAseq.features$rank.columns ])
  r.B <- as.numeric(rowB[ RNAseq.features$rank.columns ])
  return(
    sum((r.A - r.B) * (r.A - r.B))
  )
}

# Combine multiple distance metrics to complement each other.
distance.metrics <- list("NRED" = NRED,
                         "PC"   = PC)

bkgd.individual         <- Individual_Annotation_Background(RNAseq.data,
                                                            N       = 5000,
                                                            metrics = distance.metrics,
                                                            threads = 5)

bkgd.individual.Zscores <- Calc_Z_scores(bkgd.individual, distance.metrics)

bkgd.traits             <- Random_Trait_Background(RNAseq.data,
                                                   bkgd.individual.Zscores,
                                                   N = 5000,
                                                   metrics = distance.metrics,
                                                   threads = 5)

pairwise.distances      <- Calc_Pairwise_Annotation_Distance(RNAseq.data,
                                                             RNAseq.data$features$annotation.db,
                                                             distance.metrics,
                                                             bkgd.individual.Zscores,
                                                             show.progress = F,
                                                             threads = 5)

trait.attributes        <- Identify_Trait_Attributes(RNAseq.data = RNAseq.data,
                                                     pairwise.distances = pairwise.distances,
                                                     threads = 5)

trait.attributes.pruned <- Prune_Trait_Attributes(trait.attributes, bkgd.traits,
                                                  RNAseq.data,
                                                  p.threshold = 0.05,
                                                  pairwise.distances = pairwise.distances,
                                                  bkgd.individual.Zscores = bkgd.individual.Zscores)

sbs.trait.attributes    <- Traitattributes_To_Sbsmatrix(trait.attributes.pruned,
                                                        RNAseq.data$features$bins)
ret <- list("RNAseq.data" = RNAseq.data,
            "trait.attributes" = trait.attributes,
            "trait.attributes.pruned" = trait.attributes.pruned)
save(ret, file = args[3])
