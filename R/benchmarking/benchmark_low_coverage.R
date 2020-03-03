library(TbasCO)
library(magrittr)

options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

database <- Combine_databases(kegg_brite_20191208, kegg_module_20190723)
normalization.features <- list('no_feature'   = c(9159700, 4459877, 9826273, 8171512, 9542765, 10522313),
                               'ambiguous'    = c(3940698, 2023389, 4675033, 3308789, 6446272, 5966543),
                               'library_size' = c(234232896, 183166236, 228746720, 198024002, 231567992, 259156166),
                               'not_aligned'  = c(0, 0, 0, 0, 0, 0)
)

RNAseq.data <- Pre_process_input(file.path,
                                database = database,
                                normalize.method    = T,
                                normalization.features = normalization.features,
                                filter.method       ='MAD',
                                filter.low.coverage = c(args[1], args[2]))

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
                                                            N       = 10000,
                                                            metrics = distance.metrics,
                                                            threads = 1)

bkgd.individual.Zscores <- Calc_Z_scores(bkgd.individual, distance.metrics)

bkgd.traits             <- Random_Trait_Background(RNAseq.data,
                                                   bkgd.individual.Zscores,
                                                   N = 10000,
                                                   metrics = distance.metrics,
                                                   threads = 1)

pairwise.distances      <- Calc_Pairwise_Annotation_Distance(RNAseq.data,
                                                             RNAseq.data$features$annotation.db,
                                                             distance.metrics,
                                                             bkgd.individual.Zscores,
                                                             show.progress = F,
                                                             threads = 1)

trait.attributes        <- Identify_Trait_Attributes(RNAseq.data = RNAseq.data,
                                                     pairwise.distances = pairwise.distances,
                                                     threads = 1)

trait.attributes.pruned <- Prune_Trait_Attributes(trait.attributes, bkgd.traits,
                                                  RNAseq.data,
                                                  p.threshold = 0.05,
                                                  pairwise.distances = pairwise.distances,
                                                  bkgd.individual.Zscores = bkgd.individual.Zscores)

save(list("RNAseq.data" = RNAseq.data,
          'TAs' = trait.attributes,
          "SigTAs" = trait.attributes.pruned,
          'bkgd.traits' = bkgd.traits,
          'zscores' = bkgd.individual.Zscores),
          file = args[3])
