# # General loading ----
# normalization.features <- list('no_feature'   = c(9159700, 4459877, 9826273, 8171512, 9542765, 10522313),
#                                'ambiguous'    = c(3940698, 2023389, 4675033, 3308789, 6446272, 5966543),
#                                'library_size' = c(234232896, 183166236, 228746720, 198024002, 231567992, 259156166),
#                                'not_aligned'  = c(0, 0, 0, 0, 0, 0)
# )
#
#
#
# # Testing ----
# ko_terms <- c('K01895', 'K00925')
# fake_annotation.db <- list("module.dict" = list("test_module" = ko_terms),
#                            "all annotations in a module" = ko_terms)
#
#
# #Asume pairwise distances are good
# load(url("http://jorisvansteenbrugge.com/pairwise_distances.RData"))
# pairwise.distances <- pairwise.distances[ko_terms]
#
# attributes <- Identify_Trait_Attributes(RNAseq.data, pairwise.distances,
#                                         fake_annotation.db , threads = 2)
# members <- data.frame(attributes$test_module$clusters$membership,
#                       attributes$test_module$clusters$names)
# expected_members <- data.frame(c(2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 2, 2, 1, 2, 1, 2 ,1),
#                                c('8',  '11', '16', '17', '20', '19', '22', '25', '26', '28', '29', '31', '32', '36', '39', '42', '45', '48', '53'))
# colnames(members) <- NULL
# colnames(expected_members) <- NULL
#
# # # TEST
# expect_equal(
#   compare(members,expected_members)
# )
