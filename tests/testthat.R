library(magrittr)
library(testthat)
library(TcT)

RNAseq.data <- Pre_process_input(url("http://jorisvansteenbrugge.com/sample_data.csv"),
                                 url("http://jorisvansteenbrugge.com/KEGG_modules.tsv"),
                                 normalize.method    = T,
                                 filter.method       = 'MAD',
                                 filter.low.coverage = T,
                                 normalization.features)

test_check("TcT")
