library(magrittr)
library(reticulate)

reticulate::source_python('python/parse_module_definition.py')


modules <-
  read.csv('data/kegg_modules_2019_07_23.tsv',
           stringsAsFactors = F,
           sep = '\t') %>% .$Module %>% unique


sub_modules <- lapply(modules, py$parse_module)

modules -> names(sub_modules)
