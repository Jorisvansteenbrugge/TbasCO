library(magrittr)
library(dplyr)
library(ggplot2)
library(plotly)

setwd("/home/joris/nemaNAS/steen176/TbasCO/benchmark_cutoff/")
benchmarking_files <- list.files("/home/joris/nemaNAS/steen176/TbasCO/benchmark_cutoff/", pattern = '.RData', full.names = F)

load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  return(env[[nm]])
}

# Calculate the full number of trait attributes that are created
get_tas <- function(TAs){
  count <- 0
  for(trait in TAs){
    count <- count + trait$clusters$membership %>% unique %>% length

  }
  return(count)
}


get_sig_tas <- function(SigTAs){
  return(sapply(SigTAs, length) %>% as.numeric %>% sum)
}

get_num_genomes <- function(RNAseq.data){
  return(RNAseq.data$features$bins %>%  length)
}

plot_data <- matrix(ncol = 4, nrow = 0)

for(benchmark in benchmarking_files){
  data <- load_obj(benchmark)

  ta_count <- get_tas(data$TAs)
  sig_ta_count <- get_sig_tas(data$SigTAs)
  num_genomes <- get_num_genomes(data$RNAseq.data)

  plot_data <- rbind(plot_data, c(benchmark, ta_count, sig_ta_count, num_genomes))

}

plot_data %<>% as.data.frame(stringsAsFactors = F) %>% transmute(Label = V1,
                                                                 TAs = as.numeric(V2),
                                                                 SigTas = as.numeric(V3),
                                                                 NumGenomes =  as.numeric(V4))
plot_data$Label %<>% factor(., levels = sort(.))

library(gridExtra)
a <- ggplot(plot_data, aes(x = TAs, y = SigTas)) + geom_point()
d <- ggplot(plot_data, aes(x = NumGenomes, y = SigTas)) + geom_point()
b <- ggplot(plot_data, aes(y = SigTas, x = Label)) + geom_bar( stat = 'identity')
c <- ggplot(plot_data, aes(y = NumGenomes, x = Label)) + geom_bar(stat = 'identity')
grid <- gridExtra::grid.arrange(a,d, ncol = 2)

grid.arrange(b,c, ncol = 1)

