library(magrittr)
library(dplyr)
library(ggplot2)


setwd("/home/joris/scratch/tbasco/")
benchmarking_files <- list.files("/home/joris/scratch/tbasco/", pattern = '.RData', full.names = F)

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

get_num_genes <- function(RNAseq.data){
    return (RNAseq.data$table %>% nrow)
}

get_sig_tas <- function(SigTAs){
    return(sapply(SigTAs, length) %>% as.numeric %>% sum)
}

get_num_genomes <- function(RNAseq.data){
    return(RNAseq.data$features$bins %>%  length)
}


get_genomes_per_attribute <- function(SigTas){
    numbers <- c()
    
    for(trait in SigTas){
        for(attri in trait){
            g.len <- attri$genomes %>% length
            numbers <- c(numbers,g.len)
        }
    }
    return(numbers)
}






plot_data <- matrix(ncol = 4, nrow = 0)

avg_num_genomes <- data.frame(matrix(ncol=2, nrow = 0), stringsAsFactors = FALSE)

for(benchmark in benchmarking_files){
    data <- load_obj(benchmark)
    
    genes.count <- get_num_genes(data$RNAseq.data)
    if (is.null(data$trait.attributes.pruned)) {
        sig_ta_count <- 0
        ta_count <- 0
        
    } else {
        sig_ta_count <- get_sig_tas(data$trait.attributes.pruned)
        ta_count <- get_tas(data$trait.attributes)

        counts <- get_genomes_per_attribute(data$trait.attributes.pruned)
        for (c in counts){
        avg_num_genomes <- rbind(avg_num_genomes, c(benchmark, c), stringsAsFactors = FALSE)
        }
        
    }
    
    num_genomes <- get_num_genomes(data$RNAseq.data)
    
  
    
    
    plot_data <- rbind(plot_data, c(benchmark, num_genomes, genes.count, sig_ta_count ))

}


    plot_data %<>% as.data.frame(stringsAsFactors = F) %>% transmute(Label = V1,
                                                                     NumGenomes = as.numeric(V2),
                                                                     NumGenes = as.numeric(V3),
                                                                     SigTas = as.numeric(V4))

colnames(avg_num_genomes) <- c("Benchmark", "Count")
avg_num_genomes$Count %<>% as.numeric
#pdf(file = '~/scratch/tbasco/sigTasPrelim.pdf')
ggplot(data = avg_num_genomes) + geom_violin(aes(x = Benchmark, y = Count))
#dev.off()