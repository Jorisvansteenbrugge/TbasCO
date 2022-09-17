# TbasCO: Trait-BASed Comparative 'Omics

![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/jorisvsteenbrugge/tbasco?style=for-the-badge)

<b> Table of Contents </b> <br>
[Overview](https://github.com/Jorisvansteenbrugge/TbasCO#overview) <br>
[Installation](https://github.com/Jorisvansteenbrugge/TbasCO#installation) <br>
[Preprocessing Metatranscriptomes and Functional Annotations](https://github.com/Jorisvansteenbrugge/TbasCO#preprocessing-metatranscriptomes-and-functional-annotation) <br>
[TbasCO Quick-Start](https://github.com/Jorisvansteenbrugge/TbasCO#tbasco-quick-start) <br>

## Overview
TbasCO is an R package for statistically comparing transcriptional patterns of trait across different organisms. It was developed using a genome-resolved time-series metatranscriptomics experiment of an enrichment microbial community simulating phosphorus removal described in the preprint: 

[TbasCO: Trait-based Comparative ‘Omics Identifies Ecosystem-Level and Niche-Differentiating Adaptations of an Engineered Microbiome.](https://www.biorxiv.org/content/10.1101/2021.12.04.471239v1) <b> McDaniel E.A.#, van Steenbrugge J.J.M.# </b>, Noguera D.R., McMahon K.D., Raaijmakers J.M., Medema M.H., Oyserman B.O. <i> bioRxiv. </i> Dec. 2021. DOI: 10.1101/2021.12.04.471239.
<i> <br> # these authors contributed equally </i>

Installation, preprocessing of metatranscriptomes and functional annotation, and a quick-start in TbasCO are all described below. An example of the full TbasCO pipeline is also extensively describe in the `EBPR.md` vignette [vignettes/EBPR.Rmd](vignettes/EBPR.Rmd). This vignette includes using the example datasets provided by the authors in [data/sample_data.csv](data/sample_data.csv). The below tutorial and the vignettes use the KEGG database as the default trait database. However, it is possible and a simple process to add user-defined traits to the library. This is documented in [vignettes/Custom_Traits.Rmd](vignettes/Custom_Traits.Rmd).

## Installation
The package can be installed as followed from within an R environment:
```{r]
library(devtools)
install_github("jorisvansteenbrugge/TbasCO")
```

## Preprocessing Metatranscriptomes and Functional Annotation 
This tutorial for preprocessing your MAGs to use in TbasCO assumes that you have assembled your MAGs, quality-checked them, and have (optionally) confirmed the taxonomic assignment of each MAG. 

Once you have assembled and curated your representative, non-redundant MAG database from your community, you will need to: 

1. Perform functional annotations with KofamKOALA
2. Competitively map metatranscriptomic reads to your MAGs

### Dependent Tools: 

- Prokka
- KofamKOALA
- Kallisto
- R packages: tidyverse, tximport

By default, TbasCO uses the KEGG database for functional comparisons. Annotation with the KEGG database can either be done through the online [GhostKOALA](https://www.kegg.jp/ghostkoala/) portal or through the [KofamKOALA](https://www.genome.jp/tools/kofamkoala/) distribution that uses HMMER for searches against HMM profiles. KofamKOALA can either be run through the online portal or ran locally through the command line.

### Formatting Files

First, you will want to format your fasta files in a way where locus tags for each annotated coding region/protein will match between your functional predictions and the metatranscriptome mapping. One way to do this is to get all the necessary file formats (.faa, .ffn, .gbk, etc.) for each assembled genome using Prokka: 

```
for file in *.fna; do 
    N=$(basename $file .fna);
    prokka --outdir $N --prefix $N --cpus 15 $file; 
done
```

Then, concatenate the collection of .ffn and .faa files to use for competitively mapping your metatranscriptomes to and performing functional annotations, respectively: 

_Coding regions:_

```
for filename in */*.ffn; 
    do GENNAME=`basename ${filename%.ffn}`; 
    sed "s|^>|>${GENNAME}_|" $filename; 
done > all-coding-regions.ffn
```

_Proteins:_

```
for filename in */*.faa; 
    do GENNAME=`basename ${filename%.faa}`; 
    sed "s|^>|>${GENNAME}_|" $filename; 
done > all-proteins.faa
```

What these commands will do is create a concatenated file of all coding regions/proteins for all genomes, and appends the name of the genome to each locus tag. So make sure the names of your original fasta files/genome names are something informative or you know what they are (tied to the taxonomic identity for example). 

### Functional Annotation with KofamKOALA

Then annotate the concatenated protein FASTA file with KofamKOALA (download the executables and the required HMM profiles): 

```
./exec_annotation all-proteins.faa;
    -p profiles/; 
    -k ko_list; 
    -o all-proteins-kofamkoala-annots.txt; 
    --cpu 8;
```

What is nice about this pipeline and HMMs over the BLAST-based method that GhostKOALA annotates datasets with the KEGG database is that the output will give you confidence scores according to the threshold cutoff of each HMM and if the annotation for a certain protein is significant with a `*`, from which you can parse out only the significant hits and non-duplicates (i.e. taking the annotation with the highest confidence score for a particular protein). For our purposes, we want the highest confident annotation for a particular protein. You can parse this out with a series of `awk` commands: 

```
grep '*' all-proteins-kofamkoala-annots.txt | awk '{print $2"\t"$3"\t"$7" "$8" "$9" "$10" "$11" "$12" "$13" "$14}' > kofamkoala-annots-sig-mod.txt
sort -u -k1,1 kofamkoala-annots-sig-mod.txt > kofamkoala-annots-sig-mod-nodups.txt
awk '{print $1"\t"$2}' kofamkoala-annots-sig-mod-nodups.txt > ebpr-kofamkoala-annots-sig-mod-nodups-ko-list.txt
```

This will give you both a tab-delimited file of significant annotations, their scores and functional annotations, and then a list of each locus tag and KEGG identifier to use for downstream purposes, and is exactly the output that GhostKOALA gives, except without really knowing quantitatively why things were annotated the way they were. 

### Competitvely Mapping Metatranscriptome Samples to MAGs

These preprocessing steps with kallisto assume that you have performed any necessary quality filtering of metatranscriptomic reads and/or removed rRNA reads. There are always ways post-mapping to remove rRNA reads/loci, which you can do by removing genes that have abnormally high amounts of reads mapping to them or removing genes that are predicted to be ribosomal RNA by Prokka/barrnap when you performed functional annotations. 

For metatranscriptomic mapping, we used the tool `kallisto` for performing pseudoalignment of reads and quantification to coding regions. This software was originally [developed for mapping single-cell RNA-seq data](https://www.nature.com/articles/nbt.3519), but has also been demonstrated to work well with microbial community metatranscriptomics data. Using the concatenated coding regions file, you will first want to index your reference database: 

```
kallisto index -i mag-database-index all-coding-regions.ffn
```

The command for mapping for each sample is formatted as such: 

```
kallisto quant -i index -o outdir fastqfiles
```

You can also create a tab-delimited metadata file, where the first two columns are the names of the R1 and the R2 fastq files, and the 3rd column is the sample name. Such as: 

```
# SRX4072504.qced.R1.fastq	SRX4072504.qced.R2.fastq	sample1
```

Then you can pass that metadata file to a bash script that queues the mapping for each sample, so you can complete the mapping in one step for many samples. For example the bash script looks like: 

```
#! /usr/bin/bash
# Queue metatranscriptomic mapping with kallisto to reference MAG index with a metadata file specifying sample name to filenames of R1 and R2 files

index="$1"
MetadataFile="$2"

echo "Running kallisto for $MetadataFile !"

while read -r a b c; do
	echo $c $a $b;
	kallisto quant -i "$index" -o "$c" "$a" "$b";
done < "$MetadataFile"

echo "Finished running kallisto for $MetadataFile !"
```

And the command to run the script is `bash queue-kallisto.sh $index $metadataFile`. Where you input the name of your index file you created and the metadata file formatted with the fastq filenames and the sample name. This will create a directory for each sample, and in each of those directories are the kallisto mapping results that you will then parse. 

### Processing Kallisto Output 

The following steps can be performed using R. You will need an R package called `tximport` for parsing the kallisto output files. Either place the script within the directory where the output folders are (in each directory is a .tsv and .h5 file for every sample) or if you ran the jobs on a server, `scp` back all the resulting directories and lead your R script to that path. You can find the original documentation for `tximport` for parsing `kallisto` output files and other mapping programs [here](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto).

```
library(tximport)

dir <- "/path/to/your/transcriptomes/"
samples <- read.table(file.path(dir, "samples.txt"), header=TRUE)
files <- file.path(dir, samples$experiment, "abundance.h5")
# number of samples you have
names(files) <- paste0("sample", 1:6)
txi.kallisto <- tximport(files, type="kallisto", txOut = TRUE)
counts <- as.data.frame(txi.kallisto)
finalcounts <- rownames_to_column(counts, var="ID")

```

From the txi dataframe, you can then create tables or raw and normalized counts. TbasCO performs its' own normalization calculation by default to control for the abundance/differing transcript counts mapping back to each genome. The example below includes creating a `Bin` column based off of the locus tag provided to kallisto, and renaming columns based on the samples names. 

```
# raw
rawcounts <- as.data.frame(txi.kallisto$counts)
rawTable <- rownames_to_column(rawcounts, var="ID")
rawTable.split <- rawTable %>% separate(ID, c("genome"), sep='_') %>% cbind(rawTable$ID)
colnames(rawTable.split) <- c('Bin', 'B_15min_Anaerobic', 'D_52min_Anaerobic', 'F_92min_Anaerobic', 'H_11min_Aerobic', 'J_51min_Aerobic', 'N_134min_Aerobic', 'Locus_Tag')
rawOut <- rawTable.split %>% select(Bin, Locus_Tag, B_15min_Anaerobic, D_52min_Anaerobic, F_92min_Anaerobic, H_11min_Aerobic, J_51min_Aerobic, N_134min_Aerobic)
# write csv 

# normalized
normcounts <- as.data.frame(txi.kallisto$abundance)
normTable <- rownames_to_column(normcounts, var='ID')
norm.split <- normTable %>% separate(ID, c("genome"), sep='_') %>% cbind(normTable$ID)
colnames(norm.split) <- c('Bin', 'B_15min_Anaerobic', 'D_52min_Anaerobic', 'F_92min_Anaerobic', 'H_11min_Aerobic', 'J_51min_Aerobic', 'N_134min_Aerobic', 'Locus_Tag')
normOut <- norm.split %>% select(Bin, Locus_Tag, B_15min_Anaerobic, D_52min_Anaerobic, F_92min_Anaerobic, H_11min_Aerobic, J_51min_Aerobic, N_134min_Aerobic)
# write csv
```

### Merging Count Table with Functional Annotations

You can then merge your count table with your deduplicated functional annotation table with `tidyverse`. Input either the raw or normalized counts table, and the KofamKOALA annotations. The below changes the sample names to numbers instead of the specific sample names above. 

```
library(tidyverse)
# merge counts and annotations

countsFile = read.csv("raw-data/ebpr-raw-counts-names.tsv", header=FALSE, sep="\t")
annotFile = read.csv("raw-data/ebpr-kofamAnnotations.txt", header=FALSE, sep="\t")
colnames(countsFile) <- c("Bin", "Locus_Tag", "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6")
colnames(annotFile) <- c("Locus_Tag", "Annotation")
rawTable <- left_join(countsFile, annotFile)
countsTable <- rawTable[,c(2:8,9,1)]

# check to see merged correctly 
annotCounts <- countsTable %>% select(Bin, Annotation) %>% group_by(Bin) %>% mutate(Annotation.count = n()) %>% slice(1) %>% unique()

# save files
write.table(countsTable, "raw-data/full-tbasco-input-table.tsv", sep=";", row.names=FALSE, quote=FALSE)
```

By default TbasCO accepts tables with the `;` separator. This merged table now serves as your input into TbasCO.

## TbasCO Quick-Start 

#Loading and Preprocessing Data

All of the preprocessing is handled by the `Pre_process_input` function.
This function consequtively performs the following tasks:

1. Data Loading
2. Normalization
3. Filtering

### Data loading
The data set should contain the following elements:

* A column named **Gene**
* Multiple expression data columns that have at least the word **Sample** in it
* A column named **Annotation**
* A column named **Bin**

Additional columns may be present but are not used in the analysis.
Loading of data is done by assigning a file.path to the function's `filepath`
variable.

## Databases
```{r}

database <- Combine_databases(kegg_brite_20191208, kegg_module_20190723)
```

### Normalization
All of the preprocessing is handled by the `Pre_process_input` function.
This function consecutively performs the following tasks:

1. Data Loading
2. Normalization -> Default normalization method between samples and by bin
3. Filtering -> the choice between Mean Absolute Deviation ('MAD') and based on Standard Deviation ('stdev')

```{r Preprocessing}
RNAseq.data <- Pre_process_input(file.path,
                                 database = annotation.db.path,
                                 normalize.method    = T,
                                 filter.method       ='MAD',
                                 filter.low.coverage = T)
 ```                                

### Define the distance metrics to be used
In this vignette, the Pearson Correlation and the Normalized Rank Euclidean Distance are used.
The user is free to implement different metrics, as long as the function parameters
are identical to the ones in the example PC and NRED functions.
In these examples, rowA and rowB are rows from RNAseq.data$table.

```{r Defining Distance Metrics}
# Calculates the Pearson Correlation
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
```

### Main Pipeline
```{r Main Analysis}
bkgd.individual         <- Individual_Annotation_Background(RNAseq.data, 
                                                            N       = 5000, 
                                                            metrics = distance.metrics,
                                                            threads = 3)
bkgd.individual.Zscores <- Calc_Z_scores(bkgd.individual, distance.metrics)
bkgd.traits             <- Random_Trait_Background(RNAseq.data, 
                                                   bkgd.individual.Zscores,
                                                   N = 5000, 
                                                   metrics = distance.metrics,
                                                   threads = 4)
pairwise.distances      <- Calc_Pairwise_Annotation_Distance(RNAseq.data,
                                                             RNAseq.data$features$annotation.db,
                                                             distance.metrics,
                                                             bkgd.individual.Zscores,
                                                             show.progress = F,
                                                             threads = 4)
trait.attributes        <- Identify_Trait_Attributes(RNAseq.data = RNAseq.data, 
                                                     pairwise.distances = pairwise.distances,
                                                     threads = 4)
trait.attributes.pruned <- Prune_Trait_Attributes(trait.attributes, bkgd.traits, 
                                                  RNAseq.data,
                                                  p.threshold = 0.05,
                                                  pairwise.distances = pairwise.distances,
                                                  bkgd.individual.Zscores = bkgd.individual.Zscores)
sbs.trait.attributes    <- Traitattributes_To_Sbsmatrix(trait.attributes.pruned, RNAseq.data$features$bins)
```

This `sbs.trait.attributes` table is then used for downstream exploration of expression profiles and trait attributes across organisms. 
