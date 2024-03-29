% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{Pre_process_input}
\alias{Pre_process_input}
\title{Pre-processing of RNAseq data}
\usage{
Pre_process_input(
  file.path,
  database,
  normalize.method = T,
  filter.method = "MAD",
  filter.low.coverage = T,
  filter.genes.zero = F,
  normalization.features = NULL,
  taxon_file = NULL
)
}
\arguments{
\item{normalize.method}{User defined method to normalize data}

\item{filter.low.coverage}{boolean expression to wether we should filter out}

\item{filter.genes.zero}{False or a number. The number indicates the number of zeros that can max be present. Genes with more zeros than the cut-off get pruned.
genomes with a low overall expression level.}

\item{filepath}{RNAseq data file}
}
\value{
RNAseq.data A list object containin the RNAseq table with ranked expression collumns possible filtered and normalized,
together with some features (position of expression collumns, rank collumns, annotation database, a list of the genomes, and
a presence absence table of all annotations)
}
\description{
Wrapper of preprocess functions, including: Identifying Matrix
features, Normalization, Filtering for non informative lines
}
\examples{
Pre_process_input(file.path, normalize.method = FALSE, filter.method = FALSE) #No normalization or filtering
Pre_process_input(file.path, normalize.method = TRUE, filter.method = 'MAD') # Default Normalization, and MAD filtering
Pre_process_input(file.path, normalize.method = custom_function_name, filter.method = FALSE) # Custom Normalization, no filtering
}
\author{
JJM van Steenbrugge
}
