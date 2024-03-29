% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utility.R
\name{Model_Module}
\alias{Model_Module}
\title{Model Module}
\usage{
Model_Module(
  RNAseq.data,
  trait.attributes,
  Model_Bin,
  Module_Names,
  bkgd.traits
)
}
\arguments{
\item{RNAseq.data}{Collection of multple components, include RNA seq data,
annotations, etc.}

\item{trait}{character string with the name of a trait.}

\item{Yrange}{It is optional to provide a Yrange}
}
\description{
Creates a bar plot of pairwise comparisons between the traits of a model genome.
If a distance is less than the 99% quantile of the background distances for a module of similar size,
The bar is presented in dark grey.
}
\examples{
Module_Model_List <- Model_Module(RNAseq.data, trait.attributes, '39', Module_Names, c(-2.5,2.5), bkgd.traits)
\dontrun{
Model_Bin <- 39
Bin_Order <- c(48,32,31,29,22,11,39,16,53,45,42,28,20,25,19,8,36,26,17)
Yrange <- c(-2.5,2.5)
Module_Names_ATP_synthesis <- c("M00144","M00143","M00146","M00147","M00149","M00150","M00148","M00162","M00151","M00152","M00154","M00155","M00153","M00417","M00416","M00156","M00157","M00159")
}
}
\seealso{
\code{\link{Pre_process_input}} for the full list of parameters,
\url{https://github.com/Jorisvansteenbrugge/TcT/wiki/The-RNAseq.data-object}
for more information.
}
\author{
BO Oyserman
}
