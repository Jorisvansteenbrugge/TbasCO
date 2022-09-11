# Trait-BASed COmparative Transcriptomics - TbasCO

![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/jorisvsteenbrugge/tbasco?style=for-the-badge)

TbasCO is an R package for statistically comparing transcriptional patterns of trait across different organisms. It was developed using a genome-resolved time-series metatranscriptomics experiment of an enrichment microbial community simulating phosphorus removal described in the preprint: 

[TbasCO: Trait-based Comparative ‘Omics Identifies Ecosystem-Level and Niche-Differentiating Adaptations of an Engineered Microbiome.](https://www.biorxiv.org/content/10.1101/2021.12.04.471239v1) <b> McDaniel E.A.#, van Steenbrugge J.J.M.# </b>, Noguera D.R., McMahon K.D., Raaijmakers J.M., Medema M.H., Oyserman B.O. <i> bioRxiv. In Review.. </i> Dec. 2021. DOI: 10.1101/2021.12.04.471239.
<i> <br> # these authors contributed equally </i>

Installation, preprocessing of metatranscriptomes and functional annotation, and a quick-start in TbasCO are all described below. An example of the full TbasCO pipeline is also extensively describe in the `EBPR.md` vignette [vignettes/EBPR.Rmd](vignettes/EBPR.Rmd). This vignette includes using the example datasets provided by the authors in [data/sample_data.csv](data/sample_data.csv). The below tutorial and the vignettes use the KEGG database as the default trait database. However, it is possible and a simple process to add user-defined traits to the library. This is documented in [vignettes/Custom_Traits.Rmd](vignettes/Custom_Traits.Rmd).

## Installation
The package can be installed as followed from within an R environment:
```{r]
library(devtools)
install_github("jorisvansteenbrugge/TbasCO")
```

## Preprocessing Metatranscriptomes and Functional Annotation 


## TbasCO Quick-Start 
