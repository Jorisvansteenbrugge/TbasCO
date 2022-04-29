# Trait-BASed COmparative Transcriptomics - TbasCO


![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/jorisvsteenbrugge/tbasco?style=for-the-badge)


## Installation
The package can be installed as followed from within an R environment:
```{r]
library(devtools)
install_github("jorisvansteenbrugge/TbasCO")
```

An example of the full pipeline is extensively desiribed in the EBPR.rmd vignette [vignettes/EBPR.Rmd](vignettes/EBPR.Rmd). This vignette includes a way to use an example dataset provided by the authors provided in [data/sample_data.csv](data/sample_data.csv).

In the EBPR.Rmd vignette, the default trait database based on KEGG is used. It is however a simple process to add user defined traits to the library. This is documented in [vignettes/Custom_Traits.Rmd](vignettes/Custom_Traits.Rmd).
