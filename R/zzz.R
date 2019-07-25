.onLoad <- function(libname, pkgname){
  library(magrittr)
  library(reticulate)
  os <- .Platform$OS.type
  if(os == "unix"){
    assign('file.path', paste0(libname,"/TbasCO/data/sample_data.csv"), envir = .GlobalEnv)
    assign('kegg_categories', paste0(libname,"/TbasCO/data/kegg_categories.keg"), envir = .GlobalEnv)
    assign('kegg_categories_script', paste0(libname,"/TbasCO/data/get_module_categories.py"), envir = .GlobalEnv)

    assign('annotation.db.path', paste0(libname,'/TbasCO/data/kegg_modules_2019_07_23.tsv'), envir = .GlobalEnv)
    assign('ko.db.path', paste0(libname,'/TbasCO/data/KO_identifiers.keg'), envir = .GlobalEnv)

    load(paste0(libname,'/TbasCO/data/sub_modules.RData'), envir = .GlobalEnv)
  }else if(os == 'windows'){
    assign('file.path', paste0(libname,"\\TbasCO\\data\\sample_data.csv"), envir = .GlobalEnv)
    assign('kegg_categories', paste0(libname,"\\TbasCO\\data\\kegg_categories.keg"), envir = .GlobalEnv)
    assign('kegg_categories_script', paste0(libname,"\\TbasCO\\data\\get_module_categories.py"), envir = .GlobalEnv)
    assign('annotation.db.path', paste0(libname,'\\TbasCO\\data\\kegg_modules_2019_07_23.tsv'), envir = .GlobalEnv)
    assign('ko.db.path', paste0(libname,'\\TbasCO\\data\\KO_identifiers.keg'), envir = .GlobalEnv)
    load(paste0(libname,'\\TbasCO\\data\\sub_modules.RData'), envir = .GlobalEnv)
  }



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

  assign("distance.metrics", distance.metrics, envir = .GlobalEnv)
}

.onDetach <- function(libname, pkgname){
  rm(file.path, envir = .GlobalEnv)
  rm(annotation.db.path, envir = .GlobalEnv)
}
# myPackageEnvironment <- new.env()
# aa <- "test"
