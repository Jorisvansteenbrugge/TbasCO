.onLoad <- function(libname, pkgname){
  assign('file.path', paste0(libname,"/TcT/data/sample_data.csv"), envir = .GlobalEnv)
  assign('annotation.db.path', paste0(libname,'/TcT/data/KEGG_modules.tsv'), envir = .GlobalEnv)
}

.onDetach <- function(libname, pkgname){
  rm(file.path, envir = .GlobalEnv)
  rm(annotation.db.path, envir = .GlobalEnv)
}
# myPackageEnvironment <- new.env()
# aa <- "test"
