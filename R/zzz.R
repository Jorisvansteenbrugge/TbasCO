.onLoad <- function(libname, pkgname){
  os <- .Platform$OS.type
  if(os == "unix"){
    assign('file.path', paste0(libname,"/TcT/data/sample_data.csv"), envir = .GlobalEnv)
    assign('annotation.db.path', paste0(libname,'/TcT/data/KEGG_modules.tsv'), envir = .GlobalEnv)
    assign('ko.db.path', paste0(libname,'/TcT/data/KO_identifiers.keg'), envir = .GlobalEnv)
    load(paste0(libname,'/TcT/data/sub_modules.RData'), envir = .GlobalEnv)
  }else if(os == 'windows'){
    assign('file.path', paste0(libname,"\\TcT\\data\\sample_data.csv"), envir = .GlobalEnv)
    assign('annotation.db.path', paste0(libname,'\\TcT\\data\\KEGG_modules.tsv'), envir = .GlobalEnv)
    assign('ko.db.path', paste0(libname,'\\TcT\\data\\KO_identifiers.keg'), envir = .GlobalEnv)
    load(paste0(libname,'\\TcT\\data\\sub_modules.RData'), envir = .GlobalEnv)
  }

}

.onDetach <- function(libname, pkgname){
  rm(file.path, envir = .GlobalEnv)
  rm(annotation.db.path, envir = .GlobalEnv)
}
# myPackageEnvironment <- new.env()
# aa <- "test"
