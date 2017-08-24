
#' Pre-processing of RNAseq data
#' @export
Pre_Process <- function(RNAseq.table){
  RNAseq.table <- Normalization(RNAseq.table)
}



Normalization <- function(RNAseq.table, normalize_function){

  SS <- matrix_features@SS
  SE <- matrix_features@SE
  no_feature <- matrix_features@no_feature
  ambiguous <- matrix_features@ambiguous
  not_aligned <- matrix_features@not_aligned

  if(method == "default"){
    return(defaultRNA_Normalize(SS, SE, RNAseq_Annotated_Matrix, no_feature, ambiguous,not_aligned))
  }
  else if(method == "TMM" || method == "RLE"){
    return(edgeRmethods(SS, SE, method, RNAseq_Annotated_Matrix))
  }
}


edgeRmethods <- function(SS, SE, method_name, RNAseq_Annotated_Matrix){
  library(edgeR)
  RNAseq_Annotated_Matrix1<-DGEList(counts=RNAseq_Annotated_Matrix[, SS:SE],lib.size=library_size)
  norm_factors <- calcNormFactors(RNAseq_Annotated_Matrix1, method=method_name)
  RNAseq_Annotated_Matrix[, SS:SE] <- t(t(as.matrix(norm_factors)) / norm_factors$samples[,3])
  return(RNAseq_Annotated_Matrix)
}
