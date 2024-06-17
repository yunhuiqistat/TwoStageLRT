#' Estimate Dispersion Parameters Utilizing Null Families in A Gene Cluster.
#'
#' This function uses edgeR to estimate the tagwise dispersion parameter for multiple genes (from a gene cluster) with multiple family (usually the null families detected by Pois_3mean_LRT()).
#' @import edgeR
#' @param counts  a count matrix with genes in rows, and three*n_families columns, female, male and hybrid which passed Pois_3mean_LRT() and can be treated as null families.
#' @param group a vector indicating the group information of three*n_families columns, usually one null family is one group.
#' @return a vector containing gene-wise dispersion estimation.
#' @export
#' @examples
#' library(TwoStageLRT)
#' counts <- matrix(rpois(150, 10), nrow = 10, ncol = 15)
#' group <- rep(1:5, each = 3)
#' disp_est <- Disp_est_null_family(counts = counts, group = group)
#' summary(disp_est)


Disp_est_null_family <- function(counts, group = NULL){
  if(is.null(group)){
    y <- DGEList(counts = counts)
    y <- calcNormFactors(y, method = "RLE")
    y <- estimateCommonDisp(y)
    y <- estimateTrendedDisp(y)
    y <- estimateTagwiseDisp(y)
  }
  else{
    y <- DGEList(counts = counts, group = group)
    y <- calcNormFactors(y, method = "RLE")
    design <- model.matrix(~ group)
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTrendedDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
  }
  return(y$tagwise.dispersion)
}
