#' Poisson Likelihood Ratio Test for Equal Family Mean
#' This function conducts Poisson LRT test for detection of null-family in a single family and for a single gene with single replication per genotype.
#' Denote the mean for gene g and female parent as mu_g1, male parent as mu_g2 and hybrid as mu_g3, the null hypothesis is Hg0: mu_g1 = mu_g2 = mu_g3.
#' @param y_g1 raw count for female.
#' @param y_g2 raw count for male.
#' @param y_g3 raw count for hybrid.
#' @param C_1 normalization factor for female sample.
#' @param C_2 normalization factor for male sample.
#' @param C_3 normalization factor for hybrid sample.
#' @return p-value.
#' @export
#' @examples
#' library(TwoStageLRT)
#' Pois_3mean_LRT(y_g1 = 10, y_g2 = 15, y_g3 = 20, C_1 = 1, C_2 = 1.2, C_3 = 0.8)



Pois_3mean_LRT <- function(y_g1, y_g2, y_g3, C_1 = 1, C_2 = 1, C_3 = 1){

  y_set <- c(y_g1, y_g2, y_g3)
  C_set <- c(C_1, C_2, C_3)
  mu_g_hat0 <- c(sum(y_set)/sum(C_set), sum(y_set)/sum(C_set), sum(y_set)/sum(C_set))
  mu_g_hat1 <- c(y_g1/C_1, y_g2/C_2, y_g3/C_3)
  logL <- function(y_set, C_set, mu_set){
    # if y_gi != 0, mu_gi != 0, for i=1,2,3
    if(mu_set[1] == 0){first <- 0}
    else{first <- y_set[1]*log(mu_set[1]) - C_set[1]*mu_set[1]}
    if(mu_set[2] == 0){second <- 0 }
    else{second <- y_set[2]*log(mu_set[2]) - C_set[2]*mu_set[2]}
    if(mu_set[3] == 0){third <- 0}
    else{third <- y_set[3]*log(mu_set[3]) - C_set[3]*mu_set[3]}
    return(first + second + third)
    #- log(factorial(y_set[1])) - log(factorial(y_set[2])) - log(factorial(y_set[3]))
    #+ y_set[1]log(C_set[1]) + y_set[2]log(C_set[2]) + y_set[3]log(C_set[3])
    # Since the factorial and ylogC terms will cancel even when y=0, but it can be infinity to make computation fail, so we omit them.
  }
  LR <- -2 * (logL(y_set = y_set, C_set = C_set, mu_set = mu_g_hat0) - logL(y_set = y_set, C_set = C_set, mu_set = mu_g_hat1))
  pval <- pchisq(LR, df = 2, lower.tail = FALSE)
  return(pval)
}
