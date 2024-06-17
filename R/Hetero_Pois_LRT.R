#' Poisson Likelihood Ratio Test for Detection of Mid-parent Heterosis
#'
#' This function conducts Poisson LRT for detection of gene expression mid-parent heterosis (MPH) in a single family with two parental genotypes and one hybrid genotype and for a single gene with one sample per genotype.
#' Denote the mean for gene g and female parent as mu_g1, male parent as mu_g2 and hybrid as mu_g3, the null hypothesis is Hg0: mu_g1 + mu_g2 = 2 * mu_g3.
#' @param y_g1 raw count for female.
#' @param y_g2 raw count for male.
#' @param y_g3 raw count for hybrid.
#' @param C_1 normalization factor for female sample.
#' @param C_2 normalization factor for male sample.
#' @param C_3 normalization factor for hybrid sample.
#' @export
#' @return p-value.
#' @examples
#' library(TwoStageLRT)
#' Hetero_Pois_LRT(y_g1 = 10, y_g2 = 15, y_g3 = 20, C_1 = 1, C_2 = 1.2, C_3 = 0.8)



Hetero_Pois_LRT <- function(y_g1, y_g2, y_g3, C_1 = 1, C_2 = 2, C_3 = 1){
  # 1. Get MLE under Hg0
  if(C_1 > C_2){ # we solve mu_g2 first and then mu_g1
    a <- C_1 - C_2
    b <- 2*C_2 + C_3
    A <- a*b
    B <- b*(y_g1 + y_g2) - 2*a*(y_g2 + y_g3)
    C <- -2*y_g1*y_g2 - 2*y_g2^2 - 2*y_g2*y_g3
    # take the positive root
    mu_g2_hat0 <- ifelse((-B+sqrt(B^2 - 4*A*C)) / (2*A) > 0,
                         (-B+sqrt(B^2 - 4*A*C)) / (2*A) ,(-B-sqrt(B^2 - 4*A*C)) / (2*A))
    mu_g1_hat0 <- y_g1*mu_g2_hat0 / (y_g2 + a*mu_g2_hat0)
  }
  else if(C_1 < C_2){ # we solve mu_g1 first and then mu_g2
    a <- C_1 - C_2
    b <- 2*C_1 + C_3
    A <- a*b
    B <- -b*(y_g1 + y_g2) - 2*a*(y_g1 + y_g3)
    C <- 2*y_g1*y_g2 + 2*y_g1^2 + 2*y_g1*y_g3
    # take the positive root
    mu_g1_hat0 <- ifelse((-B-sqrt(B^2 - 4*A*C)) / (2*A) > 0,
                         (-B-sqrt(B^2 - 4*A*C)) / (2*A), (-B+sqrt(B^2 - 4*A*C)) / (2*A))
    mu_g2_hat0 <- y_g2*mu_g1_hat0 / (y_g1 - a*mu_g1_hat0)
  }
  else{ # when C_1 = C_2
    # if y_g1 + y_g2 = 0, this means y_g1=y_g2=0, means numerators are all 0, thus y_g1+y_g2 in denominator replaced 1 does not affect the MLE of mu_g12
    mu_g1_hat0 <- 2*(y_g1^2 + y_g1*y_g2 + y_g1*y_g3) / ((2*C_1 + C_3)*ifelse(y_g1 + y_g2 == 0, 1, y_g1 + y_g2))
    mu_g2_hat0 <- 2*(y_g2^2 + y_g1*y_g2 + y_g2*y_g3) / ((2*C_2 + C_3)*ifelse(y_g1 + y_g2 == 0, 1, y_g1 + y_g2))
  }

  mu_g3_hat0 <- (mu_g1_hat0 + mu_g2_hat0) / 2
  mu_g_hat0 <- c(mu_g1_hat0, mu_g2_hat0, mu_g3_hat0)

  # 2. Get MLE under Hg1
  mu_g1_hat1 <- y_g1/C_1
  mu_g2_hat1 <- y_g2/C_2
  mu_g3_hat1 <- y_g3/C_3
  mu_g_hat1 <- c(mu_g1_hat1, mu_g2_hat1, mu_g3_hat1)

  # 3. Get test statistic
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
  y_set = c(y_g1, y_g2, y_g3)
  C_set = c(C_1, C_2, C_3)
  LR <- -2 * (logL(y_set = y_set, C_set = C_set, mu_set = mu_g_hat0) - logL(y_set = y_set, C_set = C_set, mu_set = mu_g_hat1))
  pval <- pchisq(LR, df = 1, lower.tail = FALSE)
  return(pval)
}
