#' Negative Binomial Likelihood Ratio Test with Known Dispersion for Detection of Mid-parent Heterosis
#'
#' This function conducts negative binomial LRT with known dispersion for detection of gene expression mid-parent heterosis (MPH) in a single family with two parental genotypes and one hybrid genotype and for a single gene with one sample per genotype.
#' Denote the mean for gene g and female parent as mu_g1, male parent as mu_g2 and hybrid as mu_g3, the null hypothesis is Hg0: mu_g1 + mu_g2 = 2 * mu_g3.
#' Refer the paper for details of this test, since there is no closed form for the maximum likelihood estimate under null hypothesis, numerical algorithms are used with the global MLE as the starting points. Algorithms will fail if any of the raw counts is 0.
#' @param y_g1 raw count for female.
#' @param y_g2 raw count for male.
#' @param y_g3 raw count for hybrid.
#' @param C_1 normalization factor for female sample.
#' @param C_2 normalization factor for male sample.
#' @param C_3 normalization factor for hybrid sample.
#' @param dispersion value of dispersion parameter.
#' @import nleqslv
#' @import rootSolve
#' @export
#' @return p-value.
#' @examples
#' library(TwoStageLRT)
#' Hetero_NB_LRT(y_g1 = 10, y_g2 = 15, y_g3 = 20, C_1 = 1, C_2 = 1.2, C_3 = 0.8, dispersion = 0.1)

Hetero_NB_LRT <- function(y_g1, y_g2, y_g3, C_1 = 1, C_2 = 2, C_3 = 1, dispersion = 1){


  # rate = 1/dispersion
  if(dispersion == 0){return(NA)}
  r <- 1/dispersion

  # 1. Get MLE under global
  mu_g_hat1 <- c(y_g1/C_1, y_g2/C_2, y_g3/C_3)

  # 2. Get MLE under null
  y_set <- c(y_g1, y_g2, y_g3)
  C_set <- c(C_1, C_2, C_3)
  if(sum(y_set) == 0){return(1)}
  else{
    equations <- function(mu_set) {
      eq1 <- y_set[1]/mu_set[1] + y_set[3]/(mu_set[1]+mu_set[2]) - C_1*(r + y_set[1])/(C_1*mu_set[1] + r) - C_3*(r + y_set[3])/(C_3*(mu_set[1] + mu_set[2]) + 2*r)
      eq2 <- y_set[2]/mu_set[2] + y_set[3]/(mu_set[1]+mu_set[2]) - C_2*(r + y_set[2])/(C_2*mu_set[2] + r) - C_3*(r + y_set[3])/(C_3*(mu_set[1] + mu_set[2]) + 2*r)
      return(c(eq1, eq2))
    }

    init_list <- list(c(ifelse(y_g1 == 0, 10, y_g1), ifelse(y_g2 == 0, 10, y_g2)),
                      c(y_g1/C_1, y_g2/C_2), c(y_g1/2, y_g2/2),
                      c(10, 10), c(100, 100), c(2, 2), c(100, 10), c(10, 100),
                      c(100, 2), c(2, 100))

    # function to check the validity of the solution
    check_solution <- function(sol) {
      all(sol$x >= 0) && all(sol$x <= max(c(y_g1, y_g2, y_g3)) * 100) && sol$termcd < 4
    }

    # loop through the starting points to find a valid solution
    sol_found <- FALSE
    for (inits in init_list) {
      solution <- nleqslv(x = inits, fn = equations, control = list(maxit = 500))
      if (check_solution(solution)) {
        mu_g_hat0 <- c(solution$x, mean(solution$x))
        sol_found <- TRUE
        break
      }
    }

    # if no solution is found with nleqslv, use multiroot
    if (!sol_found) {
      root <- multiroot(f = equations, start = c(ifelse(y_g1 == 0, 10, y_g1/C_1), ifelse(y_g2 == 0, 10, y_g2/C_2)))
      if (all(root$root >= 0) && all(root$root <= max(c(y_g1, y_g2, y_g3)) * 100)) {
        mu_g_hat0 <- c(root$root, mean(root$root))
      }
      else {
        return(NA) # Or handle the failure to find a solution as appropriate
      }
    }

    logL <- function(y_set, C_set, mu_set, r){
      # under H1, if y_gi!=0, mu_gi!=0, under H0, its very rare to get 0 as MLE
      if(mu_set[1] == 0){first <- -(r + y_set[1])*log(r + C_set[1]*mu_set[1])}
      else{first <- y_set[1]*log(mu_set[1]) - (r + y_set[1])*log(r + C_set[1]*mu_set[1])}
      if(mu_set[2] == 0){second <- -(r + y_set[2])*log(r + C_set[2]*mu_set[2])}
      else{second <- y_set[2]*log(mu_set[2]) - (r + y_set[2])*log(r + C_set[2]*mu_set[2])}
      if(mu_set[3] == 0){third <- -(r + y_set[3])*log(r + C_set[3]*mu_set[3])}
      else{third <- y_set[3]*log(mu_set[3]) - (r + y_set[3])*log(r + C_set[3]*mu_set[3])}
      return(first + second + third)
      # Since the choose terms and rlogr and ylogC terms will cancel even with y=0, but it can be greater than the computation can afford  to make computation fail, so we omit them right now
    }
    LR <- -2 * (logL(y_set=y_set, C_set = C_set, mu_set = mu_g_hat0, r=r) - logL(y_set=y_set, C_set = C_set, mu_set = mu_g_hat1, r=r))
    pval <- pchisq(LR, df = 1, lower.tail = FALSE)
    return(pval)
  }
}
