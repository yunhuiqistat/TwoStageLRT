#' Simulation dispersion parameters as a function of mid-parent mean
#'
#' disp = (mp)^(-1/6) - 0.3 + rnorm(1, 0, sd = 0.002)
#' @param mp mid-parent mean.
#' @return simulated dispersion parameter value.
#' @export
#' @examples
#' library(TwoStageLRT)
#' inverse_func(100)

inverse_func <- function(mp){
  disp <- (mp)^(-1/6) - 0.3 + rnorm(1, 0, sd = 0.002)
  return(disp)
}
