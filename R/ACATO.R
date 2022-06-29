#' ACAT function by Yaowu Liu, with changes
#'
#' combine the test results of mBAT and fastBAT by the aggregated Cauchy combination method (ACAT)
#'
#' @param p: integrated in function of mBAT_SVD.R, using P-value of mBAT and fastBAT
#'
#' @return P-value of mBAT-combo
#' @export

ACATO <- function(p){
  if (all(is.na(p))) return(NA)
  p <- p[!is.na(p)]
  p[p == 1] <- 1 - 1e-16  ##this line is just added, April4 by Ang(https://rdrr.io/cran/sumFREGAT/src/R/ACATO.R)
  #### check if there are very small non-zero p values
  is.small <- (p < 1e-16)
  if (sum(is.small) == 0) {
    cct.stat <- sum(tan((0.5 - p) * pi))/length(p)
  } else {
    cct.stat <- sum((1 / p[is.small]) / pi)
    cct.stat <- cct.stat + sum(tan((0.5 - p[!is.small]) * pi))
    cct.stat <- cct.stat/length(p)
  }
  #### check if the test statistic is very large.
  if (cct.stat > 1e+15){
    pval <- (1 / cct.stat) / pi
  } else {
    pval <- 1 - pcauchy(cct.stat)
  }
  pval
}
