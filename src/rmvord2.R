#' Fast construction of intermediate correlation matrix for mean mapping 
#' Author: Emma Davies Smith, esmith at hsph.harvard.edu
#' Date: 23 OCT 2024
#' 
#' @details
#' To generate data for multivariate cluster trial, must stack participants 
#' within cluster leading to very large correlation matrix. To use mean mapping,
#' must numerically identify intermediate correlation matrix which yields desired
#' correlation matrix upon discretization of normal deviates. Existing rmvord
#' implementation within orddata library was slow due to size of matrices. The 
#' following function exploits block exchangeable correlation structure to 
#' reduce number of redundant intermediate searches and speed up data gen.
#' 
#' 
#' @param n Desired sample size
#' @param probs Bivariate success probabilities
#' @param Cor Desired correlation matrix
#'
#' @return Intermediate correlation structure
getIntermediate <- function (n, probs, Cor) 
{
  q = length(probs[1:4])
  categ_probs = 0
  cumul_probs = list(0)
  quant_probs = list(0)
  means = 0
  vars = 0
  var.wt = function(x, w) {
    m = weighted.mean(x = x, w = w)
    sum((x[1:length(x)] - m)^2 * w[1:length(x)])
  }
  for (i in 1:q) {
    categ_probs[i] = length(probs[[i]])
    cumul_probs[[i]] = cumsum(1:categ_probs[i]/10^12 + probs[[i]])
    cumul_probs[[i]][categ_probs[i]] = 1
    quant_probs[[i]] = qnorm(p = cumul_probs[[i]], mean = 0, 
                             sd = 1)
    means[i] = weighted.mean(x = 1:categ_probs[i], w = probs[[i]])
    vars[i] = var.wt(x = 1:categ_probs[i], w = probs[[i]])
  }
  Cor_norm = Cor
  for (i in 1:(q - 1)) {
    for (j in (i + 1):q) {
      gridd = rep(0, times = 201)
      for (ci in 1:(categ_probs[i] - 1)) {
        for (cj in 1:(categ_probs[j] - 1)) {
          for (steps in -100:100) {
            gridd[101 + steps] = gridd[101 + steps] + 
              pmvnorm(upper = c(quant_probs[[i]][ci], 
                                quant_probs[[j]][cj]), corr = matrix(2, 
                                                                     2, data = c(1, steps/100, steps/100, 
                                                                                 1)))[1]
          }
        }
      }
      f = suppressWarnings(approxfun(y = -100:100/100, 
                                     x = gridd))
      Cor_norm[i, j] = Cor_norm[j, i] = f(Cor[i, j] * sqrt(vars[i] * 
                                                             vars[j]) + means[i] * means[j] - categ_probs[i] * 
                                            categ_probs[j] + categ_probs[j] * sum(cumul_probs[[i]][1:(categ_probs[i] - 
                                                                                                        1)]) + categ_probs[i] * sum(cumul_probs[[j]][1:(categ_probs[j] - 
                                                                                                                                                          1)]))
    }
  }
  
  # Intermediate block matrices
  intOmega <- Cor_norm[1:2, 1:2]
  intPhi   <- Cor_norm[1:2, 3:4]
  
  # Unit matrices
  I <- diag(1, nrow=n, ncol=n)
  J <- matrix(1, nrow=n, ncol=n)
  
  # Intermediate block exchangeable correlation matrix
  intR <- I %x% (intOmega-intPhi) + J %x% intPhi
  
  retval = intR
  retval
}

#' Modified rmvord from orddata to use fast intermediate correlation matrix
#' Author: Emma Davies Smith, esmith at hsph.harvard.edu
#' Date: 23 OCT 2024
#' 
#' @details
#' Simple modification of the rmvord function within the orddata library to use
#' the fast intermediate correlation function above. This is the mean mapping
#' algorithm (i.e., mv normal deviates discretized into ordinal vars).
#' 
#' @param n Desired sample size
#' @param probs Bivariate success probabilities
#' @param Cor Desired correlation matrix
#'
#' @return Clustered multivariate ordinal data
rmvord2 <- function (n = 1, probs, Cor) 
{
  q = length(probs)
  categ_probs = 0
  cumul_probs = list(0)
  quant_probs = list(0)
  means = 0
  vars = 0
  var.wt = function(x, w) {
    m = weighted.mean(x = x, w = w)
    sum((x[1:length(x)] - m)^2 * w[1:length(x)])
  }
  for (i in 1:q) {
    categ_probs[i] = length(probs[[i]])
    cumul_probs[[i]] = cumsum(1:categ_probs[i]/10^12 + probs[[i]])
    cumul_probs[[i]][categ_probs[i]] = 1
    quant_probs[[i]] = qnorm(p = cumul_probs[[i]], mean = 0, 
                             sd = 1)
    means[i] = weighted.mean(x = 1:categ_probs[i], w = probs[[i]])
    vars[i] = var.wt(x = 1:categ_probs[i], w = probs[[i]])
  }
  
  Cor_norm <- getIntermediate(n=length(probs)/2, probs=probs, Cor=Cor)
  retval = rmvnorm(n = n, sigma = Cor_norm)
  for (i in 1:q) {
    retval[, i] = cut(x = retval[, i], breaks = c(-1/0, quant_probs[[i]]), 
                      right = FALSE)
  }
  retval
}