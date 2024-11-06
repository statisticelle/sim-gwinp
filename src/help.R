#' Author: Emma Davies Smith, esmith at hsph.harvard.edu
#' Date: 23 OCT 2024

#' Logit transformation
logit <- function(x){
  log(x/(1-x))
}

#' Expit transformation (logit back-transformation)
expit <- function(x){
  exp(x)/(1+exp(x))
}

#' Check true value of parameter falls within confidence limits
checkCoverage <- function(theta, ci){
  ifelse(theta >= ci[1] & theta <= ci[2], 1, 0)
}

#' Check null value (0.5) falls outside confidence limits
checkReject <- function(ci){
  ifelse(ci[1] <= 0.5 & ci[2] >= 0.5, 0, 1)
}

# Check if truth falls to left of lower limit
checkLeft <- function(theta, ci){
  ifelse(theta <= ci[1], 1, 0)
}

# Check if truth falls to right of upper limit
checkRight <- function(theta, ci){
  ifelse(theta >= ci[2], 1, 0)
}

# Calculate value of win probability given binomial params in each arm
calcOrdWinP <- function(k, pi0, pi1){
  # Control distribution function
  F0 <- dbinom(x=0:(k-1), size=(k-1), prob=pi0)
  
  # See J.H. Klotz (JASA 1966): The Wilcoxon, Ties, and the Computer
  Lambda <- diag(x=0.5, nrow=k)
  Lambda[upper.tri(Lambda)] <- 1
  
  F1 <- dbinom(x=0:(k-1), size=(k-1), prob=pi1)
  t(F0) %*% Lambda %*% F1
}


#' Construct target block exhangeable correlation matrix
#' @details
#' For more details on block exchangeable correlation matrix see:
#' Wang et al. (CCT 2021) A flexible sample size solution for longitudinal and
#' crossover cluster randomized trials with continuous outcomes.
#' 
#' @param n   Number of subjects per cluster
#' @param w12 Within-subject between-endpoint correlation
#' @param p11 Intracluster correlation of 1st endpoint
#' @param p22 Intracluster correlation of 2nd endpoint
#' @param p12 Within-cluster between-subject between-endpoint correlation
#'
#' @return Target correlation matrix R
#'
#' @examples
#' getCorrelation(n=3, w12=0.3, p11=0.10, p22=0.05, p12=0.025)
getCorrelation <- function(n, w12, p11, p22, p12){
  
  # Endpoint correlation
  Omega <- matrix(c(1, w12, w12, 1), ncol=2)
  
  # Intracluster correlations
  Phi <- matrix(c(p11, p12, p12, p22), ncol=2)
  
  # Unit matrices
  I <- diag(1, nrow=n, ncol=n)
  J <- matrix(1, nrow=n, ncol=n)
  
  # Desired block exchangeable correlation matrix
  R <- I %x% (Omega-Phi) + J %x% Phi
  return(R)
}

#' Get treatment arm parameters to obtain desired endpoint win probability
#' @details
#' For details on calculation see:  
#' J.H. Klotz (JASA 1966) The Wilcoxon, Ties, and the Computer
#' 
#' @param k Total number of categories for ordinal variable (fixed)
#' @param pi0 Success probability in control arm (fixed)
#' @param theta Desired win probability 
#'
#' @return Success probability in treatment arm
#'
#' @examples
#' getParams(k=5, pi0=0.3, theta=0.5)  # return should be equal to pi0
#' getParams(k=7, pi0=0.5, theta=0.64) # return should be equal to ~0.605
getParams <- function(k, pi0, theta){
  
  # Control distribution function
  F0 <- dbinom(x=0:(k-1), size=(k-1), prob=pi0)
  
  # See J.H. Klotz (JASA 1966): The Wilcoxon, Ties, and the Computer
  Lambda <- diag(x=0.5, nrow=k)
  Lambda[upper.tri(Lambda)] <- 1
  
  # Find value of pi1 which produces win probability (theta) closest to target
  pi_eval <- seq(0, 1, 0.0001)
  diff <- c()
  for (pi in pi_eval){
    F1 <- dbinom(x=0:(k-1), size=(k-1), prob=pi)
    diff <- c(diff, abs(t(F0) %*% Lambda %*% F1 - theta))
  }
  
  pi1 <- pi_eval[which.min(diff)]
  return(pi1)
}

