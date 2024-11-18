#' Generate bivariate ordinal responses for all clusters and replicates
#' Author: Emma Davies Smith, esmith at hsph.harvard.edu
#' Date:   23 OCT 2024
#' 
#' @details
#' Correlated ordinal responses generated using mean mapping procedure discussed
#' by Kaiser, Trager, and Leisch (2011). i.e. by discretizing mv normal deviates.
#' See: https://epub.ub.uni-muenchen.de/12157/1/kaiser-tr-94-ordinal.pdf
#' 
#' @param c Number of clusters per replicate
#' @param n Number of subjects per cluster
#' @param B Number of replicates
#' @param R Desired correlation matrix
#' @param pi1 Success probability for first ordinal endpoint (5-class)
#' @param pi2 Success probability for second ordinal endpoint (7-class)
#'
#' @return Dataframe with one bivariate record per row
genResponse <- function(C, n, n_type, B, R, pi1, pi2){
  
  # Adjustments for unequal sample size
  if (n_type == 'Unequal'){
    n_set <- n   # desired cluster size
    n <- 4*n_set     # generate more subjects than needed, then delete at random
    
    # Re-size correlation matrix accordingly 
    Omega <- R[1:2, 1:2]
    Phi   <- R[1:2, 3:4]
    I <- diag(1, nrow=n, ncol=n)
    J <- matrix(1, nrow=n, ncol=n)
    R <- I %x% (Omega-Phi) + J %x% Phi
  }
  
  # Set marginal distributions for each of two endpoints
  marginalsOne <- list(
    X1 = dbinom(x=0:4, size=4, prob=pi1),
    X2 = dbinom(x=0:6, size=6, prob=pi2)
  )
  
  # To generate clustered data, must stack subjects
  marginalsAll <- unlist(rep(list(marginalsOne), n), recursive=F)
  
  # Use modified rmvord from orddata package to generate data (mean mapping alg)
  # Each cluster generated as a row, bivariate participant responses as cols
  X <- rmvord2(
                n=C*B, # Generate C*B total clusters, or C clusters/rep
                probs=marginalsAll, # According to provide marginals
                Cor=R # And target correlation matrix
              )
  
  # Column pairs represent bivariate response for subject within cluster
  # Need to separate into one row per subject
  odd_ix   <- seq(1, ncol(X), 2)
  even_ix  <- seq(2, ncol(X), 2)
  
  # Generated outcomes from 1 to 5, rather than 0 to 4... shifting here.
  X1 <- X[ ,odd_ix]-1
  X2 <- X[,even_ix]-1
  
  # Turn into vectors (in row-major order)
  X1 <- unlist(as.list(t(X1)))
  X2 <- unlist(as.list(t(X2)))
  
  # Combine into labeled data frame
  tX <- data.frame(cbind(X1, X2))
  tX$B <- rep(1:ceiling(nrow(tX)/(n*C)), each=(n*C))
  tX$C <- rep(rep(1:C, each=n), B)
  
  # Add participant IDs
  tX <- tX %>%
    group_by(B, C) %>%
    mutate(pid=1:n)
  
  if(n_type=='Unequal'){
    # generate cluster size from truncated log-normal distribution
    # cv = 0.65; mu = desired n
    # n truncated below by 3 and above by 4*mu (rarely >=100)
    mu <- n_set
    var <- (mu*0.65)^2
    mu_ln <- log(mu^2/sqrt(mu^2+var))
    var_ln <- log(1+var/mu^2)
    
    # Generate number of subjects per cluster at random for each rep + cluster
    keep <- expand.grid(B=1:B, C=1:C)
    nkeep <- nrow(keep)
    keep$n <- ceiling(
      rlnormTrunc(n=nkeep, 
                  meanlog=mu_ln, sdlog=sqrt(var_ln), 
                  min=3, max=4*mu)
    )
  
    # Randomly sample subject ids to retain (out of original 100)
    keep <- keep %>%
      group_by(B, C) %>%
      mutate(pid=list(sample(1:(4*mu), size=n, replace=F))) %>% 
      unnest(pid) %>%
      select(-n)
    
    # Return responses for sampled ids
    keep <- keep %>% arrange(B, C, pid)
    tX <- tX %>% inner_join(keep, by=c('B'='B', 'C'='C', 'pid'='pid'))
  }
  
  return(tX)
  
}
