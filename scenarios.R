#' Author: Emma Davies Smith, esmith at hsph.harvard.edu
#' Date: 23 OCT 2024

# run dependencies.R first or load following library to run
# library(tidyverse) 

# Parameters -------------------------------------------------------------------

# Proportion of clusters allocated to treatment arm
r <- c(0.5, 0.7)

# Total number of clusters
C <- c(10, 20, 30, 50)

# Average number of subjects per cluster
n <- c(10, 30)

# Cluster size distribution
n_type <- c('Equal', 'Unequal')

# Within-subject between-endpoint correlation
w12 <- c(0.2, 0.5, 0.8)

# Intracluster correlations
p11 <- 0.10
p22 <- 0.05
p12 <- 0.025

# Global win probability (null, small, moderate, large)
theta <- c(0.50, 0.56, 0.64, 0.71)

# Effect heterogeneity (difference in endpoint effects)
heterog <- c('None', 'Large')

# Control success probabilities for each endpoint
pi01 <- c(0.3, 0.5)
pi02 <- c(0.3, 0.5)

# Number of simulation replicates
B <- 5000

# Scenarios --------------------------------------------------------------------

# Generates all scenarios based on input parameter values
scenarios <- expand.grid(
    B=B, r=r, C=C, n=n, n_type=n_type,
    w12=w12, p11=p11, p22=p22, p12=p12,
    theta=theta, heterog=heterog, 
    pi01=pi01, pi02=pi02
  ) %>% 
  # Due to large number of scenarios, restrict to same event probability
  filter(pi01==pi02) %>%
  # Effect heterogeneity
  mutate(theta1=case_when(heterog=='Large'    ~ theta-0.105,
                          heterog=='Moderate' ~ theta-0.07,
                          heterog=='Small'    ~ theta-0.03,
                          heterog=='None'     ~ theta),
         theta2=case_when(heterog=='Large'    ~ theta+0.105,
                          heterog=='Moderate' ~ theta+0.07,
                          heterog=='Small'    ~ theta+0.03,
                          heterog=='None'     ~ theta)) 

# Attach treatment arm parameters to obtain desired endpoint effects (slow)
scenarios$pi11 <- mapply(FUN=function(x,y) getParams(k=5, pi0=x, theta=y),
                         x=scenarios$pi01, y=scenarios$theta1)

scenarios$pi12 <- mapply(FUN=function(x,y) getParams(k=7, pi0=x, theta=y),
                         x=scenarios$pi02, y=scenarios$theta2)

# Arrange scenarios and add unique IDs
scenarios <- scenarios %>%
  filter(pi01==pi02) %>%
  arrange(r, n_type, heterog, theta, C, n, w12, pi01, pi02) %>%
  mutate(sid=row_number()) %>%
  select(sid, r, n_type, heterog, theta, theta1, theta2, C, n, 
         w12, p11, p22, p12, pi01, pi11, pi02, pi12, B)

# Print total number of scenarios and preview
cat('Total scenarios:', nrow(scenarios), '\n')
head(scenarios)

# Output scenarios to dataframe for use in simulation
write.csv(scenarios, './out/scenarios.csv', row.names=F)
