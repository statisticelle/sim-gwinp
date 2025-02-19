#' Author: Emma Davies Smith, esmith at hsph.harvard.edu
#' Date: 23 OCT 2024

# Scripts
source('src/help.R')            # helper functions
source('src/genResponse.R')     # fxn to generate responses given scenario
source('src/rmvord2.R')         # faster ordinal data gen fxn (mod from orddata)
source('src/simScenarios.R')    # fxn to simulate given scenarios

# Libraries
library(tidyverse) # dataframe manipulation
library(nlme)      # linear mixed models
library(EnvStats)  # truncated log normal distribution
library(mvtnorm)   # multivariate normal deviates
