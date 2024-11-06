#' Author: Emma Davies Smith, esmith at hsph.harvard.edu
#' Date: 23 OCT 2024

# Scripts
source('src/help.R')            # helper functions
source('src/genResponse.R')     # fxn to generate responses given scenario
source('src/rmvord2.R')         # faster ordinal data gen fxn (mod from orddata)
source('src/simScenarios.R')    # fxn to simulate given scenarios

# Libraries
library(orddata)   # ordinal data generation
library(tidyverse) # dataframe manipulation
library(nlme)      # linear mixed models
library(EnvStats)  # truncated log normal distribution

# orddata library can be installed using:
#install.packages("orddata", repos="http://R-Forge.R-project.org", dependencies=T)
