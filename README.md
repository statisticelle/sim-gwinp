# sim-gwinp
R code to reproduce simulation study of global win probability for cluster randomized trials with multiple endpoints of ordinal or different types

Corresponding Pre-print (will update with pub once available): 
[Davies Smith, E., Jairath, V., & Zou, G. (2024). Rank-based estimators of global treatment effects for cluster randomized trials with multiple endpoints. arXiv preprint arXiv:2401.12420.](https://arxiv.org/abs/2401.12420)

Steps for running the simulation study are as follows:
1. scenarios.R - Input parameter values of interest and generate all scenarios for data generation/investigation. Save to scenarios.csv. 
2. main.R - Load simulation scripts and other dependencies (stored in dependencies.R), import scenarios.csv, and simulate all scenarios using simScenarios() function.

Scripts within src/
1. help.R - Helper functions (e.g., logit, correlation matrix construction, identification of parameter values for given win probability).
2. rmvord2.R - Modifications of rmvord function from orddata package to speed up data generation in presence of block exchangeable correlation structure.
3. genResponse.R - Function to generate bivariate ordinal responses given simulation scenario parameters.
4. simScenarios.R - Function which runs all simulation scenarios, and returns results. 

Files within out/
1. scenarios.csv - Dataframe of scenarios output by scenarios.R
2. simResults.csv - Dataframe of results output by simScenarios() within main.R
