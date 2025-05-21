# sim-gwinp
R code to reproduce simulation study of global win probability for cluster randomized trials with multiple endpoints of different scales

Corresponding Paper: 
[Davies Smith E, Jairath V, Zou G. Rank-based estimators of global treatment effects for cluster randomized trials with multiple endpoints on different scales. Statistical Methods in Medical Research. 2025;0(0).](https://journals.sagepub.com/doi/abs/10.1177/09622802251338387)

Steps for running the simulation study are as follows:
1. dependencies.R - Load simulation scripts and other dependencies.
2. scenarios.R - Input parameter values of interest and generate all scenarios for data generation/investigation. Save to scenarios.csv. 
3. main.R - Import scenarios.csv and simulate all scenarios using simScenarios() function.

Scripts within src/
1. help.R - Helper functions (e.g., logit, correlation matrix construction, identification of parameter values for given win probability).
2. rmvord2.R - Modifications of rmvord function from orddata package to speed up data generation in presence of block exchangeable correlation structure.
3. genResponse.R - Function to generate bivariate ordinal responses given simulation scenario parameters.
4. simScenarios.R - Function which runs all simulation scenarios, and returns results. 

Files within out/
1. scenarios.csv - Dataframe of scenarios output by scenarios.R
2. simResults.csv - Dataframe of results output by simScenarios() within main.R
