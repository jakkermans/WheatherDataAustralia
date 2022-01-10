# WheatherDataAustralia

This repository belongs to the WheatherDataAustralia project. It contains the manuscript (both in pdf and RMarkdown), references used, plots used in the manuscript, the actual data, and all the intermediate workspaces. In order to get the results in the manuscript and to enable a reproducible workflow, please follow the following order of R scripts:

File                                    Purpose
-----                                   --------
1. Preparing & Imputing Data.R          Pre-processing the data, imputing missing values and sampling observations to be used.
2. Model Selection Using DIC.R          Selection of the best regression model using the DIC.
3. Model Selection Using BF.R           Selection of the best regression model using the Bayes Factor.
4. Convergence.R                        Inspection of the convergence of the sampled values using the Gibbs sampler.
5. Posterior Predictive Check.R         Evaluation whether the residuals of the sampled values are normally distributed.
6. Final Model Output.R                 Print the output of the final model as presented in the manuscript.
-----                                   --------

Each R script will save the current image in the Workspace directory and that image is subsequently loaded in the next R script. Generated images are stored in the Plots directory.
