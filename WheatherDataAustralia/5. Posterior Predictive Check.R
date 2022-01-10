if (!require(foreign)) {install.packages("foreign")}
if (!require(foreign)) {install.packages("mice")}
if (!require(MASS)) {install.packages("MASS")}
if (!require(tseries)) {install.packages("tseries")}
if (!require(Hmisc)) {install.packages("Hmisc")}
if (!require(knitr)) {install.packages("knitr")}
if (!require(e1071)) {install.packages("e1071")}
if (!require(bain)) {install.packages("bain")}
if (!require(jpeg)) {install.packages("jpeg")}

load(file = "Workspaces\\Convergence.RData")
set.seed(123)

#############################################################################################
#Posterior predictive p-value
#Run these lines to produce the results described in section 7: Posterior predictive check
#############################################################################################

sim_datasets <- array(data = NA, dim = c(nrow(full_matrix), nrow(wheather.complete))) #create space for the simulated datasets
obs_datasets <- array(data = NA, dim = c(nrow(full_matrix), nrow(wheather.complete))) #create space for the observed datasets
for (i in 1:nrow(sim_datasets)) {
  sim_datasets[i,] <- rnorm(ncol(sim_datasets), mean = full_matrix3[i,1] + full_matrix3[i,2]*wheather.complete$MaxTemp.cen + full_matrix3[i,3]*wheather.complete$MinTemp.cen + full_matrix3[i,4]*(wheather.complete$MaxTemp.cen*wheather.complete$MinTemp.cen), sd = sqrt(full_matrix3[i,5])) #create simulated datasets by sampling values from a normal distribution where the mean and variance are based on the values of the current iteration
}

for (i in 1:nrow(sim_datasets)) {
  sim_datasets[i,] <- sim_datasets[i,] - full_matrix3[i,1] - full_matrix3[i,2]*wheather.complete$MaxTemp.cen - full_matrix3[i,3]*wheather.complete$MinTemp.cen - full_matrix3[i,4]*(wheather.complete$MaxTemp.cen*wheather.complete$MinTemp.cen) #calculate the residuals for each simulated dataset
}

for (i in 1:nrow(obs_datasets)) {
  obs_datasets[i,] <- wheather.complete$Rainfall - full_matrix3[i,1] - full_matrix3[i,2]*wheather.complete$MaxTemp.cen - full_matrix3[i,3]*wheather.complete$MinTemp.cen - full_matrix3[i,4]*(wheather.complete$MaxTemp.cen*wheather.complete$MinTemp.cen) #for each set of regression coefficients, calculate the residuals for an observed datasets using those coefficients
}

sd_dataset <- rep(0,nrow(sim_datasets)) #create space for the simulated test statistics
obs_sd <- rep(0, nrow(obs_datasets)) #create space for the observed test statistics

for (j in 1:nrow(sim_datasets)){
  sd_dataset[j] <-  (mean(sim_datasets[j,] <= mean(sim_datasets[j,]) - 2*sd(sim_datasets[j,])) + mean(sim_datasets[j,] >= mean(sim_datasets[j,]) + 2*sd(sim_datasets[j,]))) <= 0.05 #if there are less than 5% or precisely 5% of the residuals in the tails, assign that simulated residuals dataset a 1, otherwise 0
  obs_sd[j] <- (mean(obs_datasets[j,] <= mean(obs_datasets[j,]) - 2*sd(obs_datasets[j,])) + mean(obs_datasets[j,] >= mean(obs_datasets[j,]) + 2*sd(obs_datasets[j,]))) <= 0.05 #if there are less than 5% or precisely 5% of the residuals in the tails, assign that observed residuals dataset a 1, otherwise 0
}

mean(ifelse(sd_dataset  >= obs_sd, 1, 0)) #calculate the proportion of cases where the simulated test statistic is equal to or greater than the observed test statistic as the posterior predictive p-value

save.image(file = "Workspaces\\PPC.RData")
