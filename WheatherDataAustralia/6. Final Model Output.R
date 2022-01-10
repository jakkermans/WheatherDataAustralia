if (!require(foreign)) {install.packages("foreign")}
if (!require(foreign)) {install.packages("mice")}
if (!require(MASS)) {install.packages("MASS")}
if (!require(tseries)) {install.packages("tseries")}
if (!require(Hmisc)) {install.packages("Hmisc")}
if (!require(knitr)) {install.packages("knitr")}
if (!require(e1071)) {install.packages("e1071")}
if (!require(bain)) {install.packages("bain")}
if (!require(jpeg)) {install.packages("jpeg")}

load(file = "Workspaces\\PPC.RData")

#############################################################################################
#Model output
#Run these lines to produce the results described in section 8: Interpretation of estimates
#and intervals
#############################################################################################

output <- matrix(0, nrow = 4, ncol = 4) #create a matrix to store the parameter estimates
rownames(output) <- c("Intercept", "MaxTemp.cen", "MinTemp.cen", "MaxTemp.cen*MinTemp.cen") #change the rownames to the names of the predictors
colnames(output) <- c("Mean", "SE", "Q025", "Q975") #change the column names to the statistics presented in the final table
output[,1] <- apply(full_matrix3[,1:4], 2, mean) #calculate the mean for the intercept and each predictor
output[,2] <- apply(full_matrix3[,1:4], 2, function(x) {sd(x)/sqrt(nrow(wheather.complete))}) #calculate the standard error for the intercept and each predictor
output[,3] <- apply(full_matrix3[,1:4], 2, quantile, 0.025) #calculate the 2.5th percentile for the intercept and each predictor
output[,4] <- apply(full_matrix3[,1:4], 2, quantile, 0.975) #calculate the 97.5th percentile for the intercept and each predictor
kable(output, digits = 5) #print the output table

save.image(file = "Workspaces\\Finaloutput.RData")
