if (!require(foreign)) {install.packages("foreign")}
if (!require(foreign)) {install.packages("mice")}
if (!require(MASS)) {install.packages("MASS")}
if (!require(tseries)) {install.packages("tseries")}
if (!require(Hmisc)) {install.packages("Hmisc")}
if (!require(knitr)) {install.packages("knitr")}
if (!require(e1071)) {install.packages("e1071")}
if (!require(bain)) {install.packages("bain")}

##############################################
# Preparing the data
##############################################

set.seed(123)
wheatherdat <- read.csv("weatherAUS.csv", sep = ",")
used_indices <- sample(seq(1,nrow(wheatherdat),1),2000) #sample 2000 indices to be used for the project
wheatherdat2 <- wheatherdat[used_indices,c(3,4,5,9)] #sample observations using the 2000 sampled indices

imp <- mice(wheatherdat2, m = 5)
wheather.complete <- complete(imp)
wheather.complete$MinTemp.cen <- scale(wheather.complete$MinTemp, center = TRUE, scale = FALSE)
wheather.complete$MaxTemp.cen <- scale(wheather.complete$MaxTemp, center = TRUE, scale = FALSE)
wheather.complete$WindGustSpeed.cen <- scale(wheather.complete$WindGustSpeed, center = TRUE, scale = FALSE)

save.image(file = "Workspaces\\ImputedData.RData")
