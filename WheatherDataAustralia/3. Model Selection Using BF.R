if (!require(foreign)) {install.packages("foreign")}
if (!require(foreign)) {install.packages("mice")}
if (!require(MASS)) {install.packages("MASS")}
if (!require(tseries)) {install.packages("tseries")}
if (!require(Hmisc)) {install.packages("Hmisc")}
if (!require(knitr)) {install.packages("knitr")}
if (!require(e1071)) {install.packages("e1071")}
if (!require(bain)) {install.packages("bain")}

load(file = "Workspaces\\DICSelection.RData")

#############################################################################################
#Model selection using Bayes Factor
#Run these lines of code to reproduce the results under section 5: Model selection using
#the Bayes Factor
#############################################################################################
test_data <- wheather.complete
test_data$MaxTemp.cen <- -1*test_data$MaxTemp.cen
bf_lm1 <- lm(Rainfall ~ MaxTemp.cen + WindGustSpeed.cen, data = test_data)
bf_lm2 <- lm(Rainfall ~ MinTemp.cen + WindGustSpeed.cen, data = test_data)
bf_lm3 <- lm(Rainfall ~ MaxTemp.cen*MinTemp.cen, data = test_data)
summary(bf_lm1)
summary(bf_lm2)
summary(bf_lm3)

(hyp_1 <- bain(bf_lm1, 'MaxTemp.cen > WindGustSpeed.cen'))
(hyp_2 <- bain(bf_lm2, 'MinTemp.cen > WindGustSpeed.cen'))
(hyp_3 <- bain(bf_lm3, 'MaxTemp.cen:MinTemp.cen > MaxTemp.cen*MinTemp.cen'))

save.image(file = "Workspaces\\BFSelection.RData")
