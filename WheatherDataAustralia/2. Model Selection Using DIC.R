if (!require(foreign)) {install.packages("foreign")}
if (!require(foreign)) {install.packages("mice")}
if (!require(MASS)) {install.packages("MASS")}
if (!require(tseries)) {install.packages("tseries")}
if (!require(Hmisc)) {install.packages("Hmisc")}
if (!require(knitr)) {install.packages("knitr")}
if (!require(e1071)) {install.packages("e1071")}
if (!require(bain)) {install.packages("bain")}

#############################################################################################
# Functions
#############################################################################################

gibbs_sampler <- function(pred1, pred2, out, n.iter, n.chains, burn_in) {
  ##This function conducts Gibb's sampling based on a regression model with one outcome variable and 2 predictors. #It samples new values for the regression coefficients in each iteration and subsequently stores those in a matrix. For each chain, the full matrix is stored so that convergence can be evaluated afterwards.
  
  #create storage for the results of the different chains
  matrices <- list()
  
  #Specify the initial values
  initial_betas <- list(c(2.4,-0.10,0.09), c(-2,0.30,-0.25))  
  m0.prior <- 0
  m1.prior <- 0
  m2.prior <- 0
  s0.prior <- 1000
  s1.prior <- 1000
  s2.prior <- 1000
  a.prior <- 0.001
  b.prior <- 0.001
  tm.prior <- 0
  ts.prior <- 1000
  df <- 1
  
  for (chain in 1:n.chains) {
    #create storage for the results of the current chain
    
    #Assign the initial values to the parameters
    var <- rep(0,n.iter)
    b0 <- rep(0,n.iter)
    b1 <- rep(0,n.iter)
    b2 <- rep(0,n.iter)
    var[1] <- 1
    b0[1] <- initial_betas[[chain]][1]
    b1[1] <- initial_betas[[chain]][2]
    b2[1] <- initial_betas[[chain]][3]
    
    for (i in 2:n.iter) {
      m0.post <- ((sum(out - b1[i-1]*pred1 - b2[i-1]*pred2)/var[i-1]) + (m0.prior/s0.prior)) / ((length(out) / var[i-1]) + (1 / s0.prior)) #Posterior mean of b0
      s0.post <- 1 / ((length(out) / var[i-1]) + (1 / s0.prior)) #Posterior variance of b0
      b0[i] <- rnorm(1, m0.post, sqrt(s0.post)) #Sample new value for b0
      
      m1.post <- ((sum(pred1*(out - b0[i] - b2[i-1]*pred2)) / var[i-1]) + (m1.prior / s1.prior)) / ((sum(pred1^2) / var[i-1]) + (1 / s1.prior)) #Posterior mean of b1
      s1.post <- 1 / ((sum(pred1^2) / var[i-1]) + (1 / s1.prior)) #Posterior variance of b1
      b1[i] <- rnorm(1, m1.post, sqrt(s1.post)) #Sample new value for b1
      
      #m2.post <- ((sum(pred2*(out - b0[i] - b1[i]*pred1)) / var[i-1]) + (m2.prior / s2.prior)) / ((sum(pred2^2) / var[i-1]) + (1 / s2.prior)) #Posterior mean of b2
      #s2.post <- 1 / ((sum(pred2^2) / var[i-1]) + (1 / s2.prior)) #Posterior variance of b2
      #b2[i] <- rnorm(1, m2.post, sqrt(s2.post)) #Sample new value for b2
      b2_star <- rnorm(1, mean = b2[i-1], sd = 0.1)
      u <- runif(1,0,1)
      
      b2_post <- exp(-(b2[i-1]^2)*(sum(pred2^2)/(2*var[i-1]))) + (b2[i-1]*(sum(pred2*(out - b0[i] - b1[i]*pred1))/var[i-1])) * (1 + ((b2[i-1] - tm.prior)^2/(df*(ts.prior^2))))^(-((df+1)/2))
      
      b2star_post <- exp((-(b2_star^2)*(sum(pred2^2)/(2*var[i-1]))) + (b2_star*(sum(pred2*(out - b0[i] - b1[i]*pred1))/var[i-1]))) * (1 + ((b2_star - tm.prior)^2/(df*(ts.prior^2))))^(-((df+1)/2))
      
      r <- (b2star_post/b2_post)*(b2[i-1]/b2_star)
      b2[i] <- ifelse(u <= r, b2_star, b2[i-1])
      
      a.post <- (length(out) / 2) + a.prior #Posterior shape of variance
      b.post <- b.prior + (sum((out - b0[i] + b1[i]*pred1 + b2[i]*pred2)^2) / 2) #Posterior rate of variance
      var[i] <- 1/rgamma(1, shape = a.post, rate = b.post) #Sample new value for variance
      
    }
    matrices[[chain]] <- cbind(b0[-(1:burn_in)], b1[-(1:burn_in)], b2[-(1:burn_in)], var[-(1:burn_in)]) #Store all the sampled values for b0, b1 and b2 of the current chain but remove the burn-in period
  }
  return(matrices) #Return the sampled values for b0, b1, and b2 of the different chains
}

gibbs_sampler2 <- function(pred1, pred2, out, n.iter, n.chains, burn_in) {
  ##This function conducts Gibb's sampling based on a regression model with one outcome variable and 2 predictors. #It samples new values for the regression coefficients in each iteration and subsequently stores those in a matrix. For each chain, the full matrix is stored so that convergence can be evaluated afterwards.
  
  #create storage for the results of the different chains
  matrices <- list()
  
  #Specify the initial values
  initial_betas <- list(c(2.4,-0.10,0.09, -0.009), c(-2,0.30,-0.25, 0.05))  
  m0.prior <- 0
  m1.prior <- 0
  m2.prior <- 0
  m3.prior <- 0
  s0.prior <- 1000
  s1.prior <- 1000
  s2.prior <- 1000
  s3.prior <- 1000
  a.prior <- 0.001
  b.prior <- 0.001
  tm.prior <- 0
  ts.prior <- 1000
  df <- 1
  
  for (chain in 1:n.chains) {
    #create storage for the results of the current chain
    
    #Assign the initial values to the parameters
    var <- rep(0,n.iter)
    b0 <- rep(0,n.iter)
    b1 <- rep(0,n.iter)
    b2 <- rep(0,n.iter)
    b3 <- rep(0,n.iter)
    var[1] <- 1
    b0[1] <- initial_betas[[chain]][1]
    b1[1] <- initial_betas[[chain]][2]
    b2[1] <- initial_betas[[chain]][3]
    b3[1] <- initial_betas[[chain]][4]
    
    for (i in 2:n.iter) {
      m0.post <- ((sum(out - b1[i-1]*pred1 - b2[i-1]*pred2 - b3[i-1]*(pred1*pred2))/var[i-1]) + (m0.prior/s0.prior)) / ((length(out) / var[i-1]) + (1 / s0.prior)) #Posterior mean of b0
      s0.post <- 1 / ((length(out) / var[i-1]) + (1 / s0.prior)) #Posterior variance of b0
      b0[i] <- rnorm(1, m0.post, sqrt(s0.post)) #Sample new value for b0
      
      m1.post <- ((sum(pred1*(out - b0[i] - b2[i-1]*pred2 - b3[i-1]*(pred1*pred2))) / var[i-1]) + (m1.prior / s1.prior)) / ((sum(pred1^2) / var[i-1]) + (1 / s1.prior)) #Posterior mean of b1
      s1.post <- 1 / ((sum(pred1^2) / var[i-1]) + (1 / s1.prior)) #Posterior variance of b1
      b1[i] <- rnorm(1, m1.post, sqrt(s1.post)) #Sample new value for b1
      
      m2.post <- ((sum(pred2*(out - b0[i] - b1[i]*pred1)) / var[i-1]) + (m2.prior / s2.prior)) / ((sum(pred2^2) / var[i-1]) + (1 / s2.prior)) #Posterior mean of b2
      s2.post <- 1 / ((sum(pred2^2) / var[i-1]) + (1 / s2.prior)) #Posterior variance of b2
      b2[i] <- rnorm(1, m2.post, sqrt(s2.post)) #Sample new value for b2
      
      m3.post <- ((sum((pred1*pred2)*(out - b0[i-1] - b1[i-1]*pred1 - b2[i-1]*pred2))/var[i-1]) + (m3.prior/s3.prior)) / ((sum((pred1*pred2)^2)/var[i-1]) + (1/s3.prior)) #Posterior mean of b3
      s3.post <- 1 / ((sum((pred1*pred2)^2)/var[i-1]) + (1/s3.prior)) #Posterior variance of b3
      b3[i] <- rnorm(1, m3.post, sqrt(s3.post)) #Sample new value for b3
      
      a.post <- (length(out) / 2) + a.prior #Posterior shape of variance
      b.post <- b.prior + (sum((out - b0[i] + b1[i]*pred1 + b2[i]*pred2 + b3[i-1]*(pred1*pred2))^2) / 2) #Posterior rate of variance
      var[i] <- 1/rgamma(1, shape = a.post, rate = b.post) #Sample new value for variance
      
    }
    matrices[[chain]] <- cbind(b0[-(1:burn_in)], b1[-(1:burn_in)], b2[-(1:burn_in)], b3[-(1:burn_in)], var[-(1:burn_in)]) #Store all the sampled values for b0, b1, b2 and b3 of the current chain but remove the burn-in period
  }
  return(matrices) #Return the sampled values for b0, b1, and b2 of the different chains
}

dic_calculation <- function(full_matrix, out, pred1, pred2) {
  obs_mean <- colMeans(full_matrix)[1] + colMeans(full_matrix)[2]*pred1 + colMeans(full_matrix)[3]*pred2 #calculate mean for the normal density to be sampled from
  dhat <- -2*sum(dnorm(out, mean = obs_mean, sd = sqrt(colMeans(full_matrix)[4]), log = T)) #sample from a normal density with the mean calculated on the previous line and as standard deviation the square root of the mean variance 
  dbars <- rep(0,nrow(full_matrix)) #create space for the likelihoods of the different Gibbs sampler iterations
  for (i in 1:nrow(full_matrix)) { 
    new_mean <- full_matrix[i,1] + full_matrix[i,2]*pred1 + full_matrix[i,3]*pred2 #calculate mean for the normal density to be sampled from based on the coefficients of the current iteration
    dbars[i] <- -2*sum(dnorm(out, mean = new_mean, sd = sqrt(full_matrix[i,4]), log = T)) #sample from a normal density with the mean calculated on the previous line and as standard deviation the square root of the variance of the current iteration 
  }
  pd <- abs(mean(dbars) - dhat) #calculate the estimated number of effective parameters
  DIC <- dhat + 2*pd #calculate the DIC
  return(DIC) #return the DIC
}

dic_calculation_interaction <- function(full_matrix, out, pred1, pred2, pred3) {
  obs_mean <- colMeans(full_matrix)[1] + colMeans(full_matrix)[2]*pred1 + colMeans(full_matrix)[3]*pred2 + colMeans(full_matrix)[4]*(pred1*pred2)# calculate mean for the normal density to be sampled from
  dhat <- -2*sum(dnorm(out, mean = obs_mean, sd = sqrt(colMeans(full_matrix)[5]), log = T)) #sample from a normal density with the mean calculated on the previous line and as standard deviation the square root of the mean variance 
  dbars <- rep(0,nrow(full_matrix)) #create space for the likelihoods of the different Gibbs sampler iterations
  for (i in 1:nrow(full_matrix)) {
    new_mean <- full_matrix[i,1] + full_matrix[i,2]*pred1 + full_matrix[i,3]*pred2 + full_matrix[i,4]*(pred1*pred2) #calculate mean for the normal density to be sampled from based on the coefficients of the current iteration
    dbars[i] <- -2*sum(dnorm(out, mean = new_mean, sd = sqrt(full_matrix[i,5]), log = T)) #sample from a normal density with the mean calculated on the previous line and as standard deviation the square root of the variance of the current iteration 
  }
  pd <- abs(mean(dbars) - dhat) #calculate the estimated number of effective parameters
  DIC <- dhat + 2*pd #calculate the DIC
  return(DIC) #calculate the DIC
}

#############################################################################################
#Model selection using DIC
#Run these lines of code to reproduce the results under section 4: Model selection using 
#the DIC
#############################################################################################

load(file = "Workspaces\\ImputedData.RData")

n.iters <- 10000 #set the number of iterations
samples2 <- gibbs_sampler(pred1 = wheather.complete$MaxTemp.cen, pred2 = wheather.complete$WindGustSpeed.cen, out = wheather.complete$Rainfall, n.iter = n.iters, n.chains = 2, burn_in = 1000) #sample values from the self-programmed Gibbs sampler using MaxTemp.cen and WindGustSpeed.cen as predictors
full_matrix <- rbind(samples2[[1]], samples2[[2]]) #row-bind the values of the different chains
dic_model1 <- dic_calculation(full_matrix, wheather.complete$Rainfall, wheather.complete$MaxTemp.cen, wheather.complete$WindGustSpeed.cen) #calculate the DIC

samples3 <- gibbs_sampler(pred1 = wheather.complete$MinTemp.cen, pred2 = wheather.complete$WindGustSpeed.cen, out = wheather.complete$Rainfall, n.iter = 10000, n.chains = 2, burn_in = 1000) #sample values from the self-programmed Gibbs sampler using MinTemp.cen and WindGustSpeed.cen as predictors
full_matrix2 <- rbind(samples3[[1]], samples3[[2]]) #row-bind the values of the different chains
dic_model2 <- dic_calculation(full_matrix2, out = wheather.complete$Rainfall, pred1 = wheather.complete$MinTemp.cen, pred2 = wheather.complete$WindGustSpeed.cen) #calculate the DIC

samples4 <- gibbs_sampler2(pred1 = wheather.complete$MaxTemp.cen, pred2 = wheather.complete$MinTemp.cen, out = wheather.complete$Rainfall, n.iter = n.iters, n.chains = 2, burn_in = 1000) #sample values from the self-programmed Gibbs sampler using MaxTemp.cen, MinTemp.cen and an interaction between them as predictors
full_matrix3 <- rbind(samples4[[1]], samples4[[2]]) #row-bind the values of the different chains
colnames(full_matrix3) <- c("b0", "b1", "b2", "b3", "var")
dic_model3 <- dic_calculation_interaction(full_matrix3, out = wheather.complete$Rainfall, pred1 = wheather.complete$MaxTemp.cen, pred2 = wheather.complete$MinTemp.cen) #calculate the DIC

dic_matrix <- matrix(0, nrow = 3, ncol = 1) #create a matrix to store the DIC values of the different models
colnames(dic_matrix) <- c("DIC") #set the name of the column to 'DIC'
rownames(dic_matrix) <- c("Model1: MaxTemp.cen & WindGustSpeed.cen", "Model2: MinTemp.cen & WindGustSpeed.cen", "Model3: MaxTemp.cen & MinTemp.cen & MaxTemp.cen*MinTemp.cen") #give the rows useful names so that the different compared models are clear
dic_matrix[1,] <- dic_model1 #store the DIC of the first imposed model
dic_matrix[2,] <- dic_model2 #store the DIC of the second imposed model
dic_matrix[3,] <- dic_model3 #store the DIC of the third imposed model
kable(dic_matrix, digits = 2) #print the matrix

save.image(file = "Workspaces\\DICSelection.RData")
