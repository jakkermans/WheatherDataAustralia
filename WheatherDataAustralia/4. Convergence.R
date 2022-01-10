if (!require(foreign)) {install.packages("foreign")}
if (!require(foreign)) {install.packages("mice")}
if (!require(MASS)) {install.packages("MASS")}
if (!require(tseries)) {install.packages("tseries")}
if (!require(Hmisc)) {install.packages("Hmisc")}
if (!require(knitr)) {install.packages("knitr")}
if (!require(e1071)) {install.packages("e1071")}
if (!require(bain)) {install.packages("bain")}
if (!require(jpeg)) {install.packages("jpeg")}

load(file = "Workspaces\\BFSelection.RData")

#############################################################################################
#Functions
#############################################################################################

tr_plot <- function(matrices) {
  bparams <- c("b0", "b1", "b2", "b3", "var")
  plot_list <- list()
  #This function produces trace plots, one of the four ways in which convergence can be assessed. For each parameter of interest, the values of the two chains are overlayed to see if they create the desired 'fat caterpillars'. 
  
  par(mfrow = c((ncol(matrices[[1]])-1),1)) #Create room so all plots can be placed underneath each other
  label_vals <- colnames(matrices[[1]]) #derive the column names from the input matrix which can be used as labels for the plots
  
  for (i in 1:(ncol(matrices[[1]]))) {
    jpeg(paste("Plots\\", "tr_plot_", i, ".jpeg", sep = ""))
    plot(matrices[[1]][,i], type = 'l', col = 'red', ylab = bparams[i]) #Plot the values of the current regression coefficient for chain 1
    lines(matrices[[2]][,i]) #Overlay the values of the current regression coefficient for chain 2
    dev.off()
  }
}

acf_calculation <- function(full_matrix, n) {
  bparams <- c("b0", "b1", "b2", "b3", "var")
  for (coef in 1:ncol(full_matrix)) {
    jpeg(paste("Plots\\", "acf_plot_", coef, ".jpeg", sep = ""))
    y <- full_matrix[,coef] #retrieve the data of a coefficient
    correlations <- rep(0,n+1) #create space for the autocorrelations
    for (i in 1:(n+1)) {
      correlations[i] <- cor(y, Lag(y, shift = i-1), use = "complete.obs") #calculate the autocorrelation at lag i
    } 
    plot(0, xlim = c(0,n), ylim = c(-1,1), type = "n", ylab = bparams[coef]) #create the plot template for the autocorrelation
    abline(h = 0) #create a horizontal line at y = 0
    lags <- seq(0,n+1,1) #create labels for the different lags
    for (lag in 1:(length(lags)+1)) {
      segments(x0 = lag, y0 = 0, x1 = lag, y1 = correlations[lag], col = "blue") #plot a vertical line at lag t where the high is equal to the value of the autocorrelation
    }
    dev.off()
  }
}

gelman_calculation <- function(samples, n.chains, L) {
  bparams <- c("b0", "b1", "b2", "b3", "var")
  for (j in 1:ncol(samples[[1]])) {
    jpeg(paste("Plots\\", "gr_plot_", j, ".jpeg", sep = ""))
    n.iter <- L/70 #determine the number of iterations for the Gelman-Rubin calculations
    Bs <- rep(0,n.iter) #create space for the between variances
    Ws <- rep(0,n.iter) #create space for the within variances
    GRs <- rep(0,n.iter) #create space for the Gelman-Rubin statistics
    indices <- rep(0,n.iter) #create space for the labels for the x-axis
    for (i in 1:n.iter) {
      x1 <- samples[[1]][1:(i*70),j] #extract values 1 to i*70 from first chain
      x2 <- samples[[2]][1:(i*70),j] #extract values 1 to i*70 from second chain
      x1_mean <- mean(x1) #calculate mean of the values from chain 1
      x2_mean <- mean(x2) #calculate mean of the values from chain 2
      grand_mean <- mean(c(x1_mean, x2_mean)) #calculate grand mean
      
      x1_var <- sum((x1_mean - grand_mean)^2) #variance in chain 1
      x2_var <- sum((x2_mean - grand_mean)^2) #variance in chain 2
      B <- L*(x1_var+x2_var) #between chain variance
      Bs[i] <- B #store the current between variance
      
      x1_samplevar <- (1/L) * sum((x1 - x1_mean)^2) #sample variance chain 1
      x2_samplevar <- (1/L) * sum((x2 - x2_mean)^2) #sample variance chain 2
      W <- (1/n.chains) * (x1_samplevar+x2_samplevar) #within variance
      Ws[i] <- W #store the current within variance
      
      GR <- ((((L-1)/L)*W) + ((1/L)*B)) / W #Gelman-Rubin Statistic
      GRs[i] <- GR #store the current Gelman-Rubin statistic.
      indices[i] <- (i*70) #store the current number of used iterations to be used as label for the x-axis
    }
    plot(GRs, type = "l", xaxt = "n", ylim = c(0,1.2), xlab = "Iteration", ylab = bparams[j]) #plot the Gelman-Rubin statistics in a line plot
    axis(1, at = 1:n.iter, labels = indices) #create an x-axis that represents the number of iterations
    lines(Bs, col = "blue") #add the between variances to the Gelman-Rubin plot as a blue line
    lines(Ws, col = "red") #add the within variances to the Gelman-Rubin plot as a red line
    dev.off()
  }
}

#############################################################################################
#Evaluating convergence
#Run these lines to produce the plots described in section 6: Convergence
#############################################################################################

tr_plot(samples4) #plot the traceplots of the regression coefficients
autocorrelations <- acf_calculation(full_matrix3, 50) #plot the autocorrelation plots of the regression coefficients
gelman_calculation(samples4, length(samples4), nrow(samples4[[1]])) #plot the Gelman-Rubin plots of the regression coefficients

#############################################################################################
#MC Error
#Run these lines to produce the MC errors described in section 6: Convergence
#############################################################################################

mc_output <- matrix(0,5,3) #create an output matrix
mc_output[,1] <- apply(full_matrix3, 2, sd) #calculate the standard deviation of the regression coefficients
mc_output[,2] <- apply(full_matrix3, 2, function(x) {sd(x)/sqrt(n.iters)}) #calculate the MC error of the regression coefficients
mc_output[,3] <- apply(full_matrix3, 2, function(x) {(sd(x)/sqrt(n.iters))/sd(x)}) #calculate the ratio MC error / standard deviation
rownames(mc_output) <- c("b0", "b1", "b2", "b3", "var")
kable(mc_output, col.names = c("SD","MC error", "MC error/SD")) #print the MC error results

save.image(file = "Workspaces\\Convergence.RData")
