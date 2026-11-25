Sys.setenv(JAVA_HOME = "C:/Programmi/Eclipse Adoptium/jdk-25.0.1.8-hotspot")
library(rJava)
.jinit()

options(java.parameters = "-Xmx5000m")
library(bartMachine)
library(MASS)
set_bart_machine_num_cores(3)
library(mvtnorm)

generate_response_model = function(n, Sigma_e = sigma_e, Mu_b = mu_b, Sigma_b = sigma_b, Sigma = sigma
                                   , mu_vec = mu, Gamma = gamma,Lambda = lambda, Min = min, Max = max){
  
  Xs = rmvnorm(n, mu_vec, Sigma)
  b = rnorm(n, Mu_b, Sigma_b)
  s = runif(n,min = Min, max = Max)
  error = rnorm(n, 0, Sigma_e)
  X1 = Xs[, 1]
  X2 = Xs[, 2]
  X3 = Xs[, 3]
  X4 <- rpois(n,Lambda)
  y = X1 + X2 + 2 * X3 - X1^2 + X2^2 + X1 * X2 + X3*X4 + sqrt(X4) + error
  
  
  for (i in 1 : n){
    if (X2[i] > 0 & runif(1) < Gamma){
      X2[i] = NA
    }		
  }
  
  for (i in 1 : n){
    if (X1[i] > 0 & runif(1) < Gamma){
      X3[i] = NA
    }		
  }	
  
  for (i in 1 : n){
    if (runif(1) < 0.15){
      X4[i] = NA
    }		
  }
  b = ifelse(is.na(X3), b, 0)
  s = ifelse(is.na(X4), s, 0)
  y = y + b + s
  
  return(data.frame(X1 = X1, X2 = X2, X3 = X3, X4 = X4, y = y))
}

plot_hist_of_posterior = function(pred, expectation, cis){
  hist(pred$y_hat_posterior_samples[1,], 
       br = 50,
       main = "",
       xlab = ""
  )
  cat("yhat =", round(mean(pred$y_hat_posterior_samples[1,]), 2), "se(yhat) =", round(sd(pred$y_hat_posterior_samples[1,]), 3), "\n")
  abline(v = expectation, col = "blue", lwd = 3)
  abline(v = pred$y_hat[1], col = "green", lwd = 3)
  
  abline(v = cis[1], col = "orange", lwd = 3)
  abline(v = cis[2], col = "orange", lwd = 3)
}

n = 2000
sigma <- matrix(c(1,0.2,0.4,0.2,1.5,0.1,0.4,0.1,1),nrow=3)
mu <- c(0,1,2)
lambda <- 1
min <- -20
max <- -10
mu_b <- 10
sigma_b <- 1
sigma_e <- 1
gamma <- 0.3

set.seed(110902)
Xy = generate_response_model(n)

bart_machine = build_bart_machine(Xy = Xy, use_missing_data = TRUE, use_missing_data_dummies_as_covars = TRUE, num_burn_in = 1000, run_in_sample = FALSE)

#(0,1,NA,NA)
new_data = as.data.frame(t(as.matrix(c(0, 1, NA, NA))))
colnames(new_data) = colnames(Xy[, 1 : 4])
pred = bart_machine_get_posterior(bart_machine, new_data)
expe = 3.77
plot_hist_of_posterior(pred, expe, calc_credible_intervals(bart_machine, new_data))

#(0.5,NA,NA,2)
new_data = as.data.frame(t(as.matrix(c(0.5, NA, NA, 2))))
colnames(new_data) = colnames(Xy[, 1 : 4])
pred = bart_machine_get_posterior(bart_machine, new_data)
expe = 25.37
plot_hist_of_posterior(pred, expe, calc_credible_intervals(bart_machine, new_data))

#(0,NA,1.8,NA)
new_data = as.data.frame(t(as.matrix(c(0, NA, 1.8, NA))))
colnames(new_data) = colnames(Xy[, 1 : 4])
pred = bart_machine_get_posterior(bart_machine, new_data)
expe = -4.51
plot_hist_of_posterior(pred, expe, calc_credible_intervals(bart_machine, new_data))

#(0.2,NA,NA,NA)
new_data = as.data.frame(t(as.matrix(c(0.2, NA, NA, NA))))
colnames(new_data) = colnames(Xy[, 1 : 4])
pred = bart_machine_get_posterior(bart_machine, new_data)
expe = 6.89
plot_hist_of_posterior(pred, expe, calc_credible_intervals(bart_machine, new_data))
