Sys.setenv(JAVA_HOME = "C:/Programmi/Eclipse Adoptium/jdk-25.0.1.8-hotspot")
library(rJava)
.jinit()

options(java.parameters = "-Xmx5000m")
library(bartMachine)
library(MASS)
set_bart_machine_num_cores(3)
library(mvtnorm)

n = 500
sigma <- matrix(c(1,0.2,0.4,0.2,1.5,0.1,0.4,0.1,1),nrow=3)
mu <- c(0,1,2)
lambda <- 1
sigma_e <- 1
generate_response_model = function(n, sigma_e = 1, Sigma = sigma, mu_vec = mu, Lambda = lambda){

  Xs = mvrnorm(n, mu_vec, Sigma)
  error = rnorm(n, 0, sigma_e)
  X1 = Xs[, 1]
  X2 = Xs[, 2]
  X3 = Xs[, 3]
  X4 <- rpois(n,Lambda)
  
  y_crazy = X1 + X2 + 2 * X3 - X1^2 + X2^2 + X1 * X2 + X3*X4 + sqrt(X4) + error
  
  data.frame(Xs,X4, y_crazy)
}

Nsim = 10
n = 150

###MCAR


generate_mcar_model = function(Xy, gamma){
  for (i in 1 : nrow(Xy)){
    for (j in 1 : (ncol(Xy) - 1)){
      if (runif(1) < gamma){
        Xy[i, j] = NA
      }
    }
  }
  Xy
}

approx_prop_missing = seq(from = 0, to = 0.8, by = 0.2)
gammas = 1 - (1 - approx_prop_missing)^(1/4)


results_bart_all_all_mcar = matrix(NA, nrow = length(approx_prop_missing), ncol = Nsim)
results_bart_all_cc_mcar = matrix(NA, nrow = length(approx_prop_missing), ncol = Nsim)
results_bart_cc_all_mcar = matrix(NA, nrow = length(approx_prop_missing), ncol = Nsim)
results_bart_cc_cc_mcar = matrix(NA, nrow = length(approx_prop_missing), ncol = Nsim)
rownames(results_bart_all_all_mcar) = approx_prop_missing
rownames(results_bart_all_cc_mcar) = approx_prop_missing
rownames(results_bart_cc_all_mcar) = approx_prop_missing
rownames(results_bart_cc_cc_mcar) = approx_prop_missing

for (nsim in 1 : Nsim){
  cat("nsim = ", nsim, "\n")
  
  for (g in 1 : length(gammas)){
    gamma = gammas[g]
    cat("g = ", g, "\n")
    
    Xy_train_all = generate_mcar_model(generate_response_model(n), gamma)
    Xy_train_cc = na.omit(Xy_train_all)
    Xy_test_cc = generate_response_model(n)
    Xy_test_all = generate_mcar_model(Xy_test_cc, gamma)
    
    #train models
    bart_machine_all = build_bart_machine(Xy = Xy_train_all, use_missing_data = TRUE, use_missing_data_dummies_as_covars = TRUE, verbose = FALSE, run_in_sample = FALSE)
    
    #jet if we have exceedingly few rows in Xycc
    if (nrow(Xy_train_cc) > 5){
      bart_machine_cc = build_bart_machine(Xy = Xy_train_cc, use_missing_data = TRUE, verbose = FALSE, run_in_sample = FALSE)
    }
    #test models
    results_bart_all_all_mcar[g, nsim] = bart_predict_for_test_data(bart_machine_all, Xtest = Xy_test_all[, 1 : 4], ytest = Xy_test_all[, 5])$rmse
    results_bart_all_cc_mcar[g, nsim] = bart_predict_for_test_data(bart_machine_all, Xtest = Xy_test_cc[, 1 : 4], ytest = Xy_test_cc[, 5])$rmse
    if (nrow(Xy_train_cc) > 5){ #jet if we have exceedingly few rows in Xycc
      results_bart_cc_all_mcar[g, nsim] = bart_predict_for_test_data(bart_machine_cc, Xtest = Xy_test_all[, 1 : 4], ytest = Xy_test_all[, 5])$rmse
      results_bart_cc_cc_mcar[g, nsim] = bart_predict_for_test_data(bart_machine_cc, Xtest = Xy_test_cc[, 1 : 4], ytest = Xy_test_cc[, 5])$rmse
    }
  }	
  #fine ciclo gamma, abbiamo calcolato rmse_s per i due modelli e per i due test per ogni valore di gamma (iterazione nsim)
  
  avgs_mcar_all_all = apply(results_bart_all_all_mcar, 1, mean, na.rm = TRUE)	
  rel_mcar_avgs_all_all = avgs_mcar_all_all / avgs_mcar_all_all[1]	
  sd_mcar_all_all = apply(results_bart_all_all_mcar / avgs_mcar_all_all[1], 1, sd, na.rm = TRUE)
  
  avgs_mcar_all_cc = apply(results_bart_all_cc_mcar, 1, mean, na.rm = TRUE)	
  rel_mcar_avgs_all_cc = avgs_mcar_all_cc / avgs_mcar_all_all[1]
  sd_mcar_all_cc = apply(results_bart_all_cc_mcar / avgs_mcar_all_all[1], 1, sd, na.rm = TRUE)
  
  avgs_mcar_cc_all = apply(results_bart_cc_all_mcar, 1, mean, na.rm = TRUE)
  rel_mcar_avgs_cc_all = avgs_mcar_cc_all / avgs_mcar_all_all[1]
  sd_mcar_cc_all = apply(results_bart_cc_all_mcar / avgs_mcar_all_all[1], 1, sd, na.rm = TRUE)
  
  avgs_mcar_cc_cc = apply(results_bart_cc_cc_mcar, 1, mean, na.rm = TRUE)
  rel_mcar_avgs_cc_cc = avgs_mcar_cc_cc / avgs_mcar_all_all[1]
  sd_mcar_cc_cc = apply(results_bart_cc_cc_mcar / avgs_mcar_all_all[1], 1, sd, na.rm = TRUE)
  
  #Figure 2a
  par(mar = c(4.2,4,0.2,0.2))
  plot(approx_prop_missing, 
       rel_mcar_avgs_all_all, 
       col = "blue", 
       type = "o", 
       ylim = c(1, max(rel_mcar_avgs_all_all, rel_mcar_avgs_all_cc, rel_mcar_avgs_cc_all, rel_mcar_avgs_cc_cc, na.rm = TRUE)),
       xlab = "Proportion Missing",
       ylab = "Multiple of Baseline Error")

  points(approx_prop_missing, rel_mcar_avgs_all_cc, col = "blue", type = "o", lty = 3)

  points(approx_prop_missing, rel_mcar_avgs_cc_all, col = "red", type = "o")

  points(approx_prop_missing, rel_mcar_avgs_cc_cc, col = "red", type = "o", lty = 3)

  
  
  save.image("sec_4.2_mcar.RData")
}

