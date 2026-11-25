Sys.setenv(JAVA_HOME = "C:/Programmi/Eclipse Adoptium/jdk-25.0.1.8-hotspot")
library(rJava)
.jinit()
options(java.parameters = "-Xmx5000m")
library(bartMachine)
library(MASS)
library(mvtnorm)
set_bart_machine_num_cores(3)
# Funzione di generazione dati (la tua)
generate_response_model = function(n, Sigma_e = sigma_e, Mu_b = mu_b,
Sigma_b = sigma_b, Sigma = sigma,
mu_vec = mu, Gamma = gamma, Lambda = lambda,
Min = min, Max = max){
Xs = rmvnorm(n, mu_vec, Sigma)
b = rnorm(n, Mu_b, Sigma_b)
s = runif(n, min = Min, max = Max)
error = rnorm(n, 0, Sigma_e)
X1 = Xs[, 1]
X2 = Xs[, 2]
X3 = Xs[, 3]
X4 <- rpois(n, Lambda)
y = X1 + X2 + 2 * X3 - X1^2 + X2^2 + X1 * X2 + X3*X4 + sqrt(X4) + error
for (i in 1:n){
if (X2[i] > 0 & runif(1) < Gamma){
X2[i] = NA
}
}
for (i in 1:n){
if (X1[i] > 0 & runif(1) < Gamma){
X3[i] = NA
}
}
for (i in 1:n){
if (runif(1) < 0.15){
X4[i] = NA
}
}
b = ifelse(is.na(X3), b, 0)
s = ifelse(is.na(X4), s, 0)
y = y + b + s
return(data.frame(X1 = X1, X2 = X2, X3 = X3, X4 = X4, y = y))
}
# Parametri del modello
n = 2000
sigma <- matrix(c(1, 0.2, 0.4, 0.2, 1.5, 0.1, 0.4, 0.1, 1), nrow = 3)
mu <- c(0, 1, 2)
lambda <- 1
min <- -20
max <- -10
mu_b <- 10
sigma_b <- 1
sigma_e <- 1
gamma <- 0.3
set.seed(110902)
Xy = generate_response_model(n)
bart_machine = build_bart_machine(
Xy = Xy,
use_missing_data = TRUE,
use_missing_data_dummies_as_covars = TRUE,
num_burn_in = 1000,        # Burn-in iterations
num_iterations_after_burn_in = 1000,  # Post burn-in samples
run_in_sample = FALSE,
verbose = FALSE
)
plot_bart_diagnostics <- function(bart_machine) {
# Prepara layout per multipli grafici
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
# -------------------------------------------------------------------------
# 1. SIGMA² NEL TEMPO (LA DIAGNOSTICA PIÙ IMPORTANTE)
# -------------------------------------------------------------------------
sigsqs <- bart_machine$sigsqs_after_burn_in
burn_in <- bart_machine$num_burn_in
total_iter <- length(sigsqs) + burn_in
# Include anche il burn-in se disponibile
if (!is.null(bart_machine$sigsqs_during_burn_in)) {
all_sigsqs <- c(bart_machine$sigsqs_during_burn_in, sigsqs)
} else {
all_sigsqs <- sigsqs
}
plot(1:length(all_sigsqs), all_sigsqs,
type = "l",
col = "darkblue",
xlab = "Iterazione MCMC",
ylab = expression(sigma^2),
main = expression(paste("Trace Plot di ", sigma^2)),
lwd = 1.5)
# Linea verticale per il burn-in
if (burn_in > 0 && length(all_sigsqs) > burn_in) {
abline(v = burn_in, col = "red", lwd = 2, lty = 2)
text(burn_in, max(all_sigsqs) * 0.9, "Fine Burn-in",
pos = 4, col = "red", cex = 0.8)
}
# Media post burn-in
mean_sigsq <- mean(sigsqs)
abline(h = mean_sigsq, col = "darkgreen", lwd = 2, lty = 2)
text(length(all_sigsqs) * 0.7, mean_sigsq * 1.1,
paste0("Media = ", round(mean_sigsq, 3)),
col = "darkgreen", cex = 0.8)
# -------------------------------------------------------------------------
# 2. ISTOGRAMMA DELLA DISTRIBUZIONE DI SIGMA² POST BURN-IN
# -------------------------------------------------------------------------
hist(sigsqs,
breaks = 30,
col = "lightblue",
border = "darkblue",
main = expression(paste("Distribuzione Posteriore di ", sigma^2)),
xlab = expression(sigma^2),
ylab = "Frequenza")
abline(v = mean_sigsq, col = "red", lwd = 2, lty = 2)
abline(v = quantile(sigsqs, c(0.025, 0.975)),
col = "orange", lwd = 2, lty = 2)
legend("topright",
legend = c("Media", "95% CI"),
col = c("red", "orange"),
lty = 2, lwd = 2, cex = 0.7)
# -------------------------------------------------------------------------
# 3. AUTOCORRELATION FUNCTION (ACF) DI SIGMA²
# -------------------------------------------------------------------------
acf(sigsqs,
main = expression(paste("ACF di ", sigma^2, " (Post Burn-in)")),
col = "darkblue",
lwd = 2)
# -------------------------------------------------------------------------
# 4. RUNNING MEAN DI SIGMA² (CONVERGENZA)
# -------------------------------------------------------------------------
running_mean <- cumsum(sigsqs) / seq_along(sigsqs)
plot(1:length(running_mean), running_mean,
type = "l",
col = "purple",
lwd = 2,
xlab = "Iterazione (post burn-in)",
ylab = expression(paste("Running Mean di ", sigma^2)),
main = "Convergenza della Media")
abline(h = mean_sigsq, col = "red", lwd = 2, lty = 2)
# Band di stabilità (±5%)
abline(h = mean_sigsq * c(0.95, 1.05),
col = "orange", lwd = 1, lty = 3)
par(mfrow = c(1, 1))
# -------------------------------------------------------------------------
# STAMPA STATISTICHE RIASSUNTIVE
# -------------------------------------------------------------------------
cat("\n========================================\n")
cat("DIAGNOSTICHE MCMC - RIASSUNTO\n")
cat("========================================\n\n")
cat("Iterazioni burn-in:", burn_in, "\n")
cat("Iterazioni post burn-in:", length(sigsqs), "\n\n")
cat("SIGMA² - Statistiche Posteriori:\n")
cat("  Media:", round(mean_sigsq, 4), "\n")
cat("  Mediana:", round(median(sigsqs), 4), "\n")
cat("  SD:", round(sd(sigsqs), 4), "\n")
cat("  95% CI: [", round(quantile(sigsqs, 0.025), 4), ",",
round(quantile(sigsqs, 0.975), 4), "]\n\n")
# Effective Sample Size (approssimato)
acf_vals <- acf(sigsqs, plot = FALSE)$acf
# Trova il primo lag dove ACF < 0.05
first_low_acf <- which(abs(acf_vals) < 0.05)[1]
if (!is.na(first_low_acf)) {
ess_approx <- length(sigsqs) / first_low_acf
cat("Effective Sample Size (approssimato):", round(ess_approx, 0), "\n")
}
cat("\n========================================\n\n")
}
plot_bart_diagnostics(bart_machine)
# Esegui trace plots per alcune predizioni
# plot_prediction_traces(bart_machine, n_obs = 4)   # funzione non definita in questo script
plot_tree_diagnostics <- function(bart_machine) {
cat("\n=== Diagnostiche Struttura Alberi ===\n\n")
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
# Numero medio di nodi terminali per albero nel tempo
# Nota: bartMachine non espone direttamente questa info,
# ma possiamo ottenere info aggregate
cat("Numero di alberi (m):", bart_machine$num_trees, "\n")
cat("Numero di iterazioni post burn-in:",
bart_machine$num_iterations_after_burn_in, "\n\n")
# Plot della variabile inclusion proportions
var_props <- get_var_props_over_chain(bart_machine)
barplot(sort(colMeans(var_props), decreasing = TRUE),
las = 2,
col = "steelblue",
main = "Variable Inclusion Proportions",
ylab = "Proporzione",
cex.names = 0.8)
# Heatmap delle inclusion proportions nel tempo
image(t(var_props[, order(colMeans(var_props), decreasing = TRUE)]),
main = "Inclusion Proportions nel Tempo",
xlab = "Iterazione MCMC",
ylab = "Variabile (ordinate per importanza)",
col = heat.colors(100))
par(mfrow = c(1, 1))
}
plot_tree_diagnostics(bart_machine)
cat("\n\n=== CONFRONTO: Effetto di Burn-in Insufficiente ===\n\n")
# Modello con burn-in molto corto
bart_machine_short_burnin = build_bart_machine(
Xy = Xy,
use_missing_data = TRUE,
use_missing_data_dummies_as_covars = TRUE,
num_burn_in = 50,  # TROPPO CORTO!
num_iterations_after_burn_in = 1000,
run_in_sample = FALSE,
verbose = FALSE
)
# Confronto visivo
par(mfrow = c(1, 2))
# Burn-in adeguato
plot(bart_machine$sigsqs_after_burn_in,
type = "l",
col = "darkblue",
lwd = 1.5,
main = "Burn-in = 1000 (Adeguato)",
xlab = "Iterazione",
ylab = expression(sigma^2),
ylim = range(c(bart_machine$sigsqs_after_burn_in,
bart_machine_short_burnin$sigsqs_after_burn_in)))
abline(h = mean(bart_machine$sigsqs_after_burn_in),
col = "red", lwd = 2, lty = 2)
# Burn-in insufficiente
plot(bart_machine_short_burnin$sigsqs_after_burn_in,
type = "l",
col = "purple",
lwd = 1.5,
main = "Burn-in = 50 (Insufficiente!)",
xlab = "Iterazione",
ylab = expression(sigma^2),
ylim = range(c(bart_machine$sigsqs_after_burn_in,
bart_machine_short_burnin$sigsqs_after_burn_in)))
abline(h = mean(bart_machine_short_burnin$sigsqs_after_burn_in),
col = "red", lwd = 2, lty = 2)
par(mfrow = c(1, 1))
cat("\nDifferenza nelle medie di sigma²:\n")
cat("  Burn-in adeguato:",
round(mean(bart_machine$sigsqs_after_burn_in), 4), "\n")
cat("  Burn-in breve:",
round(mean(bart_machine_short_burnin$sigsqs_after_burn_in), 4), "\n")
cat("  Differenza:",
round(abs(mean(bart_machine$sigsqs_after_burn_in) -
mean(bart_machine_short_burnin$sigsqs_after_burn_in)), 4), "\n\n")
plot_hist_of_posterior = function(pred, expectation, cis, title = ""){
hist(pred$y_hat_posterior_samples[1,],
breaks = 50,
main = title,
xlab = "Predizione",
col = "lightblue",
border = "darkblue"
)
cat("yhat =", round(mean(pred$y_hat_posterior_samples[1,]), 2),
"se(yhat) =", round(sd(pred$y_hat_posterior_samples[1,]), 3), "\n")
abline(v = expectation, col = "blue", lwd = 3)
abline(v = pred$y_hat[1], col = "green", lwd = 3)
abline(v = cis[1], col = "orange", lwd = 3)
abline(v = cis[2], col = "orange", lwd = 3)
legend("topright",
legend = c("Valore Atteso", "Media Post.", "95% CI"),
col = c("blue", "green", "orange"),
lwd = 3,
cex = 0.8)
}
# Test predizioni con missing data
cat("\n\n=== PREDIZIONI CON MISSING DATA ===\n\n")
# Scenario 1: (0, 1, NA, NA)
new_data = as.data.frame(t(as.matrix(c(0, 1, NA, NA))))
colnames(new_data) = colnames(Xy[, 1:4])
pred = bart_machine_get_posterior(bart_machine, new_data)
cat("\nScenario 1: X = (0, 1, NA, NA)\n")
plot_hist_of_posterior(pred, 3.77,
calc_credible_intervals(bart_machine, new_data),
title = "X = (0, 1, NA, NA)")
# Scenario 2: (0.5, NA, NA, 2)
new_data = as.data.frame(t(as.matrix(c(0.5, NA, NA, 2))))
colnames(new_data) = colnames(Xy[, 1:4])
pred = bart_machine_get_posterior(bart_machine, new_data)
cat("\nScenario 2: X = (0.5, NA, NA, 2)\n")
plot_hist_of_posterior(pred, 25.37,
calc_credible_intervals(bart_machine, new_data),
title = "X = (0.5, NA, NA, 2)")
# Scenario 3: (0, NA, 1.8, NA)
new_data = as.data.frame(t(as.matrix(c(0, NA, 1.8, NA))))
colnames(new_data) = colnames(Xy[, 1:4])
pred = bart_machine_get_posterior(bart_machine, new_data)
cat("\nScenario 3: X = (0, NA, 1.8, NA)\n")
plot_hist_of_posterior(pred, -4.51,
calc_credible_intervals(bart_machine, new_data),
title = "X = (0, NA, 1.8, NA)")
# Scenario 4: (0.2, NA, NA, NA)
new_data = as.data.frame(t(as.matrix(c(0.2, NA, NA, NA))))
colnames(new_data) = colnames(Xy[, 1:4])
pred = bart_machine_get_posterior(bart_machine, new_data)
cat("\nScenario 4: X = (0.2, NA, NA, NA)\n")
plot_hist_of_posterior(pred, 6.89,
calc_credible_intervals(bart_machine, new_data),
title = "X = (0.2, NA, NA, NA)")
