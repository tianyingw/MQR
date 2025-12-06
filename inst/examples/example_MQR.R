#' ============================================================================
#' MQR Package Example
#' Joint Time-varying Quantile Regressions with Multi-outcome Latent Groups
#' ============================================================================
#'
#' This example demonstrates the basic usage of the MQR package for:
#' 1. Simulating data from the MQR model
#' 2. Fitting the MQR model
#' 3. Evaluating results
#'
#' Reference:
#' Wang, T., Ma, Y., and Wei, Y. (2024). Joint Time-varying Quantile Regressions
#' with Multi-outcome Latent Groups. Journal of the Royal Statistical Society.

# Load the package
library(MQR)

# Set seed for reproducibility
set.seed(2024)

# ============================================================================
# Step 1: Simulate Data
# ============================================================================

cat("=== Simulating Data ===\n\n")

# Generate simulated data with:
# - n = 200 subjects
# - K = 10 traits (phenotypes)
# - M = 3 latent groups
# - J = 1 time point per subject
# - p = 1 primary covariate (with time-varying effect)
# - q = 2 additional covariates (with constant effects)

sim_data <- simulate_MQR_data(
  n = 200,       # Number of subjects
  K = 10,        # Number of traits
  M = 3,         # Number of latent groups
  J = 1,         # Time points per subject
  p = 1,         # Primary covariates
  q = 2,         # Additional covariates
  sigma = 1,     # Error SD
  seed = 123
)

# Display data dimensions
cat("Data dimensions:\n")
cat(sprintf("  Y (outcomes): %d x %d\n", nrow(sim_data$Y), ncol(sim_data$Y)))
cat(sprintf("  X (primary covariates): %d x %d\n", nrow(sim_data$X), ncol(sim_data$X)))
cat(sprintf("  C (additional covariates): %d x %d\n", nrow(sim_data$C), ncol(sim_data$C)))
cat(sprintf("  T_vec (time points): length %d\n\n", length(sim_data$T_vec)))

# Display true group membership
cat("True group membership:\n")
print(table(Group = sim_data$true_gamma))
cat("\n")

# ============================================================================
# Step 2: Fit MQR Model
# ============================================================================

cat("=== Fitting MQR Model ===\n\n")

# Fit MQR at the median (tau = 0.5) with known M = 3
fit <- MQR(
  Y = sim_data$Y,
  X = sim_data$X,
  C = sim_data$C,
  T_vec = sim_data$T_vec,
  tau = 0.5,           # Median regression
  M = 3,               # Number of groups (known)
  df = 4,              # B-spline degrees of freedom
  max_iter = 100,      # Maximum iterations
  tol = 1e-6,          # Convergence tolerance
  standardize = TRUE,  # Standardize outcomes by rank
  verbose = TRUE       # Print progress
)

# ============================================================================
# Step 3: View Results
# ============================================================================

cat("\n=== Model Results ===\n\n")

# Print basic summary
print(fit)

# ============================================================================
# Step 4: Evaluate Classification Accuracy
# ============================================================================

cat("\n=== Classification Evaluation ===\n\n")

# Compare estimated vs true group membership
confusion <- table(Estimated = fit$gamma, True = sim_data$true_gamma)
cat("Confusion matrix:\n")
print(confusion)

# Calculate classification accuracy
# Note: Group labels may be permuted, so we find the best matching
accuracy <- sum(diag(confusion)) / sum(confusion)
cat(sprintf("\nRaw accuracy: %.2f%%\n", accuracy * 100))

# Find best permutation match
find_best_accuracy <- function(est, true, M) {
  perms <- combinat::permn(1:M)
  if (!requireNamespace("combinat", quietly = TRUE)) {
    # Simple matching if combinat not available
    return(mean(est == true))
  }

  best_acc <- 0
  for (perm in perms) {
    est_relabeled <- perm[est]
    acc <- mean(est_relabeled == true)
    if (acc > best_acc) best_acc <- acc
  }
  return(best_acc)
}

# Simple heuristic for matching (works for M=3)
match_groups <- function(est, true, M) {
  conf <- table(est, true)
  est_new <- est
  for (m in 1:M) {
    best_match <- which.max(conf[m, ])
    est_new[est == m] <- best_match
  }
  mean(est_new == true)
}

matched_accuracy <- match_groups(fit$gamma, sim_data$true_gamma, fit$M)
cat(sprintf("Matched accuracy: %.2f%%\n", matched_accuracy * 100))

# ============================================================================
# Step 5: Predict Time-varying Coefficients
# ============================================================================

cat("\n=== Predicted Time-varying Coefficients ===\n\n")

# Predict alpha(t) at a grid of time points
time_grid <- seq(0.05, 0.95, by = 0.05)
pred <- predict(fit, newT = time_grid)

cat("Predicted alpha(t) at selected time points:\n")
cat("Time\t", paste(sprintf("Group%d", 1:fit$M), collapse = "\t"), "\n")
for (i in seq(1, length(time_grid), by = 4)) {
  cat(sprintf("%.2f\t", time_grid[i]))
  cat(paste(sprintf("%.4f", pred$alpha_T[i, ]), collapse = "\t"), "\n")
}

# ============================================================================
# Step 6: Visualization (if graphics available)
# ============================================================================

cat("\n=== Visualization ===\n\n")

# Check if we can create plots
if (interactive() || !is.null(dev.list())) {

  # Plot 1: Time-varying coefficients
  par(mfrow = c(1, 2))

  # Estimated coefficients
  matplot(time_grid, pred$alpha_T, type = "l", lty = 1, lwd = 2,
          col = 1:fit$M, xlab = "Time (t)", ylab = expression(alpha(t)),
          main = "Estimated Time-varying Coefficients")
  legend("topright", paste("Group", 1:fit$M), col = 1:fit$M,
         lty = 1, lwd = 2, bty = "n")

  # True vs Estimated group membership
  plot(1:sim_data$K, sim_data$true_gamma, pch = 19, col = "blue",
       ylim = c(0.5, fit$M + 0.5), xlab = "Trait", ylab = "Group",
       main = "Group Membership")
  points(1:sim_data$K, fit$gamma, pch = 17, col = "red")
  legend("topright", c("True", "Estimated"),
         pch = c(19, 17), col = c("blue", "red"), bty = "n")

  par(mfrow = c(1, 1))

  cat("Plots generated successfully.\n")
} else {
  cat("Skipping plots (non-interactive mode).\n")
}

# ============================================================================
# Step 7: Full Summary
# ============================================================================

cat("\n=== Full Model Summary ===\n\n")
summary(fit)

cat("\n=== Example Complete ===\n")
