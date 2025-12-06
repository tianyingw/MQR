#' @title Check Loss Function (Rho)
#'
#' @description Compute the check loss function used in quantile regression.
#'
#' @param u Residual vector
#' @param tau Quantile level
#'
#' @return Vector of check loss values
#' @keywords internal
rho_tau <- function(u, tau) {
  abs(u) * ((1 - tau) * (u < 0) + tau * (u >= 0))
}


#' @title Compute Total Loss
#'
#' @description Compute the total check loss across all traits.
#'
#' @param Y Outcome matrix
#' @param X_B Covariate-spline interaction matrix
#' @param C Additional covariates matrix
#' @param alpha Group-specific coefficients
#' @param xi Trait-specific coefficients
#' @param intercept Trait-specific intercepts
#' @param gamma Group membership
#' @param tau Quantile level
#'
#' @return Total loss value
#' @keywords internal
compute_loss <- function(Y, X_B, C, alpha, xi, intercept, gamma, tau) {

  K <- ncol(Y)
  total_loss <- 0

  for (k in 1:K) {
    m <- gamma[k]
    residual <- Y[, k] - intercept[k] - X_B %*% alpha[m, ]

    if (!is.null(C) && ncol(xi) > 0) {
      residual <- residual - C %*% xi[k, ]
    }

    total_loss <- total_loss + sum(rho_tau(residual, tau))
  }

  return(total_loss)
}


#' @title Compute Information Criterion
#'
#' @description Compute the BIC-type information criterion for model selection.
#'
#' @param loss Total loss value
#' @param n Number of observations
#' @param K Number of traits
#' @param M Number of groups
#' @param df B-spline degrees of freedom
#' @param q Number of additional covariates
#'
#' @return IC value
#' @keywords internal
compute_IC <- function(loss, n, K, M, df, q) {
  if (is.null(q)) q <- 0
  n_params <- M * df + K * (q + 1)  # alpha params + xi params + intercepts
  IC <- loss / (n * K) + log(n * K) / (2 * n * K) * n_params
  return(IC)
}


#' @title Create X-B Spline Interaction
#'
#' @description Create interaction between covariates X and B-spline basis.
#'
#' @param X Covariate matrix
#' @param B B-spline basis matrix
#'
#' @return Matrix of X*B interactions
#' @keywords internal
create_XB_interaction <- function(X, B) {

  if (!is.matrix(X)) X <- as.matrix(X)

  n <- nrow(X)
  p <- ncol(X)
  df <- ncol(B)

  X_B <- matrix(0, nrow = n, ncol = p * df)

  for (j in 1:p) {
    idx <- ((j - 1) * df + 1):(j * df)
    X_B[, idx] <- X[, j] * B
  }

  return(X_B)
}


#' @title Standardize Outcomes by Marginal Rank
#'
#' @description Transform outcomes to marginal ranks using quantile regression.
#'
#' @param Y Outcome matrix (n x K)
#' @param T_vec Time vector
#' @param df Degrees of freedom for B-spline in quantile estimation
#'
#' @return Standardized outcome matrix
#' @keywords internal
standardize_by_rank <- function(Y, T_vec, df = 4) {

  K <- ncol(Y)
  n <- nrow(Y)
  tau_seq <- seq(0.01, 0.99, length.out = 99)

  Y_std <- matrix(NA, nrow = n, ncol = K)

  for (k in 1:K) {
    # Fit quantile regression across tau grid
    tryCatch({
      fit <- quantreg::rq(Y[, k] ~ splines::bs(T_vec, df = df, intercept = FALSE),
                          tau = tau_seq)
      fitted_vals <- fit$fitted.values

      # Find rank for each observation
      for (i in 1:n) {
        y_i <- Y[i, k]
        # Find smallest tau where y_i <= fitted quantile
        idx <- which(y_i <= fitted_vals[i, ])
        if (length(idx) == 0) {
          Y_std[i, k] <- max(tau_seq)
        } else {
          Y_std[i, k] <- tau_seq[min(idx)]
        }
      }
    }, error = function(e) {
      # If QR fails, use simple rank
      Y_std[, k] <<- rank(Y[, k]) / (n + 1)
    })
  }

  return(Y_std)
}


#' @title Initialize MQR by K-means
#'
#' @description Initialize parameters by fitting individual QR and clustering.
#'
#' @param Y Outcome matrix
#' @param X_B Covariate-spline interaction matrix
#' @param C Additional covariates matrix
#' @param tau Quantile level
#' @param M Number of groups
#' @param n_init Number of k-means initializations
#'
#' @return List with initial alpha, xi, intercept, gamma
#' @keywords internal
initialize_MQR <- function(Y, X_B, C, tau, M, n_init = 5) {

  K <- ncol(Y)
  p_alpha <- ncol(X_B)
  q <- if (is.null(C)) 0 else ncol(C)

  # Fit individual QR for each trait
  alpha_individual <- matrix(NA, nrow = K, ncol = p_alpha)
  xi_individual <- matrix(NA, nrow = K, ncol = max(q, 1))
  intercept_individual <- rep(NA, K)

  for (k in 1:K) {
    fit <- fit_single_trait_qr(Y[, k], X_B, C, tau)
    alpha_individual[k, ] <- fit$alpha
    xi_individual[k, ] <- fit$xi
    intercept_individual[k] <- fit$intercept
  }

  # Run k-means multiple times and select best
  best_loss <- Inf
  best_result <- NULL

  for (init in 1:n_init) {
    km <- stats::kmeans(alpha_individual, centers = M, nstart = 10)

    # Compute loss for this initialization
    alpha_init <- km$centers
    gamma_init <- km$cluster

    loss <- compute_loss(Y, X_B, C, alpha_init, xi_individual,
                         intercept_individual, gamma_init, tau)

    if (loss < best_loss) {
      best_loss <- loss
      best_result <- list(
        alpha = alpha_init,
        xi = xi_individual,
        intercept = intercept_individual,
        gamma = gamma_init
      )
    }
  }

  return(best_result)
}


#' @title Fit Single Trait Quantile Regression
#'
#' @description Fit QR for a single trait.
#'
#' @param y Outcome vector for one trait
#' @param X_B Covariate-spline interaction matrix
#' @param C Additional covariates matrix
#' @param tau Quantile level
#'
#' @return List with alpha, xi, intercept
#' @keywords internal
fit_single_trait_qr <- function(y, X_B, C, tau) {

  p_alpha <- ncol(X_B)

  if (is.null(C) || ncol(C) == 0) {
    # No additional covariates
    design <- cbind(1, X_B)
    fit <- quantreg::rq(y ~ X_B, tau = tau)

    alpha <- coef(fit)[2:(p_alpha + 1)]
    xi <- numeric(0)
    intercept <- coef(fit)[1]
  } else {
    q <- ncol(C)
    fit <- quantreg::rq(y ~ X_B + C, tau = tau)

    alpha <- coef(fit)[2:(p_alpha + 1)]
    xi <- coef(fit)[(p_alpha + 2):(p_alpha + 1 + q)]
    intercept <- coef(fit)[1]
  }

  return(list(alpha = alpha, xi = xi, intercept = intercept))
}


#' @title Fit Pooled Quantile Regression Within Group
#'
#' @description Fit QR by pooling data from multiple traits in the same group.
#'
#' @param Y Outcome matrix
#' @param X_B Covariate-spline interaction matrix
#' @param C Additional covariates matrix
#' @param g_idx Indices of traits in this group
#' @param tau Quantile level
#'
#' @return List with alpha, xi, intercept
#' @keywords internal
fit_pooled_qr <- function(Y, X_B, C, g_idx, tau) {

  n_obs <- nrow(Y)
  n_traits <- length(g_idx)
  p_alpha <- ncol(X_B)
  q <- if (is.null(C)) 0 else ncol(C)

  # Stack data for traits in this group
  Y_pooled <- as.vector(Y[, g_idx])
  X_B_pooled <- do.call(rbind, replicate(n_traits, X_B, simplify = FALSE))

  # Create trait-specific intercept indicators
  intercept_mat <- matrix(0, nrow = n_obs * n_traits, ncol = n_traits)
  for (j in 1:n_traits) {
    intercept_mat[((j - 1) * n_obs + 1):(j * n_obs), j] <- 1
  }

  # Create trait-specific covariate matrices
  if (!is.null(C) && q > 0) {
    C_pooled <- matrix(0, nrow = n_obs * n_traits, ncol = n_traits * q)
    for (j in 1:n_traits) {
      col_idx <- ((j - 1) * q + 1):(j * q)
      row_idx <- ((j - 1) * n_obs + 1):(j * n_obs)
      C_pooled[row_idx, col_idx] <- C
    }
    design <- cbind(X_B_pooled, C_pooled, intercept_mat)
  } else {
    design <- cbind(X_B_pooled, intercept_mat)
  }

  # Fit pooled QR (no intercept since we have trait-specific intercepts)
  fit <- quantreg::rq(Y_pooled ~ design - 1, tau = tau)
  coefs <- coef(fit)

  # Extract coefficients
  alpha <- coefs[1:p_alpha]

  if (!is.null(C) && q > 0) {
    xi_vec <- coefs[(p_alpha + 1):(p_alpha + n_traits * q)]
    xi <- matrix(xi_vec, nrow = n_traits, ncol = q, byrow = TRUE)
    intercept <- coefs[(p_alpha + n_traits * q + 1):(p_alpha + n_traits * q + n_traits)]
  } else {
    xi <- matrix(0, nrow = n_traits, ncol = 1)
    intercept <- coefs[(p_alpha + 1):(p_alpha + n_traits)]
  }

  return(list(alpha = alpha, xi = xi, intercept = intercept))
}


#' @title Select Number of Groups by Information Criterion
#'
#' @description Select optimal M using BIC-type information criterion.
#'
#' @param Y Outcome matrix
#' @param X_B Covariate-spline interaction matrix
#' @param C Additional covariates matrix
#' @param tau Quantile level
#' @param M_max Maximum number of groups to consider
#' @param df B-spline degrees of freedom
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param n_init Number of initializations
#' @param verbose Print progress
#'
#' @return List with M_opt (optimal M) and IC_values
#' @keywords internal
select_M_by_IC <- function(Y, X_B, C, tau, M_max, df, max_iter, tol, n_init, verbose) {

  n_obs <- nrow(Y)
  K <- ncol(Y)
  q <- if (is.null(C)) 0 else ncol(C)

  IC_values <- rep(NA, M_max)

  for (M in 1:M_max) {
    if (verbose) message(sprintf("  Evaluating M = %d", M))

    # Initialize and fit
    init <- initialize_MQR(Y, X_B, C, tau, M, n_init)

    result <- fit_MQR_iterative(Y, X_B, C, tau, M, df,
                                 alpha_init = init$alpha,
                                 xi_init = init$xi,
                                 intercept_init = init$intercept,
                                 gamma_init = init$gamma,
                                 max_iter = max_iter,
                                 tol = tol,
                                 verbose = FALSE)

    loss <- compute_loss(Y, X_B, C, result$alpha, result$xi,
                          result$intercept, result$gamma, tau)

    IC_values[M] <- compute_IC(loss, n_obs, K, M, df, q)
  }

  M_opt <- which.min(IC_values)

  return(list(M_opt = M_opt, IC_values = IC_values))
}


#' @title Order Groups for Identifiability
#'
#' @description Reorder groups so that group labels are consistent.
#'
#' @param alpha Group-specific coefficients
#' @param gamma Group membership
#' @param M Number of groups
#'
#' @return List with reordered alpha and gamma
#' @keywords internal
order_groups <- function(alpha, gamma, M) {

  # Order by first appearance in gamma
  first_appear <- match(1:M, gamma)

  # Handle groups that may be empty
  valid_groups <- which(!is.na(first_appear))
  order_idx <- order(first_appear[valid_groups])

  # Create mapping
  new_order <- rep(NA, M)
  new_order[valid_groups] <- order_idx

  # Reorder alpha
  alpha_new <- alpha
  if (length(valid_groups) > 0) {
    alpha_new[order_idx, ] <- alpha[valid_groups, , drop = FALSE]
  }

  # Reorder gamma
  gamma_new <- gamma
  for (m in valid_groups) {
    gamma_new[gamma == m] <- new_order[m]
  }

  return(list(alpha = alpha_new, gamma = gamma_new))
}


#' @title Simulate Data for MQR
#'
#' @description Generate simulated data following the MQR model structure.
#'
#' @param n Number of subjects
#' @param K Number of traits
#' @param M Number of latent groups
#' @param J Number of time points per subject (default 1)
#' @param p Number of primary covariates (default 1)
#' @param q Number of additional covariates (default 2)
#' @param df B-spline degrees of freedom for true alpha (default 4)
#' @param sigma Error standard deviation (default 1)
#' @param seed Random seed (default NULL)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Y}: Outcome matrix (n*J x K)
#'   \item \code{X}: Primary covariate matrix (n*J x p)
#'   \item \code{C}: Additional covariate matrix (n*J x q)
#'   \item \code{T_vec}: Time vector (length n*J)
#'   \item \code{true_gamma}: True group membership (length K)
#'   \item \code{true_alpha}: True group-specific coefficients
#'   \item \code{true_xi}: True trait-specific coefficients
#' }
#'
#' @examples
#' # Generate data with 3 groups
#' sim_data <- simulate_MQR_data(n = 200, K = 10, M = 3)
#'
#' # Check dimensions
#' dim(sim_data$Y)  # 200 x 10
#' table(sim_data$true_gamma)  # Group sizes
#'
#' @export
simulate_MQR_data <- function(n = 200, K = 10, M = 3, J = 1, p = 1, q = 2,
                               df = 4, sigma = 1, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Total observations
  n_total <- n * J

  # Generate time points
  if (J == 1) {
    T_mat <- matrix(runif(n, 0, 1), nrow = n, ncol = 1)
  } else {
    T_mat <- matrix(NA, nrow = n, ncol = J)
    T_mat[, 1] <- runif(n, 0, 0.3)
    for (j in 2:J) {
      T_mat[, j] <- T_mat[, j-1] + runif(n, 0.05, 0.15)
    }
    T_mat <- pmin(T_mat, 1)
  }
  T_vec <- as.vector(T_mat)

  # Generate primary covariates X
  X <- matrix(rnorm(n_total * p), nrow = n_total, ncol = p)

  # Generate additional covariates C
  C <- matrix(NA, nrow = n_total, ncol = q)
  C[, 1] <- rbinom(n_total, 1, 0.5)  # Binary covariate
  if (q > 1) {
    for (j in 2:q) {
      C[, j] <- rnorm(n_total)
    }
  }

  # True group membership (roughly equal sized groups)
  group_sizes <- rep(K %/% M, M)
  remainder <- K %% M
  if (remainder > 0) {
    group_sizes[1:remainder] <- group_sizes[1:remainder] + 1
  }
  true_gamma <- rep(1:M, times = group_sizes)

  # True group-specific coefficient functions
  # Using different polynomial forms for each group
  B_basis <- splines::bs(T_vec, df = df, intercept = TRUE)

  true_alpha_funcs <- list(
    function(t) 0.3 * (1.2 - 0.4 * t^(1/3)),    # Group 1
    function(t) 0.1 + 0.06 * sqrt(t),            # Group 2
    function(t) -0.3 * (0.1 + 0.5 * t^(1/3))    # Group 3
  )

  # Extend for more groups if needed
  if (M > 3) {
    for (m in 4:M) {
      coef <- runif(1, -0.5, 0.5)
      power <- runif(1, 0.3, 1)
      true_alpha_funcs[[m]] <- function(t, c = coef, pw = power) c * t^pw
    }
  }

  # Compute true alpha coefficients (project onto B-spline basis)
  true_alpha <- matrix(NA, nrow = M, ncol = p * df)
  T_grid <- seq(0, 1, length.out = 100)
  B_grid <- splines::bs(T_grid, df = df, intercept = TRUE)

  for (m in 1:M) {
    alpha_vals <- true_alpha_funcs[[min(m, length(true_alpha_funcs))]](T_grid)
    # Least squares projection onto B-spline basis
    for (j in 1:p) {
      idx <- ((j - 1) * df + 1):(j * df)
      true_alpha[m, idx] <- solve(t(B_grid) %*% B_grid, t(B_grid) %*% alpha_vals)
    }
  }

  # True trait-specific coefficients
  true_xi <- matrix(runif(K * q, -1, 1), nrow = K, ncol = q)
  true_intercept <- runif(K, -1, 1)

  # Generate outcomes Y
  Y <- matrix(NA, nrow = n_total, ncol = K)

  X_B <- create_XB_interaction(X, B_basis)

  for (k in 1:K) {
    m <- true_gamma[k]

    # Mean component
    mu <- true_intercept[k] + X_B %*% true_alpha[m, ] + C %*% true_xi[k, ]

    # Add noise
    Y[, k] <- mu + rnorm(n_total, 0, sigma)
  }

  return(list(
    Y = Y,
    X = X,
    C = C,
    T_vec = T_vec,
    n = n,
    K = K,
    M = M,
    J = J,
    true_gamma = true_gamma,
    true_alpha = true_alpha,
    true_xi = true_xi,
    true_intercept = true_intercept
  ))
}
