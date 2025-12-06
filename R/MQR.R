#' @title Multi-phenotype Quantile Regression with Latent Groups
#'
#' @description
#' Main function to fit joint time-varying quantile regressions with multi-outcome

#' latent groups. This method simultaneously identifies latent group structures and
#' estimates time-varying covariate effects using B-spline approximation.
#'
#' @param Y A matrix of outcomes (n*J rows by K columns), where n is the number of
#'   subjects, J is the number of time points per subject, and K is the number of traits.
#' @param X A matrix of primary covariates with time-varying effects (n*J rows by p columns).
#' @param C A matrix of covariates with constant effects (n*J rows by q columns).
#' @param T_vec A vector of time points (length n*J).
#' @param tau The quantile level (default 0.5).
#' @param M The number of latent groups. If NULL, will be selected by information criterion.
#' @param M_max Maximum number of groups to consider when M is NULL (default 5).
#' @param df Degrees of freedom for B-spline basis (default 4).
#' @param max_iter Maximum number of iterations (default 100).
#' @param tol Convergence tolerance (default 1e-6).
#' @param n_init Number of random initializations for k-means (default 5).
#' @param standardize Logical, whether to standardize outcomes by marginal ranks (default TRUE).
#' @param verbose Logical, whether to print progress (default FALSE).
#'
#' @return A list of class "MQR" containing:
#' \itemize{
#'   \item \code{alpha}: Estimated group-specific time-varying coefficients (M x (p*df) matrix)
#'   \item \code{xi}: Estimated trait-specific coefficients for C (K x q matrix)
#'   \item \code{intercept}: Estimated trait-specific intercepts (length K)
#'   \item \code{gamma}: Estimated group membership for each trait (length K)
#'   \item \code{M}: Number of groups
#'   \item \code{tau}: Quantile level
#'   \item \code{loss}: Final loss value
#'   \item \code{B_basis}: B-spline basis object for prediction
#'   \item \code{converged}: Logical, whether the algorithm converged
#'   \item \code{n_iter}: Number of iterations
#'   \item \code{IC}: Information criterion value (if M was selected)
#' }
#'
#' @details
#' The model assumes:
#' \deqn{Q_{Y_{ik}}(\tau | X_{ik}, T_{ik}, C_{ik}) = \eta_k + X_{ik}^T \alpha_{\gamma_k}(T_{ik}) + C_{ik}^T \xi_k}
#'
#' where \eqn{\eta_k} are trait-specific intercepts, \eqn{\alpha_m(t)} are group-specific
#' time-varying coefficients approximated by B-splines: \eqn{\alpha_m(t) = B(t)^T \theta_m},
#' \eqn{\gamma_k \in \{1,...,M\}} is the group membership of trait k,
#' and \eqn{\xi_k} are trait-specific coefficients for additional covariates.
#'
#' The algorithm alternates between:
#' 1. Updating coefficients (\eqn{\alpha}, \eqn{\xi}, \eqn{\eta}) given group membership
#'    by pooling data within groups and fitting quantile regression
#' 2. Updating group membership (\eqn{\gamma}) given coefficients by assigning each
#'    trait to the group minimizing its check loss
#'
#' Initialization uses k-means clustering on individual trait QR coefficient estimates.
#' Model selection for M uses a BIC-type information criterion:
#' \deqn{IC(M) = L_n/(nK) + log(nK)/(nK) \times n_p(M)}
#'
#' @examples
#' # Generate simulated data
#' set.seed(123)
#' sim_data <- simulate_MQR_data(n = 200, K = 10, M = 3, J = 1)
#'
#' # Fit MQR model
#' fit <- MQR(Y = sim_data$Y, X = sim_data$X, C = sim_data$C,
#'            T_vec = sim_data$T_vec, tau = 0.5, M = 3)
#'
#' # View results
#' print(fit)
#' summary(fit)
#'
#' @references
#' Wang, T., Ma, Y., and Wei, Y. (2024). Joint Time-varying Quantile Regressions
#' with Multi-outcome Latent Groups. \emph{Journal of the Royal Statistical Society}.
#'
#' @export
#' @importFrom quantreg rq
#' @importFrom splines bs
#' @importFrom stats kmeans
MQR <- function(Y, X, C = NULL, T_vec, tau = 0.5, M = NULL,
                M_max = 5, df = 4, max_iter = 100, tol = 1e-6,
                n_init = 5, standardize = TRUE, verbose = FALSE) {

  # Input validation
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.null(C) && !is.matrix(C)) C <- as.matrix(C)

  n_obs <- nrow(Y)
  K <- ncol(Y)
  p <- ncol(X)

  if (length(T_vec) != n_obs) {
    stop("Length of T_vec must equal number of rows in Y")
  }

  # Standardize outcomes by marginal ranks if requested
  if (standardize) {
    Y <- standardize_by_rank(Y, T_vec, df = df)
  }

  # Create B-spline basis
  B_basis <- splines::bs(T_vec, df = df, intercept = TRUE)

  # Create interaction terms: X * B(T)
  X_B <- create_XB_interaction(X, B_basis)

  # Select M by information criterion if not provided
  if (is.null(M)) {
    if (verbose) message("Selecting number of groups by IC...")
    M_select <- select_M_by_IC(Y, X_B, C, tau, M_max, df, max_iter, tol, n_init, verbose, p)
    M <- M_select$M_opt
    if (verbose) message(sprintf("Selected M = %d", M))
  }

  # Get initial values via k-means on individual QR estimates
  init <- initialize_MQR(Y, X_B, C, tau, M, n_init)

  # Run iterative algorithm
  result <- fit_MQR_iterative(Y, X_B, C, tau, M, df,
                               alpha_init = init$alpha,
                               xi_init = init$xi,
                               intercept_init = init$intercept,
                               gamma_init = init$gamma,
                               max_iter = max_iter,
                               tol = tol,
                               verbose = verbose)

  # Compute final loss and IC
  result$loss <- compute_loss(Y, X_B, C, result$alpha, result$xi,
                               result$intercept, result$gamma, tau)
  result$IC <- compute_IC(result$loss, n_obs, K, M, df, ncol(C), p)

  # Store additional info
  result$M <- M
  result$tau <- tau
  result$K <- K
  result$df <- df
  result$B_basis <- B_basis
  result$call <- match.call()

  class(result) <- "MQR"
  return(result)
}


#' @title Fit MQR with Iterative Algorithm
#'
#' @description Internal function implementing the iterative optimization algorithm.
#'
#' @param Y Outcome matrix
#' @param X_B Covariate-spline interaction matrix
#' @param C Additional covariates matrix
#' @param tau Quantile level
#' @param M Number of groups
#' @param df B-spline degrees of freedom
#' @param alpha_init Initial group-specific coefficients
#' @param xi_init Initial trait-specific coefficients
#' @param intercept_init Initial intercepts
#' @param gamma_init Initial group membership
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param verbose Print progress
#'
#' @return List with estimated parameters
#' @keywords internal
fit_MQR_iterative <- function(Y, X_B, C, tau, M, df,
                               alpha_init, xi_init, intercept_init, gamma_init,
                               max_iter, tol, verbose) {

  K <- ncol(Y)
  n_obs <- nrow(Y)
  p_alpha <- ncol(X_B)

  # Initialize

  alpha <- alpha_init
  xi <- xi_init
  intercept <- intercept_init
  gamma <- gamma_init

  loss_old <- Inf
  converged <- FALSE
  n_iter <- 0

  for (iter in 1:max_iter) {
    n_iter <- iter
    gamma_old <- gamma
    alpha_old <- alpha

    # Step 1: Update coefficients given group membership
    update_coef <- update_coefficients(Y, X_B, C, gamma, tau, M)
    alpha <- update_coef$alpha
    xi <- update_coef$xi
    intercept <- update_coef$intercept

    # Step 2: Update group membership given coefficients
    gamma <- update_membership(Y, X_B, C, alpha, xi, intercept, tau, M)

    # Check convergence
    loss_new <- compute_loss(Y, X_B, C, alpha, xi, intercept, gamma, tau)

    if (verbose && iter %% 10 == 0) {
      message(sprintf("Iteration %d: loss = %.6f", iter, loss_new))
    }

    # Convergence criteria
    gamma_stable <- all(gamma == gamma_old)
    alpha_stable <- max(abs(alpha - alpha_old)) < tol
    loss_stable <- abs(loss_old - loss_new) < tol

    if (gamma_stable && (alpha_stable || loss_stable)) {
      converged <- TRUE
      if (verbose) message(sprintf("Converged at iteration %d", iter))
      break
    }

    loss_old <- loss_new
  }

  if (!converged && verbose) {
    warning("Maximum iterations reached without convergence")
  }

  # Reorder groups for identifiability (by first coefficient)
  reorder <- order_groups(alpha, gamma, M)
  alpha <- reorder$alpha
  gamma <- reorder$gamma

  return(list(
    alpha = alpha,
    xi = xi,
    intercept = intercept,
    gamma = gamma,
    converged = converged,
    n_iter = n_iter
  ))
}


#' @title Update Coefficients Given Group Membership
#'
#' @description Update alpha, xi, and intercepts by pooling data within groups.
#'
#' @param Y Outcome matrix
#' @param X_B Covariate-spline interaction matrix
#' @param C Additional covariates matrix
#' @param gamma Current group membership
#' @param tau Quantile level
#' @param M Number of groups
#'
#' @return List with updated alpha, xi, intercept
#' @keywords internal
update_coefficients <- function(Y, X_B, C, gamma, tau, M) {

  K <- ncol(Y)
  n_obs <- nrow(Y)
  p_alpha <- ncol(X_B)
  q <- if (is.null(C)) 0 else ncol(C)

  alpha <- matrix(0, nrow = M, ncol = p_alpha)
  xi <- matrix(0, nrow = K, ncol = max(q, 1))
  intercept <- rep(0, K)

  for (m in 1:M) {
    g_idx <- which(gamma == m)

    if (length(g_idx) == 0) next

    if (length(g_idx) == 1) {
      # Single trait in group - fit individual QR
      k <- g_idx[1]
      fit <- fit_single_trait_qr(Y[, k], X_B, C, tau)
      alpha[m, ] <- fit$alpha
      xi[k, ] <- fit$xi
      intercept[k] <- fit$intercept

    } else {
      # Multiple traits - pool data within group
      fit <- fit_pooled_qr(Y, X_B, C, g_idx, tau)
      alpha[m, ] <- fit$alpha
      xi[g_idx, ] <- fit$xi
      intercept[g_idx] <- fit$intercept
    }
  }

  return(list(alpha = alpha, xi = xi, intercept = intercept))
}


#' @title Update Group Membership Given Coefficients
#'
#' @description Assign each trait to the group minimizing the check loss.
#'
#' @param Y Outcome matrix
#' @param X_B Covariate-spline interaction matrix
#' @param C Additional covariates matrix
#' @param alpha Group-specific coefficients
#' @param xi Trait-specific coefficients
#' @param intercept Trait-specific intercepts
#' @param tau Quantile level
#' @param M Number of groups
#'
#' @return Updated group membership vector
#' @keywords internal
update_membership <- function(Y, X_B, C, alpha, xi, intercept, tau, M) {

  K <- ncol(Y)
  gamma <- rep(NA, K)
  loss_mat <- matrix(NA, nrow = M, ncol = K)

  for (m in 1:M) {
    X_term <- X_B %*% alpha[m, ]

    for (k in 1:K) {
      residual <- Y[, k] - intercept[k] - X_term
      if (!is.null(C)) {
        residual <- residual - C %*% xi[k, ]
      }
      loss_mat[m, k] <- sum(rho_tau(residual, tau))
    }
  }

  # Assign each trait to group with minimum loss
  for (k in 1:K) {
    gamma[k] <- which.min(loss_mat[, k])
  }

  return(gamma)
}


#' @title Predict from MQR Model
#'
#' @description Predict time-varying coefficients at new time points.
#'
#' @param object An object of class "MQR"
#' @param newT Vector of new time points for prediction
#' @param ... Additional arguments (ignored)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{alpha_T}: Predicted group-specific coefficients at newT (length(newT) x M matrix)
#'   \item \code{newT}: The time points
#' }
#'
#' @export
#' @method predict MQR
predict.MQR <- function(object, newT, ...) {

  if (!inherits(object, "MQR")) {
    stop("Object must be of class 'MQR'")
  }

  # Create B-spline basis at new time points
  B_new <- predict(object$B_basis, newT)

  # Compute alpha(t) for each group
  M <- object$M
  p <- ncol(object$alpha) / object$df  # number of X variables

  alpha_T <- list()
  for (j in 1:p) {
    idx <- ((j - 1) * object$df + 1):(j * object$df)
    alpha_T[[j]] <- B_new %*% t(object$alpha[, idx, drop = FALSE])
    colnames(alpha_T[[j]]) <- paste0("Group", 1:M)
  }

  if (p == 1) {
    alpha_T <- alpha_T[[1]]
  }

  return(list(alpha_T = alpha_T, newT = newT))
}


#' @title Print MQR Object
#'
#' @description Print method for MQR objects.
#'
#' @param x An object of class "MQR"
#' @param ... Additional arguments (ignored)
#'
#' @export
#' @method print MQR
print.MQR <- function(x, ...) {

  cat("\nMulti-phenotype Quantile Regression with Latent Groups\n")
  cat("======================================================\n\n")
  cat(sprintf("Quantile level (tau): %.2f\n", x$tau))
  cat(sprintf("Number of traits (K): %d\n", x$K))
  cat(sprintf("Number of groups (M): %d\n", x$M))
  cat(sprintf("B-spline df: %d\n", x$df))
  cat(sprintf("Converged: %s\n", x$converged))
  cat(sprintf("Iterations: %d\n", x$n_iter))
  cat(sprintf("Final loss: %.4f\n", x$loss))

  cat("\nGroup membership:\n")
  for (m in 1:x$M) {
    traits_in_group <- which(x$gamma == m)
    cat(sprintf("  Group %d: traits %s\n", m, paste(traits_in_group, collapse = ", ")))
  }

  invisible(x)
}


#' @title Summary for MQR Object
#'
#' @description Summary method for MQR objects.
#'
#' @param object An object of class "MQR"
#' @param ... Additional arguments (ignored)
#'
#' @return A summary list (invisibly)
#' @export
#' @method summary MQR
summary.MQR <- function(object, ...) {

  cat("\nSummary of MQR Model\n")
  cat("====================\n\n")

  print(object)

  cat("\nGroup-specific coefficients (alpha):\n")
  print(round(object$alpha, 4))

  cat("\nTrait-specific coefficients (xi):\n")
  print(round(object$xi, 4))

  cat("\nTrait-specific intercepts:\n")
  print(round(object$intercept, 4))

  invisible(list(
    tau = object$tau,
    K = object$K,
    M = object$M,
    gamma = object$gamma,
    alpha = object$alpha,
    xi = object$xi,
    intercept = object$intercept,
    loss = object$loss,
    converged = object$converged
  ))
}
