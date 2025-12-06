# MQR: Multi-phenotype Quantile Regression with Latent Groups



## Overview

**MQR** implements a novel quantile regression framework for investigating time-dependent covariate effects on multiple outcomes, where different outcome variables may share common covariate effects. The key innovation is simultaneously identifying latent group structures and estimating covariate effects without the need to model the multivariate distribution of the response vector.

This package implements the methodology described in:

> Wang, T., Ma, Y., and Wei, Y. (2025). **Joint Time-varying Quantile Regressions with Multi-outcome Latent Groups**. under review.

## Features

- **Joint modeling**: Simultaneously estimates time-varying coefficients and identifies latent groups
- **Distribution-free**: No parametric assumptions on outcome distributions
- **B-spline approximation**: Flexible modeling of time-varying effects
- **Model selection**: Information criterion-based selection of the number of groups
- **Scalable**: Efficient algorithm handling thousands of observations

## Installation

### From GitHub (development version)

```r
# Install devtools if not already installed
install.packages("devtools")

# Install MQR from GitHub
devtools::install_github("tianyingw/MQR")
```

### From source

```r
# Clone the repository and install
# git clone https://github.com/tianyingw/MQR.git
install.packages("path/to/MQR", repos = NULL, type = "source")
```

## Quick Start

### Basic Example

```r
library(MQR)

# Simulate data with 3 latent groups
set.seed(123)
sim_data <- simulate_MQR_data(
  n = 200,    # 200 subjects
  K = 10,     # 10 traits
  M = 3,      # 3 latent groups
  J = 1       # 1 time point per subject
)

# Fit MQR model
fit <- MQR(
  Y = sim_data$Y,
  X = sim_data$X,
  C = sim_data$C,
  T_vec = sim_data$T_vec,
  tau = 0.5,  # median regression
  M = 3       # specify 3 groups
)

# View results
print(fit)
summary(fit)

# Compare estimated vs true group membership
table(Estimated = fit$gamma, True = sim_data$true_gamma)
```

### Model Selection

If the number of groups is unknown, MQR can select it using an information criterion:

```r
# Let MQR select M automatically
fit_auto <- MQR(
  Y = sim_data$Y,
  X = sim_data$X,
  C = sim_data$C,
  T_vec = sim_data$T_vec,
  tau = 0.5,
  M = NULL,      # auto-select M
  M_max = 5,     # consider M = 1, 2, ..., 5
  verbose = TRUE
)

print(fit_auto$M)  # Selected number of groups
```

### Prediction

Predict time-varying coefficients at new time points:

```r
# Predict alpha(t) at specific time points
new_times <- seq(0.1, 0.9, by = 0.1)
pred <- predict(fit, newT = new_times)

# Plot the time-varying coefficients
matplot(new_times, pred$alpha_T, type = "l", lty = 1, lwd = 2,
        xlab = "Time", ylab = expression(alpha(t)),
        main = "Estimated Time-Varying Coefficients")
legend("topright", paste("Group", 1:fit$M), col = 1:fit$M, lty = 1, lwd = 2)
```

## Model Details

### The Model

For subject $i$ and trait $k$, the conditional quantile function is:

$$Q_{Y_{ik}}(\tau | X_{ik}, T_{ik}, C_{ik}) = X_{ik} \alpha_{\gamma_k}(T_{ik}) + C_{ik}^T \xi_k$$

where:
- $\alpha_m(t)$ are M distinct time-varying coefficient functions (group-specific)
- $\gamma_k \in \{1,...,M\}$ is the latent group membership of trait $k$
- $\xi_k$ are trait-specific coefficients for additional covariates

### Algorithm

The estimation procedure alternates between:

1. **Update coefficients**: Given group membership, estimate $\alpha$ and $\xi$ by pooling data within groups and fitting quantile regression
2. **Update membership**: Given coefficients, assign each trait to the group minimizing the check loss

The time-varying coefficients are approximated using B-splines:
$$\alpha_m(t) \approx B(t)^T \theta_m$$

### Model Selection

The number of groups M is selected by minimizing a BIC-type information criterion:
$$IC(M) = L_n(\hat{\theta}, \hat{\gamma}) + \frac{\log(nK)}{2nK} \cdot n_p(M)$$

where $L_n$ is the check loss and $n_p(M)$ is the number of parameters.

## Main Functions

| Function | Description |
|----------|-------------|
| `MQR()` | Main function to fit the MQR model |
| `simulate_MQR_data()` | Generate simulated data |
| `predict.MQR()` | Predict time-varying coefficients |
| `print.MQR()` | Print model summary |
| `summary.MQR()` | Detailed model summary |

## Arguments for `MQR()`

| Argument | Description | Default |
|----------|-------------|---------|
| `Y` | Outcome matrix (n*J x K) | Required |
| `X` | Primary covariate matrix | Required |
| `C` | Additional covariates | `NULL` |
| `T_vec` | Time points vector | Required |
| `tau` | Quantile level | `0.5` |
| `M` | Number of groups (NULL for auto-select) | `NULL` |
| `M_max` | Maximum M to consider | `5` |
| `df` | B-spline degrees of freedom | `4` |
| `max_iter` | Maximum iterations | `100` |
| `tol` | Convergence tolerance | `1e-6` |
| `standardize` | Standardize outcomes by rank | `TRUE` |
| `verbose` | Print progress | `FALSE` |

## Output

The `MQR()` function returns an object of class "MQR" containing:

- `alpha`: Estimated group-specific B-spline coefficients (M x df matrix)
- `xi`: Estimated trait-specific coefficients (K x q matrix)
- `intercept`: Estimated trait-specific intercepts (length K)
- `gamma`: Estimated group membership (length K)
- `M`: Number of groups
- `tau`: Quantile level
- `loss`: Final loss value
- `IC`: Information criterion value
- `converged`: Convergence status
- `n_iter`: Number of iterations

## Citation

If you use this package, please cite:

```bibtex
@article{wang2025joint,
  title={Joint Time-varying Quantile Regressions with Multi-outcome Latent Groups},
  author={Wang, Tianying and Ma, Yanyuan and Wei, Ying},
  year={2025}
}
```

## Dependencies

- R (>= 3.5.0)
- quantreg
- splines

## License

GPL-3

## Contact

- **Tianying Wang** - tianying.wang@colostate.edu
- Department of Statistics, Colorado State University

## Issues and Contributions

Please report bugs and feature requests on [GitHub Issues](https://github.com/tianyingw/MQR/issues).

Contributions are welcome via pull requests.
