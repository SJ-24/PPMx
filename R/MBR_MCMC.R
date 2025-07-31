# Function to compute the log of the posterior distribution for 'a'
log_a_posterior <- function(a, theta, b, alpha_a, beta_a) {

  # Log-likelihood from theta ~ Gamma(a, b)
  log_likelihood <- sum(dgamma(theta, shape = a, rate = b, log = TRUE))

  # Log-prior from a ~ Gamma(alpha_a, beta_a)
  log_prior <- dgamma(a, shape = alpha_a, rate = beta_a, log = TRUE)

  # Combine log-likelihood and log-prior
  return(log_likelihood + log_prior)
}

# Metropolis-Hastings algorithm for 'a' with lognormal proposal distribution
update_a_lognormal <- function(current_a, theta, b, alpha_a, beta_a, proposal_sd) {
  # browser()
  # Propose a new value for 'a' using a lognormal distribution
  proposed_a <- exp(rnorm(1, mean = log(current_a), sd = proposal_sd))

  # Compute the log posterior for the current and proposed values of 'a'
  log_pos_current <- log_a_posterior(current_a, theta + .Machine$double.xmin, b, alpha_a, beta_a)
  log_pos_proposed <- log_a_posterior(proposed_a, theta + .Machine$double.xmin, b, alpha_a, beta_a)

  # Compute the acceptance probability adding the log of the proposal densities
  log_acceptance_ratio <- log_pos_proposed - log(proposed_a) -
    log_pos_current + log(current_a)
  acceptance_prob <- min(1, exp(log_acceptance_ratio))

  # Accept or reject the proposal
  if (runif(1) < acceptance_prob) {
    return(proposed_a) # Accept the proposed value
  } else {
    return(current_a) # Retain the current value
  }
}

# Gibbs algorithm for 'b'
update_b <- function(current_b, theta, a, alpha_b, beta_b) {
  # browser()
  rgamma(1, shape = a * length(theta) + alpha_b, rate = sum(theta) + beta_b)
}

log_similarity_function <- function(S, similarity_matrix) {
  # browser()
  n <- length(S)
  if (n == 1) {
    return(0)  # log(1) is 0
  } else {
    # Extract the submatrix corresponding to the cluster
    sub_matrix <- similarity_matrix[S, S]

    # Compute the sum of the upper triangular part (excluding diagonal)
    similarity <- mean(sub_matrix[upper.tri(sub_matrix)])

    # Return the logarithm of the similarity
    return(log(similarity))
  }
}


#Algorithm 8 in Neal (2000) for cluster and theta
update_theta <- function(y, t, theta, a, b, U, M, m, similarity_matrix){
  # browser()

  # Step 1: update cluster label
  for (u in 1:U) {
    if(verbose){
      cat(u, "\n")
    }
    theta_u = theta[u]
    is_unique_u = TRUE

    for(u1 in 1:U){
      if(u1 != u && theta[u1] == theta_u){
        is_unique_u = FALSE
        break
      }
    }

    theta_ = theta[-u]
    unique_theta_ = unique(theta_)
    K_ = length(unique_theta_)

    prob_c = rep(0, K_ + m)

    for(h in 1:K_){
      S_h = c()
      for(u1 in 1:U){
        if(u1 != u && theta[u1] == unique_theta_[h]) S_h = c(S_h, u1)
      }

      # if(length(S_h) == 1){
      #   browser()
      # }

      prob_c[h] = dpois(y[u], t[u] * unique_theta_[h], log = TRUE) +
        log(factorial(length(S_h))) +
        log_similarity_function(c(S_h, u), similarity_matrix) -
        log(factorial(length(S_h) - 1)) -
        log_similarity_function(S_h, similarity_matrix)
    }

    theta_pro = rgamma(m, shape = a, rate = b)

    if(is_unique_u){
      theta_pro[1] = theta_u
    }

    for(h in 1:m){
      prob_c[K_ + h] = dpois(y[u], t[u] * theta_pro[h], log = TRUE) +
        log(M) + log_similarity_function(u, similarity_matrix) - log(m)
    }

    prob_c = exp(prob_c)
    prob_c = prob_c / sum(prob_c)
    h = 1
    r_u = runif(1)
    while(r_u > prob_c[h]){
      h = h + 1
      prob_c[h] = prob_c[h] + prob_c[h-1]
    }

    if(verbose){
      cat("r_u:", r_u, "\n")
      cat("prob_c:", prob_c, "\n")
    }

    if(h <= K_){
      theta[u] = unique_theta_[h]
    }else{
      theta[u] = theta_pro[h - K_]
    }
  }


  # Step 2: update theta
  unique_theta = unique(theta)
  K = length(unique_theta)
  for(k in 1:K){
    S_k = which(theta == unique_theta[k])
    theta[S_k] = rgamma(1, shape = a + sum(y[S_k]), rate = b + sum(t[S_k]))
  }

  return(theta)
}

#' MCMC Sampler for Multi-resolution Binary Regression Model
#'
#' This function implements Markov Chain Monte Carlo (MCMC) sampling for a
#' Poisson likelihood model with Gamma-distributed rate parameters (\eqn{\theta_u})
#' under a covariate-dependent product partition model (PPMx). It supports clustering
#' of experimental units using a similarity matrix derived from covariate profiles.
#' The algorithm uses Gibbs sampling for \eqn{b} and \eqn{\theta}, and Metropolis-Hastings
#' with a log-normal proposal for the hyperparameter \eqn{a}.
#'
#' @param y A numeric vector of observed adverse event counts across experimental units.
#' @param t A numeric vector of total exposure time (e.g., patient-time at risk), matching the length of \code{y}.
#' @param similarity_matrix A symmetric matrix of pairwise similarity values among units, typically computed
#'   using \code{\link{calculate_similarity}}. Guides the clustering behavior in the PPMx prior.
#' @param M Concentration parameter of the product partition model. Controls expected number of clusters.
#' @param m Number of auxiliary components used in Neal's Algorithm 8 to propose new cluster values.
#' @param n_burn Number of burn-in iterations to discard from the MCMC chain.
#' @param n_iter Total number of MCMC iterations (must be greater than \code{n_burn}).
#' @param proposal_sd Standard deviation of the log-normal proposal used to update the shape parameter \eqn{a}.
#' @param alpha_a Shape parameter for the Gamma prior on \eqn{a}, i.e., \eqn{a \sim \text{Gamma}(\alpha_a, \beta_a)}.
#' @param beta_a Rate parameter for the Gamma prior on \eqn{a}.
#' @param alpha_b Shape parameter for the Gamma prior on \eqn{b}, i.e., \eqn{b \sim \text{Gamma}(\alpha_b, \beta_b)}.
#' @param beta_b Rate parameter for the Gamma prior on \eqn{b}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{a_spls}}{Posterior samples of the shape parameter \eqn{a}, excluding burn-in.}
#'   \item{\code{b_spls}}{Posterior samples of the rate parameter \eqn{b}, excluding burn-in.}
#'   \item{\code{theta_spls}}{A matrix of posterior samples for unit-specific rates \eqn{\theta_u}
#'     (rows are samples, columns are units), excluding burn-in.}
#' }
#'
#' @details
#' The sampler implements a hierarchical Bayesian model in which unit-level AE counts
#' are modeled as Poisson distributed, \eqn{y_u \sim \text{Poisson}(t_u \cdot \theta_u)}.
#' The rates \eqn{\theta_u} are shared across clusters defined by a random partition,
#' where clustering is influenced by pairwise similarity of covariate profiles.
#'
#' The PPMx prior is enforced via similarity scores computed outside the function
#' and passed in through the \code{similarity_matrix} argument. Neal's Algorithm 8 is used
#' to update cluster allocations and \eqn{\theta_u} values simultaneously.
#'
#' @references
#' Müller, P., Quintana, F., & Rosner, G. (2011). A product partition model with regression on covariates.
#' Journal of Computational and Graphical Statistics, 20(1), 260–278.
#'
#' Dahl, D. B., Day, R., & Tsai, J. W. (2017). Random partition distribution indexed by pairwise information.
#' Journal of the American Statistical Association, 112(518), 721–732.
#'
#' Neal, R. M. (2000). Markov chain sampling methods for Dirichlet process mixture models.
#' Journal of Computational and Graphical Statistics, 9(2), 249–265.
#'
#' @examples
#' # Set hyperparameters and MCMC settings
#' M <- 2               # Total mass parameter in the PPM
#' m <- 3               # Number of auxiliary components in Algorithm 8
#' proposal_sd <- 0.2   # Proposal SD for Metropolis-Hastings update of 'a'
#' alpha_a <- 1; beta_a <- 1
#' alpha_b <- 1; beta_b <- 1
#' n_burn <- 1000
#' n_iter <- 11000
#' set.seed(123)
#'
#' # Extract observed outcomes from trial data
#' t <- indat_outcome$T       # Total exposure time per unit
#' y <- indat_outcome$SAE     # Observed SAE counts per unit
#'
#' # Assume similarity matrix has already been computed (see calculate_similarity)
#' # unblind_similarity_matrix <- calculate_similarity(...)
#'
#' # Run MCMC sampler
#' unblind_spls <- MBR_MCMC(
#'   y = y,
#'   t = t,
#'   similarity_matrix = unblind_similarity_matrix,
#'   M = M,
#'   m = m,
#'   n_burn = n_burn,
#'   n_iter = n_iter,
#'   proposal_sd = proposal_sd,
#'   alpha_a = alpha_a, beta_a = beta_a,
#'   alpha_b = alpha_b, beta_b = beta_b
#' )
#'
#' # Extract posterior samples of theta
#' unblind_theta_spls <- unblind_spls$theta_spls
#'
#' @export
MBR_MCMC = function(y, t, similarity_matrix,
                    M, m,
                    n_burn = 1000,
                    n_iter = 11000,
                    proposal_sd,
                    alpha_a, beta_a, alpha_b, beta_b){

  U = length(y)
  # Starting value of a, b, theta
  a = 1
  b = 1
  theta = rep(sum(y) / sum(t), U)

  # MCMC
  a_spls = rep(NA, n_iter)
  b_spls = rep(NA, n_iter)
  theta_spls = matrix(NA, ncol = U, nrow = n_iter)

  for(i in 1:n_iter){
    # i = 1
    a = update_a_lognormal(a, theta, b, alpha_a, beta_a, proposal_sd)

    b = update_b(b, theta, a, alpha_a, beta_a)

    theta = update_theta(y, t, theta, a, b, U, M, m, similarity_matrix)

    a_spls[i] = a
    b_spls[i] = b
    theta_spls[i,] = theta
  }

  return(list(a_spls = a_spls[-(1:n_burn)],
              b_spls = b_spls[-(1:n_burn)],
              theta_spls = theta_spls[-(1:n_burn),]))
}
