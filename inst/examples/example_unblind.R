library(AEPPMx)
# Load unblinded data
example_path <- system.file("extdata", "data_unblind.csv", package = "AEPPMx")
indat <- read.csv(example_path)

# Separate outcome and covariate data
indat_outcome <- indat[, c("N", "T", "SAE")]
indat_cov <- indat[, setdiff(names(indat), c("N", "T", "SAE"))]

# Convert dose column to numeric values
indat_cov$Dose <- as.numeric(gsub("[^0-9]", "", indat_cov$Dose))

# Define covariate types and weights
covariate_types <- c("binary", "composite", "binary", "Intervention", "Dose", "categorical", "continuous")
weights <- c(2, 2, 5, 10, NA, 2, 2)

# Initialize bandwidths
bandwidths <- rep(NA, length(weights))
continuous_var_index = which(covariate_types == "continuous")

# Extract mean and standard deviation for continuous variables
for(i in continuous_var_index){
  vals <- sub("(.*)\\(.*", "\\1", indat_cov[, i])
  bracket_vals <- sub(".*\\((.*)\\).*", "\\1", indat_cov[, i])

  vals <- as.numeric(vals)
  indat_cov[, i] <- vals
  vals_SD <- as.numeric(bracket_vals)
  bandwidths[i] <- sqrt(sum(indat_outcome$N * vals_SD^2) / sum(indat_outcome$N))
}

# Compute similarity matrix for unblinded data
unblind_similarity_matrix <-
  calculate_similarity(indat_cov, covariate_types, weights, bandwidths)

# Set hyperparameters and MCMC settings
M <- 2  # Total mass parameter in PPM
m <- 3  # Number of auxiliary components in MCMC algorithm (Neal 2000)
verbose <- FALSE

# Hyper-prior parameters
alpha_a <- 1
beta_a <- 1
alpha_b <- 1
beta_b <- 1
proposal_sd <- 0.2
# Number of MCMC samples
n_burn <- 1000
n_iter <- 11000

seed = 123
set.seed(seed)

# Extract observed outcomes from trial data
t <- indat_outcome$T
U <- length(t)
y <- indat_outcome$SAE

# Run MCMC sampler
unblind_spls = MBR_MCMC(y, t, unblind_similarity_matrix,
                        M, m,
                        n_burn = 1000,
                        n_iter = 11000,
                        proposal_sd,
                        alpha_a, beta_a, alpha_b, beta_b)
unblind_theta_spls <- unblind_spls$theta_spls

# Define indices for treatment and control arms after unblinding
unblind.ctrl.index <- 1
unblind.trt.index <- c(2, 3)

# Run analysis to compute posterior decision probabilities (E2, E3)
res = MBR_analysis(
  unblind_theta_spls,
  unblind_similarity_matrix,
  delta = 0,
  blind = FALSE,
  unblind.ctrl.index = unblind.ctrl.index,
  unblind.trt.index = unblind.trt.index
)
res$E2_post_prob
res$E3_post_prob
