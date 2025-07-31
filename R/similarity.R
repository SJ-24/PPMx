### **Similarity Functions**
similarity_continuous <- function(x, y, gamma) {
  exp(-((x - y)^2) / gamma^2)
}

similarity_ordinal <- function(x, y, range) {
  if (is.na(x) || is.na(y) || is.na(range)) {
    return(0)
  }
  1 - abs(x - y) / range
}

similarity_binary <- function(x, y) {
  if (is.na(x) || is.na(y)) {
    return(0)
  }
  as.numeric(x == y)
}


similarity_categorical <- function(x, y) {
  1 - abs(x - y)
}

similarity_composite <- function(x, y) {
  # Convert the values to sets
  set_x <- unlist(strsplit(as.character(x), ","))
  set_y <- unlist(strsplit(as.character(y), ","))

  # Compute intersection and union
  intersection_size <- length(intersect(set_x, set_y))
  union_size <- length(union(set_x, set_y))

  # Compute Jaccard similarity
  Jaccard_similarity = intersection_size / union_size
  if (union_size == 0 || Jaccard_similarity == 0) {
    return(0)  # Avoid division by zero; treat as disjoint if no union
  } else if(Jaccard_similarity == 1) {
    return(1)
  } else {
    return(0.5)
  }
}

#' Calculate Pairwise Similarity Matrix for Experimental Units
#'
#' This function computes a pairwise similarity matrix among experimental units
#' based on covariate profiles. The similarities are used in the PPMx model to define
#' a covariate-dependent partition prior for clustering units with similar profiles.
#' The function accommodates multiple covariate types, including binary, continuous,
#' ordinal, categorical, and composite covariates, as well as special handling for
#' intervention and dose combinations, including mixed interventions under blinding.
#'
#' @param cov_data A data frame where each row corresponds to an experimental unit (e.g., a study arm or cohort),
#'   and each column corresponds to a covariate. In the example, this includes variables such as binary indicators,
#'   composite labels, intervention type, dose, categorical and continuous patient characteristics.
#'
#' @param covariate_types A character vector indicating the type of each covariate in \code{cov_data}.
#'   Must match the column order of \code{cov_data}. Acceptable types are:
#'   \itemize{
#'     \item \code{"binary"}: Exact match required (e.g., sex, condition presence).
#'     \item \code{"composite"}: Jaccard-like similarity (e.g., age group labels with multiple categories).
#'     \item \code{"categorical"}: One-hot encoded categorical distance (e.g., study phase).
#'     \item \code{"continuous"}: Gaussian kernel similarity based on numeric distance.
#'     \item \code{"ordinal"}: Linear scaled distance between levels (used internally for dose).
#'     \item \code{"Intervention"}: Special handling for treatment vs. placebo, including mixed arms.
#'     \item \code{"Dose"}: Paired with \code{"Intervention"} to assess similarity of treatment intensity.
#'   }
#'
#' @param weights A numeric vector of weights (one per covariate) that determines the contribution of each
#'   covariate to the overall similarity score. Larger values assign more importance to that covariate.
#'   Note: A weight should be specified for \code{"Intervention"} covariates, as these govern the core treatment
#'   comparisons. The \code{"Dose"} covariate is used internally when \code{"Intervention"} is non-placebo,
#'   but does not require a separate weight entry; it will be skipped during direct similarity computation.
#'   In the example, \code{"Intervention"} receives a high weight (e.g., 10) to emphasize its influence,
#'   while \code{"Dose"} has \code{NA} to indicate it's excluded from the weighted sum.
#'
#' @param bandwidths A numeric vector of bandwidth parameters used only for continuous covariates.
#'   Typically estimated as the pooled standard deviation across units (as shown in the example).
#'   Set \code{NA} for non-continuous covariates.
#'
#' @param mixed_intervention A character vector specifying the components of a "Mixed" intervention arm
#'   under blinding (e.g., \code{c("A", "Placebo")}). Only required if \code{cov_data} includes "Mixed" in
#'   the \code{"Intervention"} column.
#'
#' @param mixed_dose A numeric vector giving the assumed dose levels corresponding to \code{mixed_intervention}.
#'   For "Placebo", use \code{NA}. Used when computing similarity under blinding.
#'
#' @param mixed_allocation A numeric vector (summing to 1) representing the assumed allocation proportions
#'   to each treatment in the \code{mixed_intervention}. This reflects prior knowledge of randomization.
#'   Used to compute expected similarity for a "Mixed" arm relative to other arms.
#'
#' @return A symmetric numeric matrix with dimensions equal to the number of rows in
#'   \code{cov_data}. Each entry (i, j) reflects the normalized similarity between
#'   experimental units i and j.
#'
#' @details
#' The similarity function supports multiple covariate types and accounts for missing values
#' by computing similarity only over observed covariates for each unit pair. For intervention
#' covariates, the function handles mixed arms under blinding by computing a weighted
#' similarity based on assumed allocations to known treatments and corresponding doses.
#' For non-placebo matched interventions, similarity includes both binary agreement
#' and ordinal proximity of doses.
#'
#' This function is central to computing the pairwise similarities \( s(x_u, x_{u'}) \)
#' that are averaged within clusters to define the PPMx similarity function used for
#' posterior inference and decision rules in AE monitoring.
#'
#' @examples
#' # Load unblinded data
#' indat <- read.csv("./data_unblind.csv")
#'
#' # Separate outcome and covariate data
#' indat_outcome <- indat[, c("N", "T", "SAE")]
#' indat_cov <- indat[, setdiff(names(indat), c("N", "T", "SAE"))]
#'
#' # Convert dose column to numeric values
#' indat_cov$Dose <- as.numeric(gsub("[^0-9]", "", indat_cov$Dose))
#'
#' # Define covariate types and weights
#' covariate_types <- c("binary", "composite", "binary", "Intervention", "Dose", "categorical", "continuous")
#' weights <- c(2, 2, 5, 10, NA, 2, 2)
#'
#' # Initialize bandwidths
#' bandwidths <- rep(NA, length(weights))
#' continuous_index <- which(covariate_types == "continuous")
#'
#' # Extract mean and standard deviation for continuous variables
#' for (i in continuous_index) {
#'   means <- as.numeric(sub("(.*)\\(.*", "\\1", indat_cov[, i]))
#'   sds <- as.numeric(sub(".*\\((.*)\\).*", "\\1", indat_cov[, i]))
#'   indat_cov[, i] <- means
#'   bandwidths[i] <- sqrt(sum(indat_outcome$N * sds^2) / sum(indat_outcome$N))
#' }
#'
#' # Compute similarity matrix for unblinded data
#' unblind_similarity_matrix <- calculate_similarity(
#'   cov_data = indat_cov,
#'   covariate_types = covariate_types,
#'   weights = weights,
#'   bandwidths = bandwidths
#' )
#'
#' @export
calculate_similarity <- function(cov_data, covariate_types, weights, bandwidths,
                                 mixed_intervention = NA, mixed_dose = NA, mixed_allocation = NA) {
  n <- nrow(cov_data)
  similarity_matrix <- matrix(0, n, n)  # Initialize similarity matrix
  # mixed_study_index = which(cov_data[, "Study"] == mixed_study)
  # browser()
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # cat("i=", i, "; j=", j, "\n")
      # if(i == 2){
      #   browser()
      # }
      observed_indices <- which(!is.na(cov_data[i, ]) & !is.na(cov_data[j, ]))  # Observed covariates
      Xi <- sum(weights[observed_indices], na.rm = TRUE)  # Normalization term
      if (Xi == 0) next  # Skip if no observed covariates

      similarity_sum <- 0
      for (d in observed_indices) {
        x_i <- cov_data[i, d]
        x_j <- cov_data[j, d]
        covariate_type <- covariate_types[d]

        if (covariate_type == "Dose") {
          next # dose level will be considered with intervention
        } else if (covariate_type == "Intervention") {
          if (x_i == "Mixed") {
            similarity_sum_temp = 0
            nmixed = length(mixed_dose) # number of intervention mixed in the study
            dose_j <- cov_data[j, "Dose"]
            for(nm in 1:nmixed){
              if (mixed_intervention[nm] == x_j) {
                if(x_j == "Placebo"){
                  similarity_sum_temp = similarity_sum_temp + mixed_allocation[nm]
                }else{
                  dose_range <- max(c(cov_data[cov_data[, "Intervention"] == x_j, "Dose"], mixed_dose), na.rm = TRUE)
                  similarity_sum_temp = similarity_sum_temp +
                    mixed_allocation[nm] * similarity_ordinal(mixed_dose[nm], dose_j, dose_range)
                }
              }
            }
          } else if (x_j == "Mixed"){
            similarity_sum_temp = 0
            nmixed = length(mixed_dose) # number of intervention mixed in the study
            dose_i <- cov_data[i, "Dose"]
            for(nm in 1:nmixed){
              if (mixed_intervention[nm] == x_i) {
                if(x_i == "Placebo"){
                  similarity_sum_temp = similarity_sum_temp + mixed_allocation[nm]
                }else{
                  dose_range <- max(c(cov_data[cov_data[, "Intervention"] == x_i, "Dose"], mixed_dose), na.rm = TRUE)
                  similarity_sum_temp = similarity_sum_temp +
                    mixed_allocation[nm] * similarity_ordinal(dose_i, mixed_dose[nm], dose_range)
                }
              }
            }
          } else{
            similarity_sum_temp <- similarity_binary(x_i, x_j)
            if (x_i == x_j && x_j != "Placebo") {
              dose_i <- cov_data[i, "Dose"]
              dose_j <- cov_data[j, "Dose"]
              dose_range <- max(cov_data[cov_data[, "Intervention"] == x_i, "Dose"], na.rm = TRUE)
              similarity_sum_temp <- similarity_sum_temp * similarity_ordinal(dose_i, dose_j, dose_range)
            }
          }
          # if (is.na(dose_i)) {
          #   if(x_i != "Placebo") warning(paste0("Unit ", i, " doesn't have the dose level."))
          # } else if (is.na(dose_j)) {
          #   if(x_j != "Placebo") warning(paste0("Unit ", j, " doesn't have the dose level."))
          # }
          similarity_sum <- similarity_sum + weights[d] * similarity_sum_temp
        } else if (covariate_type == "continuous") {
          similarity_sum <- similarity_sum +
            weights[d] * similarity_continuous(x_i, x_j, bandwidths[d])
        } else if (covariate_type == "binary") {
          similarity_sum <- similarity_sum +
            weights[d] * similarity_binary(x_i, x_j)
        } else if (covariate_type == "categorical"){
          similarity_sum <- similarity_sum +
            weights[d] * similarity_categorical(x_i, x_j)
        } else if (covariate_type == "composite") {
          similarity_sum <- similarity_sum +
            weights[d] * similarity_composite(x_i, x_j)
        }
      }
      similarity_matrix[i, j] <- similarity_sum / Xi
      similarity_matrix[j, i] <- similarity_matrix[i, j]  # Symmetric matrix
    }
  }

  return(similarity_matrix)
}
