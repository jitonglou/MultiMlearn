#### For calculating true propensity score ####
#' @title Calculate true propensity scores in the simulation study
#'
#' @description This function uses `equation (2)` in `Section S.1.1` of the supplementary material
#' to calculate true propensity scores for subjects in the simulated datasets.
#'
#' @param X a covariate matrix of subjects.
#' @param K the number of treatments.
#'
#' @return A matrix of true propensity scores of assigning each subject (row) to each treatment (column).
#'
#' @importFrom magrittr "%>%"
#'
#' @export
pi.true = function(X, K){
  p = ncol(X)
  apply(X, 1, function(x){
    xbeta = t(x) %*% (mapply(
      shift, n = 1:K-1,
      MoreArgs = list(x = c(rep(1,2), rep(0, p-2)), type = "lag"),
      SIMPLIFY = TRUE
    )) # result of t(x) %*% beta
    return(exp(xbeta)/sum(exp(xbeta)))
  }) %>%
    t()
}

#### For calculating main effects ####
#' @title Calculate true main effects in the simulation study
#'
#' @description This function uses `equation (5)` in `Section S.1.2` of the supplementary material
#' to calculate true main effects of covariates on the outcome
#' for subjects in the simulated datasets.
#'
#' @param X a covariate matrix of subjects.
#'
#' @return A vector of true main effects of covariates.
#'
#' @export
mu.true = function(X){
  apply(X, 1, function(x){
    1+x[1]+x[2]
  })
}


#### For calculating interaction effects ####
#' @title Lead/lag for vectors and lists
#'
#' @description `Lead` or `lag` vectors or lists.
#'
#' @param x a vector or list.
#' @param n an integer denoting the offset by which to lead or lag the input.
#' This function only supports n>=0.
#' @param type Default is `"lead"` (look `"forwards"`).
#' The other possible value is `"lag"` (look `"backwards"`).
#'
#' @return An object with the same class as `x` that contains the lead/lag of input `x`.
#'
#' @export
shift = function(x, n=1L, type="lag"){
  if (n == 0){return(x)}
  else if (n > 0){
    len = length(x)
    x_new = vector(class(x), len)
    if (type == "lag"){
      x_new[1:n] = x[(len-n+1):len]
      x_new[(n+1):len] = x[1:(len-n)]
    } else if (type == "lead"){
      x_new[1:(len-n)] = x[(n+1):len]
      x_new[(len-n+1):len] = x[1:n]
    }
    return(x_new)
  }
}

#' @title Calculate true interaction effects in the simulation study
#'
#' @description This function uses `equation (6)` in `Section S.1.2` of the supplementary material
#' to calculate true interaction effects of covariates and treatments on the outcome
#' for subjects in the simulated datasets.
#'
#' @param X A covariate matrix of subjects.
#' @param group A vector of group memberships of subjects.
#' @param shift_func A function for the method used in `equation (6)` to obtain
#' different true interaction effects for subjects in different groups.
#'
#' @return A matrix of true propensity scores.
#' Each row represents each subject.
#' Each column represents the interaction effects of covaraites and each treatment.
#'
#' @importFrom magrittr "%>%"
#'
#' @export
delta.true = function(X, group, shift_func = shift){
  lapply(1:nrow(X), function(i){
    apply(X, 1, function(x){
      return(c(
        0.2+x[1]^2+x[2]^2+x[3]^2+x[4]^2,
        0.2+x[1]^2-x[2]^2+x[3]^2-x[4]^2,
        0.2+x[1]^2-x[2]^2-x[3]^2+x[4]^2,
        0.2-x[1]^2-x[2]^2+x[3]^2-x[4]^2
      ))
    })[,i]
  }) %>%
    mapply(shift_func, x = ., n = group - 1, SIMPLIFY = TRUE) %>%
    t()
}

#### Simulate data ####
#' @title Generate a dataset for estimating individualized treatment rules
#'
#' @description This function generates a simulated dataset for estimating
#' individualized treatment rules.
#' The outcome variable is assumed to follow `equation (1)` in
#' `Section S.1.1` of the supplementary material.
#'
#' @param N the number of subjects.
#' @param p the number of covariates.
#' @param K the number of treatments.
#' @param J the number of subject groups.
#' @param propensity_func a user-defined function that calculates true propensity scores of
#' assigning each subject to each treatment. \cr
#' An example of the accepted format of this function can be found in
#' \code{\link{pi.true}}.
#' @param main_func a user-defined function that calculates main effects of
#' covariates on the outcome for each subject. \cr
#' An example of the accepted format of this function can be found in
#' \code{\link{mu.true}}.
#' @param interaction_func a user-defined function that calculates interaction effects of
#' of covariates and treatments on the outcome for each subject. \cr
#' An example of the accepted format of this function can be found in
#' \code{\link{delta.true}}.
#'
#' @return A matrix containing the following columns of all subjects: \cr
#' \itemize{
#'   \item `ID`: IDs.
#'   \item `cluster`: group memberships.
#'   \item `treatment`: observed/assigned treatments.
#'   \item `reward`: values of the (continuous) outcome.
#'   \item `feature_1`, \ldots, `feature_p`: values of the `p` covariates.
#'   \item `interaction_1`, \ldots, `interaction_K`: values of interaction effects
#'   of covariates and each of the `K` treatments.
#'   \item `pi_1`, \ldots, `pi_K`: propensity scores of assigning each of the `K`
#'   treatments.
#'   \item `eps`: values of random errors.
#'   \item `treatment_opt`: optimal treatments.
#'   \item `reward_opt`: outcome values of the optimal treatments.
#'   \item `pi_true`: propensity scores of assigning the observed treatments.
#'   \item `pi_opt`: propensity scores of assigning the optimal treatments.
#' }
#'
#'
#' @family functions used in Section S.1 of the supplementary material
#' @seealso \code{\link{pi.true}} for `propensity_func`,
#' \code{\link{mu.true}} for `main_func`,
#' and \code{\link{delta.true}} for `interaction_func`.
#'
#' @examples
#' #######################################
#' ## The simulated dataset in Section S.1.2 of
#' ## the supplementary material for N=200
#'
#' set.seed(0)
#' simdata200 = simulate_data(
#'   N = 200, p = 20, K = 4, J = 4,
#'   propensity_func = pi.true, # equation (2) in Section S.1.1
#'   main_func = mu.true, # equation (5) in Section S.1.2
#'   interaction_func = delta.true # equation (6) in Section S.1.2
#' )
#'
#' @importFrom magrittr "%>%"
#' @export
simulate_data = function(N = 200, p = 20, K = 4, J = 4,
                         propensity_func, main_func, interaction_func){
  group = sample(x = 1:J, size = N, replace = TRUE) # group membership for each subject
  eps = rnorm(N, 0, 1) # random errors
  X = matrix(runif(N*p, -1, 1), nrow = N, ncol = p) # covariate matrix, dimension: Nxp
  colnames(X) = paste0("feature_", 1:p)

  ## value of interaction effect for each treatment and each subject
  X_int = interaction_func(X, group)
  colnames(X_int) = paste0("interaction_", 1:K)

  ## optimal treatment, propensity score for each treatment class, assigned treatment
  treatment_opt = apply(X_int, 1, which.max)

  treatment_prob = propensity_func(X, K)
  colnames(treatment_prob) = paste0("pi_", 1:K)

  treatment = apply(treatment_prob, 1, FUN = sample,
                    x = 1:K, size = 1, replace = FALSE)

  treatment_idx = cbind(1:N, treatment)
  treatment_opt_idx = cbind(1:N, treatment_opt)

  return(data.frame(
    ID = 1:N, cluster = group,
    treatment = as.factor(paste0("treatment_", treatment)),
    reward = main_func(X) + X_int[treatment_idx] + eps,
    X,  X_int, treatment_prob, eps,
    treatment_opt = as.factor(paste0("treatment_", treatment_opt)),
    reward_opt = main_func(X) + X_int[treatment_opt_idx] + eps,
    pi_true = treatment_prob[treatment_idx],
    pi_opt = treatment_prob[treatment_opt_idx]
  ))
}


