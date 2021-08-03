#### Difference between outcomes ####
#' @title Calculate the difference between rewards (outcomes)
#'
#' @description This function calculates the difference between outcomes of
#' a set of subjects and the subjects in their matched sets. That is, the
#' value of `Rj - Ri` in `equation (7)` in `Section 2.2` of the manuscript.
#'
#' @param data_MS a matrix that contains outcome, treatment, and feature
#' information of a set of subjects.
#' @param idx_MS a data frame of two columns `ID` and `index`. `ID` records
#' the IDs of subjects in `data_MS`. `index` records the column indices of
#' these subjects in `dist_mat`.
#' @param max_size an integer indicating the upper limit of the sizes of all
#' matched sets. The default setting is `max_size=1` which means finding the
#' nearest neighbors of subjects.
#' @param delta a scalar, as defined in `equation (6)` in `Section 2.2` of
#' the manuscript, indicating the upper limit of distances between a subject
#' and the subjects in its matched set. \cr
#' In future versions, we will extend this argument to a vector which means
#' the upper limit can vary with subjects.
#' @param dist_mat a precalculated matrix of distances between subjects.
#' This matrix must include all subjects in `data_MS`.
#'
#' @return A vector of length `max_size*nrow(data_MS)`.
#' The `(i-1)*max_size+1` to the `i*max_size` elements are the distances,
#' in ascending order, between the ith subject and the subjects in
#' its matched set. If the size of this matched set is smaller than `max_size`,
#' some elements will be `NA`.
#'
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr slice_max
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#'
#' @export
diff_reward = function(data_MS, idx_MS, max_size=1, delta, dist_mat){
  out = matrix(NA, nrow=nrow(data_MS), ncol=max_size)

  for (i in 1:nrow(data_MS)){
    ## for each subject, extract the distances between it and other subjects
    ## who were assigned different treatments
    ms_temp = data.frame(data_MS, distance = dist_mat[idx_MS$index,i]) %>%
      filter(treatment != data_MS$treatment[i] & distance < delta)
    if(nrow(ms_temp)>0){
        ## compute difference between rewards
        out[i,] = ms_temp %>%
          arrange(distance) %>%
          slice_max(max_size, with_ties = FALSE) %>%
          mutate(reward_res_diff = data_MS$reward_res[i] - reward_res) %>%
          select(reward_res_diff) %>% t()
    }
  }
  return(c(t(out))) # return a vector of length max_size*nrow(data_MS), can have NA
}

#### Summarize recommendations ####
#' @title Summarize recommendations for multicategory treatment setup
#'
#' @description This function summarizes recommendations for multicategory
#' treatment setup from the results of binary treatment comparisons, using
#' the majority voting approach.
#'
#' @param data a data frame containing the ID (`ID`), outcome (`reward`),
#' and observed treatment (`treatment`) information of subjects.
#' @param pred a matrix with `K(K-1)/2` columns, where `K` is the total number
#' of treatments. Each column is the recommendations between two treatments
#' for all subjects.
#'
#' @return A data frame with the following 4 columns:
#' \itemize{
#'   \item `ID`: IDs of subjects, same as those in `data`.
#'   \item `reward`: outcome values of subjects, same as those in `data`.
#'   \item `treatment`: observed treatments of subjects, same as those in `data`.
#'   \item `vote`: recommended treatments of subjects that are summarized from `pred`
#'   using the majority voting approach.
#' }
#'
#' @importFrom magrittr "%>%"
#' @importFrom dplyr select
#'
#' @export
summarize_rec = function(data, pred){
  out = data %>% select(ID, reward, treatment)
  out$vote = apply(pred, 1, function(z){
    names(table(z))[which.max(table(z))]
  })
  return(out)
}



#### Empirical value function ####
#' @title Calculate empirical value functions
#'
#' @description This function calculates the empirical value functions as
#' defined in `equation (8)` in `Section 2.2` of the manuscript.
#'
#' @param data a data frame containing the ID (`ID`), outcome (`reward`),
#' and observed treatment (`treatment`) information of subjects.
#' @param pred a matrix with `K(K-1)/2` columns, where `K` is the total number
#' of treatments. Each column is the recommendations between two treatments
#' for all subjects.
#' @param propensity a data frame with `K` columns, where `K` is the total number
#' of treatments. The `k`th column is the propensity scores of assigning
#' the `k`th treatment to subjects.
#'
#' @return A scalar as defined in `equation (8)` in `Section 2.2` of
#' the manuscript.
#'
#' @importFrom magrittr "%>%"
#'
#' @export
evf = function(data, pred, propensity){
  fit_vote = apply(pred, 1, function(z){
    names(table(z))[which.max(table(z))]
  })
  row_idx = which(data$treatment == fit_vote) # TRUE/FALSE
  if (length(row_idx)>0){
    denom = rep(1, length(row_idx)) # denom has the same length as row_idx
    for (i in 1:length(row_idx)){
      denom[i] = propensity[i,fit_vote[i]] # take care about the colname of propensity
    }
    return(sum(data$reward[row_idx]/denom)/sum(1/denom)) # DO NOT USE MEAN()
  } else {return(NA)}
}

#### Weighted SVMs for binary treatment comparisons ####
#' @title Fit weighted kernel support vector machines (SVMs)
#'
#' @description This function fits a set of weighted kernel SVMs
#' for binary treatment comparisons.
#'
#' @param train_data a data frame for subjects in training fold(s), containing
#' the ID (`ID`), outcome (`reward`), outcome residual (`reward_res`), observed
#' treatment (`treatment`), and health feature information.
#' @param test_data a data frame for subjects in test fold(s), containing
#' the ID (`ID`), outcome (`reward`), outcome residual (`reward_res`), observed
#' treatment (`treatment`), and health feature information.
#' @param idx a data frame of two columns `ID` and `index`. `ID` records
#' the IDs of subjects in `train_data`. `index` records the column indices of
#' these subjects in `dist_mat`.
#' @param trts a vector of treatment names.
#' @param max_size an integer indicating the upper limit of the sizes of all
#' matched sets. The default setting is `max_size=1` which means finding the
#' nearest neighbors of subjects.
#' @param delta a scalar, as defined in `equation (6)` in `Section 2.2` of
#' the manuscript, indicating the upper limit of distances between a subject
#' and the subjects in its matched set. \cr
#' In future versions, we will extend this argument to a vector which means
#' the upper limit can vary with subjects.
#' @param dist_mat a precalculated matrix of distances between subjects.
#' This matrix must include all subjects in `train_data`.
#' @param g_func a function that transforms the differences between outcomes
#' of a set of subjects and the subjects in their matched sets to the weights
#' in SVMs. In `equation (7)` in `Section 2.2` of the manuscript,
#' `g(.) = |.|` and the weights are `|Rj-Ri|`.
#' @param kernel the kernel function used in SVMs. Supported argument values
#' can be found in \code{\link[kernlab]{ksvm}} and \code{\link[kernlab]{dots}}.
#' Default: "rbfdot".
#' @param kpar the list of hyper-parameters (kernel parameters). Valid
#' parameters for supported kernels can be found in \code{\link[kernlab]{ksvm}}
#' and \code{\link[kernlab]{dots}}. Default: "automatic".
#' @param C a scalar that is the cost of constraints violation in SVMs. This is
#' the "C"-constant of the regularization term in the Lagrange formulation.
#'
#' @return A list with 2 sublists as follows:
#' \itemize{
#'   Suppose there are `K` treatments in the input data,
#'   \item `model`: a list with `K(K-1)/2` sublists. Each sublist is a weighted
#'   SVM (trained on `train_data`) for the corresponding binary treatment
#'   comparison.
#'   \item `prediction`: a matrix with `K(K-1)/2` columns. Each column is the
#'   recommendations between the corresponding treatment pair for subjects in
#'   `test_data`.
#' }
#'
#' @seealso \code{\link[kernlab]{ksvm}} and \code{\link[kernlab]{dots}} for
#' `kernel`. \code{\link[personalized]{weighted.ksvm}} for fitting weighted SVMs.
#'
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr case_when
#' @importFrom dplyr select
#' @importFrom dplyr semi_join
#' @importFrom personalized weighted.ksvm
#'
#' @export
mlearn.wsvm = function(
  train_data, test_data, idx,
  trts, max_size, delta, dist_mat, g_func,
  kernel="rbfdot", kpar="automatic", C
){
  ## total number of treatments
  K = length(trts)

  ## get indices of columns that will be used for prediction
  col_select = which(!(colnames(train_data) %in% c("ID", "cluster", "treatment", "reward", "reward_res")))

  ## create an empty list to record fitted SVMs for all binary treatment comparisons
  model_list = vector("list", K*(K-1)/2)

  ## create an empty matrix to record predictions on test data for all binary treatment comparisons
  pred_mat = matrix(NA, ncol = K*(K-1)/2, nrow = nrow(test_data))

  ## for each binary treatment comparison, fit a weighted SVM
  for (l in 1:(K-1)){ # treatment l
    for (j in ((l+1):K)){ # treatment j
      idxc = l*(K-1)-l*(l-1)/2-(K-j) # index of the comparison
      names(model_list)[idxc] = paste(trts[l], "vs", trts[j], sep="_")

      ## extract training data and create an indicator feature, trt_bi,
      ## for treatment l and j
      data_forsvm = train_data %>%
        filter(treatment == trts[l] | treatment == trts[j]) %>%
        mutate(trt_bi = case_when(
          treatment == trts[l] ~ 1,
          TRUE ~ -1))

      idx_forsvm = idx %>%
        semi_join(data_forsvm, by="ID")
      # print(head(data_forsvm))
      # print(head(idx_forsvm))

      ## calculate differences of reward residuals between subjects and their matched sets
      reward_res_diff = diff_reward(
        data_MS = data_forsvm,
        idx_MS = idx_forsvm,
        max_size = max_size,
        delta = delta,
        dist_mat = dist_mat
      ) %>% t %>% c # a vector of length max_size*nrow(data_forsvm)

      ## exclude subjects who could not find at least one matchup or
      ## the differences between reward residuals are 0
      flag = !((is.na(reward_res_diff)) | (reward_res_diff==0))

      reward_res_diff = reward_res_diff[flag]
      idx_reward_res_diff = rep(1:nrow(data_forsvm),each=max_size)[flag]
      sizes = idx_reward_res_diff %>% table %>% c # sizes of matched sets of subjects

      data_forsvm = data_forsvm[idx_reward_res_diff,]
      idx_forsvm = idx_forsvm[idx_reward_res_diff,]

      ## save the binary classifier (M-learning for binary treatment comparison)
      model_list[[idxc]] = weighted.ksvm(
        # y=ifelse(data_forsvm$trt_bi*sign(reward_res_diff)>=0,trts[l],trts[j]),
        y=(data_forsvm$trt_bi*sign(reward_res_diff)),
        x=(data_forsvm[, col_select] %>% as.matrix()),
        weights=g_func(reward_res_diff)/sizes,
        kernel=kernel, kpar=kpar,
        C=C, nfolds=1, foldid = NULL
      )

      ## save the classifications for samples in the validate fold
      pred_mat[,idxc] = ifelse(
        predict(model_list[[idxc]],
                newx=(test_data[, col_select] %>% as.matrix()))>=0,
        trts[l], trts[j]
      )
    } # loop j ends
  } # loop l ends

  return(list(model=model_list, prediction=pred_mat))
}

#### CV for tuning parameters in weighted SVMs ####
#' @title Fit cross-validated weighted kernel support vector machines (SVMs)
#'
#' @description This function performs a cross-validation (on tuning parameters)
#' of weighted SVMs for multicategory treatment comparisons.
#'
#' @param data a data frame containing the ID (`ID`), outcome (`reward`),
#' outcome residual (`reward_res`), observed treatment (`treatment`), and health
#' feature information of subjects.
#' @param idx a data frame of two columns `ID` and `index`. `ID` records
#' the IDs of subjects in `train_data`. `index` records the column indices of
#' these subjects in `dist_mat`.
#' @param trts a vector of treatment names.
#' @param max_size an integer indicating the upper limit of the sizes of all
#' matched sets. The default setting is `max_size=1` which means finding the
#' nearest neighbors of subjects.
#' @param delta a scalar, as defined in `equation (6)` in `Section 2.2` of
#' the manuscript, indicating the upper limit of distances between a subject
#' and the subjects in its matched set. \cr
#' In future versions, we will extend this argument to a vector which means
#' the upper limit can vary with subjects.
#' @param dist_mat a precalculated matrix of distances between subjects.
#' This matrix must include all subjects in `train_data`.
#' @param g_func a function that transforms the differences between outcomes
#' of a set of subjects and the subjects in their matched sets to the weights
#' in SVMs. In `equation (7)` in `Section 2.2` of the manuscript,
#' `g(.) = |.|` and the weights are `|Rj-Ri|`.
#' @param kernel the kernel function used in SVMs. Supported argument values
#' can be found in \code{\link[kernlab]{ksvm}} and \code{\link[kernlab]{dots}}.
#' Default: "rbfdot".
#' @param kpar the list of hyper-parameters (kernel parameters). Valid
#' parameters for supported kernels can be found in \code{\link[kernlab]{ksvm}}
#' and \code{\link[kernlab]{dots}}. Default: "automatic".
#' @param nfolds_inner the number of folds in the cross-validation. Values
#' greater than or equal to 3 usually yield better results. Default: 3.
#' @param tuneGrid a data frame of tuning parameter(s). Each column for each
#' parameter. Usually, the first column is the cost of constraints violation
#' ("C"-constant) in SVMs.
#' @param propensity a data frame with `K` columns, where `K` is the total number
#' of treatments. The `k`th column is the propensity scores of assigning
#' the `k`th treatment to subjects.
#'
#' @return A list with 7 sublists as follows:
#' \itemize{
#'   \item `best_fit`: the final weighted SVM using the best tuning parameter(s).
#'   \item `params`: the list of tuning parameter(s) used to train the model.
#'   Same as `tuneGrid`.
#'   \item `best_param`: the best tuning parameter(s).
#'   \item `best_idx`: the index of the best tuning parameter(s) in `tuneGrid`/`params`.
#'   \item `cv_mat`: the matrix of the metric values for the cross-validation.
#'   \item `cv_est`: the cross-validation estimators (row means of `cv_mat`).
#'   \item `foldid_inner`: a vector recording the split of folds.
#' }
#'
#' @export
mlearn.wsvm.tune = function(
  data, idx,
  trts, max_size, delta, dist_mat, g_func,
  kernel="rbfdot", kpar="automatic",
  nfolds_inner=3, tuneGrid, propensity
){
  ## create an empty matrix of the metric (empirical value function) for cross-validation
  cv_mat = matrix(NA, nrow = nrow(tuneGrid), ncol = nfolds_inner)

  foldid_inner = sample(rep(seq(nfolds_inner), length = nrow(data)))

  ## cross-validation for tuning parameters
  for (k in 1:nfolds_inner){
    inner_test = which(foldid_inner == k) # a numeric vector that records indices of subjects in the test fold

    if (nfolds_inner == 1){
      train_idx = inner_test
    } else {
      train_idx = seq(nrow(data))[-inner_test]
    }

    # print(paste("Inner fold",k,"has",length(train_idx),"training samples."))

    ## for each tuning parameter (set), fit a weighted SVM for multicategory treatment recommendation
    for (i in 1:nrow(tuneGrid)){
      fit = try(
        mlearn.wsvm(
          train_data=data[train_idx,], test_data=data[inner_test,], idx,
          trts, max_size, delta, dist_mat, g_func,
          kernel, kpar, C=tuneGrid[i,1]),
        silent = TRUE
      ) # return a list with 2 sublists: model and prediction; or "typr-error" class; not every C can yield converged results

      ## calculate the empirical value function for the test fold
      if (class(fit) == "list"){
        cv_mat[i,k] = evf(data[inner_test,], fit$prediction, propensity[inner_test,])
      }
    } # end of for (i in 1:nrow(tuneGrid))
    # print(paste("Inner fold",k,"is done."))
  } # end of for (k in 1:nfolds_inner)

  cv_est = rowMeans(cv_mat, na.rm = TRUE) # obtain cross-validation estimators
  best_idx = which.max(cv_est) # obtain the index of the best tuning parameter(s)
  best_param = tuneGrid[best_idx,] # obtain the value(s) of the best tuning parameter(s)


  # print("C is selected.")
  # print(paste0("Best C is 2^(", log2(best_param), ")."))

  if (all(is.na(cv_est))){
    ## if no cross-validation result converges, set the final model to NULL
    best_fit = NULL
  } else if (nfolds_inner == 1){
    best_fit = fit
  } else {
    ## get the final model by training on the whole dataset using the best tuning parameter(s)
    best_fit = try(
      mlearn.wsvm(
        train_data=data, test_data=data, idx,
        trts, max_size, delta, dist_mat, g_func,
        kernel, kpar, C=best_param),
      silent = TRUE
    )

    ## even using each fold converged, using the whole dataset might not converge
    ## if using the whole dataset does not converge, train on the fold with the maximum EVF to get the final model
    if (class(best_fit) == "try-error") {
      test_idx = which(foldid_inner == which.max(cv_mat[best_idx, ]))
      best_fit = mlearn.wsvm(
        train_data=data[-test_idx,], test_data=data[test_idx,], idx,
        trts, max_size, delta, dist_mat, g_func,
        kernel, kpar, C=best_param)
    }
  } # end of !(is.null(cv_mat) | all(is.na(cv_mat)))

  # print("best_fit is done.")
  # print(summary(best_fit))

  return(list(best_fit = best_fit, params = tuneGrid,
              best_param = best_param, best_idx = best_idx,
              cv_mat = cv_mat, cv_est = cv_est, foldid_inner = foldid_inner)
  )
}


#### CV for ITR ####
#' @title Fit a nested cross-validation of weighted kernel support vector
#' machines (SVMs)
#'
#' @description This function performs a nested cross-validation
#' (on tuning parameters) of weighted SVMs for multicategory treatment
#' comparisons and estimating individualized treatment rules.
#'
#' @param data a data frame containing the ID (`ID`), outcome (`reward_res`),
#' observed treatment (`treatment`), and health feature information of subjects.
#' @param idx a data frame of two columns `ID` and `index`. `ID` records
#' the IDs of subjects in `train_data`. `index` records the column indices of
#' these subjects in `dist_mat`.
#' @param trts a vector of treatment names.
#' @param max_size an integer indicating the upper limit of the sizes of all
#' matched sets. The default setting is `max_size=1` which means finding the
#' nearest neighbors of subjects.
#' @param delta a scalar, as defined in `equation (6)` in `Section 2.2` of
#' the manuscript, indicating the upper limit of distances between a subject
#' and the subjects in its matched set. \cr
#' In future versions, we will extend this argument to a vector which means
#' the upper limit can vary with subjects.
#' @param dist_mat a precalculated matrix of distances between subjects.
#' This matrix must include all subjects in `train_data`.
#' @param g_func a function that transforms the differences between outcomes
#' of a set of subjects and the subjects in their matched sets to the weights
#' in SVMs. In `equation (7)` in `Section 2.2` of the manuscript,
#' `g(.) = |.|` and the weights are `|Rj-Ri|`.
#' @param kernel the kernel function used in SVMs. Supported argument values
#' can be found in \code{\link[kernlab]{ksvm}} and \code{\link[kernlab]{dots}}.
#' Default: "rbfdot".
#' @param kpar the list of hyper-parameters (kernel parameters). Valid
#' parameters for supported kernels can be found in \code{\link[kernlab]{ksvm}}
#' and \code{\link[kernlab]{dots}}. Default: "automatic".
#' @param nfolds_outer the number of folds in the outer layer of the nested
#' cross-validation. Default: 3.
#' @param nfolds_inner the number of folds in the inner layer (for tuning
#' parameters) of the nested cross-validation. Values greater than or equal
#' to 3 usually yield better results. Default: 3.
#' @param tuneGrid a data frame of tuning parameter(s). Each column for each
#' parameter. Usually, the first column is the cost of constraints violation
#' ("C"-constant) in SVMs.
#' @param propensity a data frame with `K` columns, where `K` is the total number
#' of treatments. The `k`th column is the propensity scores of assigning
#' the `k`th treatment to subjects.
#' @param foldit_outer (optional) a user-specified vector recording the split
#' of folds in the outer layer of the nested cross-validation. This vector
#' should match the number of rows in `data` and the number of treatments in
#' `trts`.
#'
#' @return A list with 3 sublists as follows:
#' \itemize{
#'   \item `fit`: a list with `nfolds_outer` sublists. The `j`th sublist contains
#'   the inner cross-validation result of the weighted SVM that used the `j`th
#'   fold of subjects as the test fold.
#'   \item `foldid_outer`: a vector recording the split of folds in the
#'   outer layer of the nested cross-validation.
#'   \item `prediction`: a matrix with 5 columns recording ID(`ID`),
#'   outcome (`reward`), observed treatment (`treatment`), recommended
#'   treatment (`vote`), and the fold in the cross-validation(`fold`)
#'   information of subjects.
#' }
#'
#' @importFrom dplyr arrange
#' @importFrom pracma combs
#' @export
mlearn.wsvm.cv = function(
  data, idx,
  trts, max_size,  delta, dist_mat, g_func,
  kernel="rbfdot", kpar="automatic",
  nfolds_outer=3, nfolds_inner=3, tuneGrid, propensity,
  foldid_outer=NULL
){
  ## total number of treatments
  K = length(trts)
  combs_idx = pracma::combs(1:K,2)

  ## get indices of columns that will be used for prediction
  col_select = which(!(colnames(data) %in% c("ID", "cluster", "treatment", "reward", "reward_res")))

  ## create empty objects for the output
  fit_list = vector("list", nfolds_outer)
  data_pred = NULL

  if (is.null(foldid_outer)){
    foldid_outer = sample(rep(seq(nfolds_outer), length = nrow(data)))
  } else if (length(foldid_outer != nfolds_outer) | any(sort(unique(foldid_outer$fold)) != seq(nfolds_outer))){
    print("ERROR. Check the length of foldid_outer or the number of folds in foldid_outer.")
    return(NULL)
  }

  for (k in 1:nfolds_outer){
    outer_test = which(foldid_outer == k) # a numeric vector recording index

    if (nfolds_outer == 1){
      train_idx = outer_test
    } else {
      train_idx = seq(nrow(data))[-outer_test]
    }
    print(paste("Outer fold",k,"has",length(train_idx),"training samples."))

    ## fit weighted SVMs with tuning parameters on training folds
    fit = mlearn.wsvm.tune(
      data=data[train_idx,], idx=idx[train_idx,],
      trts, max_size, delta, dist_mat, g_func,
      kernel, kpar,
      nfolds_inner, tuneGrid, propensity=propensity[train_idx,])

    # print("fit is done.")
    # print(summary(fit))

    ## obtain predictions on the test fold
    if (nfolds_outer == 1){
      data_pred = try(
        data.frame(summarize_rec(data, fit$best_fit$prediction), fold=1),
        silent = TRUE
      )
    } else {
      data_pred = try(
        fit$best_fit$model %>%
          sapply(function(model){
            predict(model, newx = data[outer_test, col_select] %>% as.matrix)
          }) %>%
          apply(1, function(pred, trt1, trt2){ifelse(pred>=0, trt1, trt2)},
                trt1 = trts[combs_idx[,1]], trt2 = trts[combs_idx[,2]]
                ) %>%
          t %>%
          summarize_rec(data=data[outer_test,], pred=.) %>%
          data.frame(fold = k) %>%
          rbind(data_pred, .),
        silent = TRUE
      )
    }
    # print("pred is done.")

    fit_list[[k]] = fit

    print(paste("Outer fold",k,"is done."))
  } # end of for (k in 1:nfolds_outer)


  names(fit_list) = paste0("fold", seq(nfolds_outer))

  return(list(fit = fit_list, foldid_outer = foldid_outer,
              prediction = try(data_pred %>% arrange(ID), silent = TRUE))
  )
}

