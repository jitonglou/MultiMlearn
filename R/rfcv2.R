#### Customized function for training random forest ####
#' @title Customized function for training random forest
#'
#' @description `rfcv2` creates a random forest model which has both
#' `mtry` and `ntree` as tuning parameters for cross-validations.
#' This function is an extension to random forest models that are currently
#' supported by the `train` function of the `caret` package as all of those
#' models just use `mtry`.
#'
#' @param type the type of the prediction problem.
#' One of `Regression` and `Classification`.
#'
#' @return A function to be used in the `train` function of the `caret` package.
#'
#' @import randomForest randomForest
#' @examples
#' library(caret)
#' library(randomForest)
#' library(mlbench)
#'
#' #######################################
#' ## Classification Example
#' data(iris)
#'
#' set.seed(0)
#' rf_class_fit = train(Species ~ .,
#'                      data=iris,
#'                      method=rfcv2("Classification"),
#'                      tuneGrid=expand.grid(
#'                        .mtry=seq(1,ncol(iris)-1, 1),
#'                        .ntree=seq(100,500,100)),
#'                      trControl=trainControl(method="cv"))
#' print(rf_class_fit)
#'
#' #######################################
#' ## Regression Example
#' data(BostonHousing)
#'
#' set.seed(0)
#' rf_reg_fit = train(medv ~ .,
#'                    data = BostonHousing,
#'                    method=rfcv2("Regression"),
#'                    tuneGrid=expand.grid(
#'                      .mtry=seq(1,sqrt(ncol(BostonHousing)-1), 1),
#'                      .ntree=seq(100,500,100)),
#'                    trControl=trainControl(method="cv"))
#' print(rf_reg_fit)
#'
#' @export
rfcv2 = function(type){
  if (type %in% c("Classification", "Regression")){
    rf = list(type = type, library = "randomForest", loop = NULL)
    rf$parameters = data.frame(
      parameter = c("mtry", "ntree"),
      class = rep("numeric", 2),
      label = c("mtry", "ntree")
    )
    rf$grid = function(x, y, len = NULL, search = "grid"){}
    rf$fit = function(x, y, wts, param, lev, last, weights, classProbs, ...){
      randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
    }
    rf$predict = function(modelFit, newdata, preProc = NULL, submodels = NULL){
      predict(modelFit, newdata)
    }
    rf$prob = ifelse(type == "Classification",
                     function(modelFit, newdata, preProc = NULL, submodels = NULL){
                       predict(modelFit, newdata, type = "prob")
                     },
                     function(modelFit, newdata, preProc = NULL, submodels = NULL){NULL}
    )
    rf$sort = function(x){x[order(x[,1]),]}
    rf$levels = function(x){x$classes}
  } else {
    rf = list()
    print("Error. Please set type to Classification or Regression.")
  }
  return(rf)
}
