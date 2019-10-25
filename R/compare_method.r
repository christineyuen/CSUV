# This file includes functions to get comparing methods and comparing method fits

# @import futile.logger

#' Get a list of variable selection methods implemented in the CSUV package
#' @export get.compare.methods
#' @return a list of functions
#' @examples
#' \dontrun{
#' X = matrix(rnorm(1000), nrow = 100)
#' Y = rowSums(X[,1:3])+rnorm(100)
#' lasso.method = get.compare.methods()$lasso
#' lasso.mod = lasso.method(X, Y, intercept = FALSE)
#' print(lasso.mod$est.b)
#' }
get.compare.methods<-function(){
  compare.methods = list(lasso = lm.lasso.min,
                         elastic = lm.elastic.half,
                         relaxo = lm.relaxo,
                         # adapt = lm.lasso.adapt,
                         mcp = lm.mcp,
                         scad = lm.scad)
  return (fit.methods.helper.add.name(compare.methods))
}

#' Get a fitted model selected by cross validation
#' @export lm.cv
#' @param X covariates (n times p matrix, n: number of entries, p: number of covariates)
#' @param Y response (vector with n entries)
#' @param intercept TRUE to fit the data with an intercept, FALSE to fit the data without an intercept
#' @param fit.percent percentage of observations used in fitting in cross validation. Each set of subsampling data will have (n times fit.percent) observations for fitting and n times (1-fit.percent) observations for calculating the mse
#' @param num.repeat number of sets of subsampling data used in cross validation
#' @param method.names vector of method names to be used in cross validation. Choose among "lasso", "elastic", "relaxo", "mcp" and "scad". Default is to use all methods listed above
#' @param num.core number of cores to use. Default is 1 (i.e. no parallel running)
#' @return a list which includes the estimated coefficients (est.b) and the corresponding ordinary least square fit from stats::lm()
#' @examples
#' \dontrun{
#' X = matrix(rnorm(1000), nrow = 100)
#' Y = rowSums(X[,1:3])+rnorm(100)
#' cv.mod = lm.cv(X, Y, intercept = FALSE, fit.percent = 0.5, num.repeat = 50)
#' print(cv.mod$est.b)
#' }
lm.cv<-function(X, Y, intercept, fit.percent, num.repeat,
                method.names = NULL, num.core = 1){
  # the function selects the fitted model by lowest delete-k-cv fit mse
  #
  # Args:
  #   X: covariates
  #   Y: response
  #   intercept: whether the model should have intercept
  #   fit.percent: percentage of observations used in fitting
  #   num.repeat: number of iterations
  #   fit.methods: fitting methods
  #   is.keep.all: whether we keep the fitted models, when ols fit is not available
  #   current.fit: current fits

  fit.methods = get.compare.methods()
  if(!is.null(method.names)){
    fit.methods = fit.methods[method.names]
  }

  cv.result = cv.for.methods(X, Y, intercept, fit.percent, num.repeat,
                             fit.methods, is.keep.all = TRUE,
                             current.fit = NULL, num.core = 1)$cv.result
  mse = lapply(cv.result, function(cv.result.by.fold)
    mean(sapply(cv.result.by.fold, function(one.cv.fold.result) one.cv.fold.result$mse)))
  method.num = which.min(mse)
  est.b = fit.methods[[method.num]](X, Y, intercept)$est.b
  return (list(est.b = est.b, mse = unlist(mse)))
  # return (lm.ols.refit.one(X, Y, intercept, est.b))
}

get.compare.fit<-function(x, y, intercept, method.names, current.compare.fit = NULL){
  if(is.null(method.names)|| method.names == ""){
    return (NULL)
  }
  futile.logger::flog.info("start fitting compare methods")
  methods = shiny::isolate(get.compare.methods())
  fit.names = shiny::isolate(method.names)[which(!(shiny::isolate(method.names)%in%names(shiny::isolate(current.compare.fit))))]
  new.fit = do.call(rbind, lapply(fit.names, function(fit.name)
    methods[[fit.name]](shiny::isolate(x),shiny::isolate(y),
                      intercept = intercept)$est.b))
  rownames(new.fit) = fit.names
  futile.logger::flog.info("finished fitting compare methods")
  return (rbind(current.compare.fit, new.fit))
}
