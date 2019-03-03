#' Process data frame
#'
#' \code{dataprocess} processes the original dataset into useful elements: Y, X, strata, weights, offset.
#' @author  Yan Wang, Tianxi Cai, Chuan Hong
#'

dataprocess = function(mf, Terms, family){
  ##### Untangling strata variables
  strats = attr(Terms, "specials")$strata
  if (length(strats)){
    stemp = survival:::untangle.specials(terms(mf),"strata",1)
    if (length(stemp$vars)==1){
      strata.keep = mf[[stemp$vars]]
    }else{
      strata.keep = strata(mf[,stemp$vars],shortlabel = T)
    }
    strats = as.numeric(strata.keep)
  }
  Y <- model.extract(mf, "response")
  ##
  adrop <- 0
  dropterms = c()
  stemp <- untangle.specials(Terms, "strata", 1)
  if (length(stemp$vars) > 0) {
    hasinteractions <- FALSE
    for (i in stemp$vars) {
      if (any(attr(Terms, "order")[attr(Terms, "factors")[i,
                                                          ] > 0] > 1))
        hasinteractions <- TRUE
    }
    if (!hasinteractions)
      dropterms <- c(dropterms, stemp$terms)
    else adrop <- c(0, match(stemp$var, colnames(attr(Terms,
                                                      "factors"))))
  }
  if (length(dropterms)) {
    temppred <- attr(terms, "predvars")
    Terms2 <- Terms[-dropterms]
    if (!is.null(temppred)) {
      attr(Terms2, "predvars") <- temppred[-(1 + dropterms)]
    }
    X <- model.matrix(Terms2, mf)
    renumber <- match(colnames(attr(Terms2, "factors")),
                      colnames(attr(Terms, "factors")))
    attr(X, "assign") <- c(0, renumber)[1 + attr(X, "assign")]
  }
  else X <- model.matrix(Terms, mf)
  if (family$family == 'Cox PH')
    X = X[,colnames(X)!='(Intercept)']
  ##
  offset <- model.offset(mf)
  if (is.null(offset) | all(offset == 0))
    offset <- rep(0, nrow(mf))
  else if (any(!is.finite(offset)))
    stop("offsets must be finite")
  weights <- model.weights(mf)
  if (is.null(weights)){
    weights <- rep(1, nrow(mf))
  }
  if (!is.null(weights) && any(!is.finite(weights)))
    stop("weights must be finite")
  return(list(Y = Y, X = X, strats = strats, weights = weights, offset = offset))
}
