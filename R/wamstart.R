#' Warm-start fit wrapper
#'
#' \code{warmfit} is to obtain an warm-start of the initial estimator using a subset of data.
#'
#' @param Y outcome vector or matrix
#' @param X design matrix
#' @param strata strata vector
#' @param weights weight for each observation
#' @param offset offset for each observation
#' @param family family of the outcome
#' @param idx indices for the subset to fit
#' @export
#' @author Yan Wang, Tianxi Cai, Chuan Hong
#'

warmfit <- function(Y, X, strata, weights, offset, family, idx){
  if (family$family == 'Cox PH'){
    if (ncol(Y) == 2){
      bini = coxph.fit(x = X[idx,], y = Y[idx,],
                       strata = strata[idx],
                       offset = offset[idx], init = NULL,
                       control = coxph.control(),
                       weights = weights[idx],
                       method = 'efron', rownames = NULL)$coefficients
    }else{
      bini = agreg.fit(x = X[idx,], y = Y[idx,],
                       strata = strata[idx],
                       offset = offset[idx], init = NULL,
                       control = coxph.control(),
                       weights = weights[idx],
                       method = 'efron', rownames = NULL)$coefficients
    }
  }else{
    stop("The current version supports only Cox model.")
  }
  return(bini)
}
