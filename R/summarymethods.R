#' Summary method for dcalasso objects
#'
#' \code{summary.dcalasso} summarizes output of dcalasso fit
#' @param object a dcalasso fit
#' @param unpen whether to print out the unpenalized result
#' @export
#' @author Yan Wang, Tianxi Cai, Chuan Hong

summary.dcalasso = function(object, unpen = F, ...){
  stopifnot(inherits(object, "dcalasso"))
  coef.u <- object$coefficients.unpen
  coef.p <- object$coefficients.pen
  s.err.u <- sqrt(diag(object$cov.unpen))
  s.err.p <- sqrt(diag(object$cov.pen))
  tvalue.u <- coef.u/s.err.u
  tvalue.p <- coef.p/s.err.p

  coef.table.p <- cbind(coef.p, s.err.p, tvalue.p, 2 * pnorm(-abs(tvalue.p)))
  coef.table.u <- cbind(coef.u, s.err.u, tvalue.u, 2 * pnorm(-abs(tvalue.u)))

  colnames(coef.table.p) = c("Penalized Est", "Std. Error", "z value", "Pr(>|z|)")
  colnames(coef.table.u) = c("Unpenalized Est", "Std. Error", "z value", "Pr(>|z|)")

  ans = list(coef.table.p = coef.table.p, coef.table.u = coef.table.u,
             cov.pen = object$cov.pen, cov.unpen = object$cov.unpen,
             K = object$K, n = object$n, n.pen = object$n.pen, family = object$family$family,
             iter = object$iter, BIC = object$BIC.opt, lambda = object$lambda.opt)
  class(ans) <- 'summary.dcalasso'
  return(ans)
}


#' Print summary for dcalasso objects
#'
#' \code{print.summary.dcalasso} summarizes output of dcalasso fit
#'
#' @param object a summary.dcalasso object
#' @param unpen whether to print out the unpenalized result
#' @export
#' @author Yan Wang, Tianxi Cai, Chuan Hong

print.summary.dcalasso = function(object, unpen = F, ...){
  stopifnot(inherits(object, "summary.dcalasso"))
  cat("\nDivide-and-conquer adaptive lasso for a ", object$family," model, n=",
      object$n,".\n")
  cat("\nInitial estimator computed for K=",object$K, "and one-step estimation with",
      object$iter,"iterations.\n")
  cat("\nPenalized summary:\n")
  printCoefmat(object$coef.table.p, signif.stars=T, na.print = ".")
  if (unpen){
    cat("\nUnpenalized summary:\n")
    printCoefmat(object$coef.table.u, signif.stars=T, na.print = "NA")
  }
  cat("\nBIC = ", object$BIC," with lambda = ",object$lambda,"\n")
}


#' Print dcalasso objects
#'
#' \code{print.dcalasso} summarizes output of dcalasso fit
#'
#' @param object a dcalasso object
#' @export
#' @author Yan Wang, Tianxi Cai, Chuan Hong
print.dcalasso = function(object, ...){
  summary.dcalasso(object)
}


#' Plot BIC paths for dcalasso objects
#'
#' \code{plot.dcalasso} summarizes output of dcalasso fit
#'
#' @param object a dcalasso object
#' @export
#' @author Yan Wang, Tianxi Cai, Chuan Hong
plot.dcalasso = function(object, ...){
  plot(log10(object$lambda), object$BIC, xlab = 'Log10 lambda', ylab = 'BIC')
}


#' Extract variance covariance from a dcalasso objects
#'
#' \code{vcov.dcalasso} extracts variance covariance objects
#'
#' @param object a dcalasso object
#' @param unpen whether to switch to the unpenalized variance covariance
#' @export
#' @author Yan Wang, Tianxi Cai, Chuan Hong
vcov.dcalasso = function(object, unpen=F, ...){
  if (unpen){
    return(object$cov.unpen)
  }else{
    return(object$cov.pen)
  }
}


#' Extract coefficients from dcalasso objects
#'
#' \code{coef.dcalasso} extracts coefficients from dcalasso objects
#'
#' @param object a dcalasso object
#' @param unpen whether to switch to the unpenalized coefficients
#' @export
#' @author Yan Wang, Tianxi Cai, Chuan Hong
coef.dcalasso = function(object, unpen=F, ...){
  if (unpen){
    return(object$coefficients.unpen)
  }else{
    return(object$coefficients.pen)
  }
}


#' Prediction of dcalasso object
#'
#' \code{predict.dcalasso} makes prediction of a dcalasso object based on the adaptive lasso estimation.
#' @param object a dcalasso object
#' @param newdata a new data frame
#' @param type "terms", "link", "response" same as predict.glm
#' @export
#' @author Yan Wang, Tianxi Cai, Chuan Hong

predict.dcalasso = function(object, newdata, type = 'link'){
  Terms = delete.response(object$Terms)
  m = model.frame(Terms, newdata)
  X = model.matrix(Terms, m)
  X = X[,names(object$coefficients.pen)]
  if (type == 'terms'){
    return(list(fit = X * object$coefficients.pen,
                se.fit = X * sqrt(diag(object$cov.unpen))))
  }else if (type == 'link'){
    return(list(fit = X %*% object$coefficients.pen,
                se.fit = X %*% sqrt(diag(object$cov.unpen))))
  }else if (type == 'response'){
    return(fit = object$family$linkinv(X %*% object$coefficients.pen))
  }
}
