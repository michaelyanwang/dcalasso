#' Computes score and information given subset
#'
#' \code{scoreinfofun} computes information and scores.
#'
#' @param Y outcome vector or matrix
#' @param X design matrix
#' @param stratas strata vector
#' @param weights weight for each observation
#' @param offset offset for each observation
#' @param family family of the outcome
#' @param bini initial coefficient estimate
#' @param idx indices for evaluation
#' @useDynLib dcalasso ScoreAFUNC ScoreAFUNC_AG
#' @export
#' @author Yan Wang, Tianxi Cai, Chuan Hong
#'

scoreinfofun <- function(Y, X, stratas, weights, offset, family, bini, idx){
  if (family$family == 'Cox PH'){
    if (ncol(Y) == 2){
      # Time independent case
      if (length(stratas) == 0){
        sorted = order(Y[idx,1])
        newstrat = as.integer(rep(0, length(idx)))
      }else{
        sorted = order(stratas[idx],Y[idx,1])
        newstrat <- as.integer(c(1 * (diff(as.numeric(stratas[idx])) != 0), 1))
      }
      # SEXP ScoreAFUNC(SEXP time2,   SEXP status2,  SEXP covar2,  SEXP betahat2,
      #                 SEXP offset2, SEXP weights2, SEXP strata2               )
      tmp = .Call("ScoreAFUNC",
                  as.double(Y[idx[sorted],1]),
                  as.integer(Y[idx[sorted],2]),
                  as.matrix(X[idx[sorted],]),
                  as.double(bini),
                  as.double(offset[idx[sorted]]),
                  as.double(weights[idx[sorted]]),
                  as.integer(newstrat))
    }else{
      # Time dependent case
      if (length(stratas) == 0){
        sort.end <- order(-Y[idx,2]) - 1L
        sort.start <- order(-Y[idx,1]) - 1L
        newstrat = length(idx)
      }else{
        sort.end <- order(stratas[idx], -Y[idx,2]) - 1L
        sort.start <- order(stratas[idx], -Y[idx,1]) - 1L
        newstrat = cumsum(table(stratas[idx]))
      }
      # SEXP ScoreAFUNC_AG(SEXP surv2,      SEXP covar2,   SEXP strata2,
      #                    SEXP weights2,   SEXP offset2,
      #                    SEXP sort12,     SEXP sort22,   SEXP betahat2)
      tmp = .Call("ScoreAFUNC_AG",
                  Y[idx,],
                  X[idx,],
                  as.integer(newstrat),
                  as.double(weights[idx]),
                  as.double(offset[idx]),
                  sort.start,sort.end,as.double(bini))
    }
    score = tmp[[1]]
    info = matrix(tmp[[2]],length(bini),length(bini))
  }else{
    stop("The current version supports only Cox model.")
  }
  return(list(score = score, info = info))
}
