#'
#'\tabular{ll}{
#'Package: \tab dcalasso\cr
#'Type: \tab Package\cr
#'Version: \tab 0.1\cr
#'Date: \tab 2018-04-09 \cr
#'License: \tab GPL (>= 2)\cr
#'LazyLoad: \tab yes\cr
#'}
#' Divide-and-conquer adaptive lasso for Cox model based on divide and conquer and least square approximation
#'
#'@name dcalasso-package
#'@aliases dcalasso
#'@docType package
#'@title Divide-and-Conquer adaptive lasso for Cox model with least square approximation
#'@author Yan Wang  \email{yaw719@mail.harvard.edu}, Tianxi Cai \email{tcai@hsph.harvard.edu}, Chuan Hong <Chuan_Hong@hms.harvard.edu>
#'@references
#'the paper
#'@import survival mgcv glmnet doParallel foreach MASS
#'@useDynLib dcalasso ScoreAFUNC ScoreAFUNC_AG
#'@keywords divide-and-conquer least square approximation
