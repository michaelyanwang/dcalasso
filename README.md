# Description

  Fast divide-and-conquer Cox proportional hazards model with adaptive lasso, developed by Yan Wang \<yaw719@mail.harvard.edu\>, Tianxi Cai \<tcai@hsph.harvard.edu\>, and Chuan Hong \<Chuan_Hong@hms.harvard.edu\>

  The dcalasso package reduces the computational burden of fitting large adaptive lasso for Cox model when n>>p, by divide and conquer, least square approximation, and one-step estimation.

# Citation

  Wang, Yan, Chuan Hong, Nathan Palmer, Qian Di, Joel Schwartz, Isaac Kohane, and Tianxi Cai. "A Fast Divide-and-Conquer Sparse Cox Regression." arXiv preprint arXiv:1804.00735 (2018).

# Installation

  <tt>require(devtools)</tt>
  
  <tt>install_github("michaelyanwang/dcalasso")</tt>
  
  <tt>require("dcalasso")</tt>

# Examples
See <tt>?dcalasso</tt>

# Key function
<tt> dcalasso(formula, family=cox.ph(), data = NULL, data.rds = NULL, weights, subsets, na.action, offset, lambda = 10^seq(-10,3,0.01), gamma = 1, K = 20, iter.os = 2, ncores = 1) </tt>

where <tt>formula</tt> is a formula of a Cox model (see <tt>coxph</tt>), <tt>family</tt> specifies the family of the outcome (<tt>cox.ph</tt> currently), <tt>data</tt> specifies the dataset, if <tt>data</tt> is too large to load into the memory, <data.rds> specifies a vector of file locations where the data are split and saved, <tt>weights, subsets, na.action, offset</tt> are the same as the corresponding arguments in <tt>coxph</tt>, <tt>lambda</tt> is the penalization parameter for adaptive lasso (see <tt>glmnet</tt>), <tt>gamma</tt> is the exponent for the penalization of the adaptive penalty (see <tt>glmnet</tt>), <tt>K</tt> is the number of split, <tt>iter.os</tt> is the number of iterations for the computation of the initial unpenalized estimator, and <tt>ncores</tt> is the number of cores used in the computation.
