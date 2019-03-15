# Description

  Fast divide-and-conquer Cox proportional hazards model with adaptive lasso, developed by [Yan Wang](https://www.researchgate.net/profile/Yan_Wang374) \<yaw719@mail.harvard.edu\>, [Tianxi Cai](https://www.hsph.harvard.edu/tianxi-cai/), and [Chuan Hong](https://dbmi.hms.harvard.edu/person/postdoctoral-fellows/chuan-hong)

  The dcalasso package aims to fit large adaptive lasso for Cox model when n>>p, by divide-and-conquer strategy, least square approximation, and one-step estimation (divide-and-conquer adaptive lasso estimator). To compute a divide-and-conquer adaptive lasso estimate, the function first finds the estimate for the corresponding unpenalized Cox proportional hazards model and then shrinks it with the adaptive lasso-penalized penalty. This divide-and-conquer adaptive lasso penalty has variable selection consistency and asymptotic normality property as the standard adaptive lasso estimator.
  
  The package can fit both adaptive lasso for time-independent Cox proportional hazards model and adaptive lasso for time-dependent Cox proportional hazards model.  
  
  The dcalasso is specialized in computing the divide-and-conquer adaptive lasso with fast computation, with the following novel features: (1) a divide-and-conquer strategy is applied for the computation of the initial unpenalized Cox proportional hazards model by splitting the observations into <tt>K</tt> chunks and processing each separately; (2) a fast linearization is applied in the estimation of the shrinkage step further reducing the computational burden.
  
  When the interest is to estimate the Cox model without adaptive lasso shrinkage, this package is also useful and the user can simply take the unpenalized estimate result. Note that when n>>p, the computation for the adaptive lasso estimate takes a very small amount of time, i.e. the dcalasso is still computationally advantageous over the standard approach to fitting Cox model for extremely large dataset when the penalized estimate is not of interest.

# Citation

  Wang, Yan, Chuan Hong, Nathan Palmer, Qian Di, Joel Schwartz, Isaac Kohane, and Tianxi Cai. "[A Fast Divide-and-Conquer Sparse Cox Regression.](https://arxiv.org/pdf/1804.00735.pdf)" arXiv preprint arXiv:1804.00735 (2018).

# Installation

  <tt>require(devtools)</tt>
  
  <tt>install_github("michaelyanwang/dcalasso")</tt>
  
  <tt>require(dcalasso)</tt>

# Examples
See <tt>?dcalasso</tt>

# Key function
<tt> dcalasso(formula, family=cox.ph(), data = NULL, data.rds = NULL, weights, subsets, na.action, offset, lambda = 10^seq(-10,3,0.01), gamma = 1, K = 20, iter.os = 2, ncores = 1) </tt>

where <tt>formula</tt> is a formula of a Cox model (see <tt>coxph</tt>), <tt>family</tt> specifies the family of the outcome (<tt>cox.ph</tt> currently), <tt>data</tt> specifies the dataset, if <tt>data</tt> is too large to load into the memory, <tt>data.rds</tt> specifies a vector of file locations where the data are split and saved, <tt>weights, subsets, na.action, offset</tt> are the same as the corresponding arguments in <tt>coxph</tt>, <tt>lambda</tt> is the penalization parameter for adaptive lasso (see <tt>glmnet</tt>), <tt>gamma</tt> is the exponent for the penalization of the adaptive penalty (see <tt>glmnet</tt>), <tt>K</tt> is the number of split which will be overwritten by <tt>length(data.rds)</tt> if <tt>data.rds</tt> is specified, <tt>iter.os</tt> is the number of iterations for the computation of the initial unpenalized estimator (default 2; a larger <tt>iter.os</tt> will result in longer computation but an unpenalized estimator closer to the results at convergence, and <tt>ncores</tt> is the number of cores used in the computation.
