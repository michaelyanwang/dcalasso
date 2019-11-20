# Package name

  dcalasso: Fast divide-and-conquer Cox proportional hazards model with adaptive lasso
  
# Questions that the package addresses

  The package answers the following questions:
  
* I want to fit a Cox proportional hazards model. But my dataset is too large to load to RAM. What should I do?
* I want to fit a Cox proportional hazards model. But my dataset is too large to save as one piece and standard software requires to load the dataset as a whole piece. What should I do?
* I want to do variable selection in Cox model. But my dataset is too large to make the computation feasible, too large to load to RAM, too large to save as one piece for Cox model variable selection method. What should I do?
  ....
  
  Essentially - what should I do when my data for Cox model w/ or w/o variable selection are too large?

# Methodology

  The dcalasso package aims to fit Cox proportional hazards model to extremely large, when both n and p are extremely large and n>>p. The method and package have the following features: 
 * The package tackles the Cox model fitting for extremely large data using the divide-and-conquer strategy, even when the data are too large to save as one file.
 * This approach is able to achieve a fast computation. Meanwhile, it returns a set of results that are close to the precision as if the model was fitted to the dataset as a whole.
 * The package could provide model fitting without variable selection as well as model fitting with variable selection. It returns results for both an unpenalized Cox model without variable selection and an adaptive LASSO-penalized variable selection for the Cox model.
 * The adaptive LASSO variable selection has variable selection consistency and asymptotic normality.
 * The method can be applied to both time-independent and time-dependent survival datasets.
 * The package is flexible in terms of multi-core or single-core computation.
  
  The method is detailed [here](https://academic.oup.com/biostatistics/advance-article-abstract/doi/10.1093/biostatistics/kxz036/5572660?redirectedFrom=fulltext). Briefly, the method first finds a divide-and-conquer Cox model estimate without adaptive LASSO penalty by applying the divide-and-conquer strategy with one-step estimation to the data that are divided into subsets. Then it finds the divide-and-conquer adaptive LASSO estimate based on the divide-and-conquer Cox estimate, using least square approximation. 
    
# Developers

The [method](https://academic.oup.com/biostatistics/advance-article-abstract/doi/10.1093/biostatistics/kxz036/5572660?redirectedFrom=fulltext) was developed by 
* [Yan Wang](https://www.researchgate.net/profile/Yan_Wang374) \<yaw719@mail.harvard.edu\>
* [Chuan Hong](https://dbmi.hms.harvard.edu/person/postdoctoral-fellows/chuan-hong)
* [Nathan Palmer](https://dbmi.hms.harvard.edu/people/nathan-patrick-palmer)
* [Qian Di](https://scholar.google.com/citations?user=BpMY1OkAAAAJ&hl=en)
* [Joel Schwartz](https://www.hsph.harvard.edu/joel-schwartz/)
* [Isaac Kohane](https://dbmi.hms.harvard.edu/people/isaac-samuel-kohane)
* [Tianxi Cai](https://www.hsph.harvard.edu/tianxi-cai/)

# Citation

  Wang, Yan, Chuan Hong, Nathan Palmer, Qian Di, Joel Schwartz, Isaac Kohane, and Tianxi Cai. "[A Fast Divide-and-Conquer Sparse Cox Regression.](https://arxiv.org/pdf/1804.00735.pdf)". 2019 Sep 23. kxz036


# Installation

```
require(devtools)
install_github("michaelyanwang/dcalasso")
require(dcalasso)
```
# Key function
```
dcalasso(formula, family=cox.ph(), data = NULL, data.rds = NULL, weights, subsets, na.action, offset, lambda = 10^seq(-10,3,0.01), gamma = 1, K = 20, iter.os = 2, ncores = 1)
```

* <tt>formula</tt> is a formula of a Cox model (see <tt>coxph</tt>)
* <tt>family</tt> specifies the family of the outcome (<tt>cox.ph</tt> currently)
* <tt>data</tt> specifies the dataset, if <tt>data</tt> is too large to load into the memory, <tt>data.rds</tt> specifies a vector of file locations where the data are split and saved
* <tt>weights, subsets, na.action, offset</tt> are the same as the corresponding arguments in <tt>coxph</tt>
* <tt>lambda</tt> is the penalization parameter for adaptive lasso (see <tt>glmnet</tt>)
* <tt>gamma</tt> is the exponent for the penalization of the adaptive penalty (see <tt>glmnet</tt>)
* <tt>K</tt> is the number of split which will be overwritten by <tt>length(data.rds)</tt> if <tt>data.rds</tt> is specified
* <tt>iter.os</tt> is the number of iterations for the computation of the initial unpenalized estimator (default 2; a larger <tt>iter.os</tt> will result in longer computation but an unpenalized estimator closer to the results at convergence
* <tt>ncores</tt> is the number of cores used in the computation.


# Key examples

1. Time-independent dataset: Fitting a Cox model for a time-independent dataset with 50 variables and 100,000 samples. 

```
# Data simulation
set.seed(1)
N = 1e5; p.x = 50; K = 100; n = N/K;  cor = 0.2;
bb = c(rep(0.4,4),rep(0.2,4),rep(0.1,4),rep(0.05,4))
beta0 = c(1, bb, rep(0, p.x - length(bb)))
dat.mat0 = as.data.frame(SIM.FUN(N, p.x = p.x, cor = cor, family='Cox',beta0 = beta0))
dat.mat0[,'strat'] = rep(1:20, each = N/20)

# Model fitting
modp = dcalasso(as.formula(paste0('Surv(u,delta)~',paste(paste0('V',3:52),collapse='+'))),
                family = 'cox.ph',data = dat.mat0,
                K = 10, iter.os = 4, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)
```
  In this case, the dataset was loaded as a whole (<tt>data=dat.mat0</tt>). The same formulaic syntax for <tt>coxph</tt> applies here. For a time-independent dataset, two arguments are required: time and event.The dcalasso function internally divides it into 10 folds (<tt>K=10</tt>). The divide-and-conquer Cox estimate was estimated using 4 iterations of one-step updates (<tt>iter.os = 4</tt>), with the process paralleled to 2 CPUs (<tt>ncores = 2</tt>).
  
  The print statement provides coefficients for both unpenalized estimate and adaptive LASSO estimate. The plot statement provides the relationship between the penalization factor <tt>lambda</tt> and model's Bayesian information criteria (BIC), which was the metric built in the package for variable selection.


2. Time-dependent dataset: Fitting a Cox model for a dataset with 50 time-dependent variables, 50 additional time-independent variables, and 100,000 samples. 

```
# Data simulation
set.seed(1)
n.subject = 1e5; p.ti = 50; p.tv = 50; K = 20; n = n.subject/K; cor = 0.2;  lambda.grid = 10^seq(-10,3,0.01);
beta0.ti = NULL
beta0.tv = NULL
dat.mat0 = as.data.frame(SIM.FUN.TVC(p.ti, p.tv, n.subject, cor, beta0.ti, beta0.tv))
dat.mat0[,'strat'] = dat.mat0[,dim(dat.mat0)[2]]%%(n.subject/20)
dat.mat0 = dat.mat0[,-(dim(dat.mat0)[2]-1)]

# Model fitting
modp = dcalasso(as.formula(paste0('Surv(t0,t1,status)~',paste(paste0('V',4:103),collapse='+'))),
                family = 'cox.ph',data = dat.mat0,
                K = 10, iter.os = 2, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)
```
  In this case, the dataset was loaded as a whole (<tt>data=dat.mat0</tt>). The same formulaic syntax for <tt>coxph</tt> applies here. For a time-dependent dataset, three arguments are required: start, end, and event. The dcalasso function internally divides it into 10 folds (<tt>K=10</tt>). The divide-and-conquer Cox estimate was estimated using 2 iterations of one-step updates (<tt>iter.os = 2</tt>), with the process paralleled to 2 CPUs (<tt>ncores = 2</tt>).
  
  The print statement provides coefficients for both unpenalized estimate and adaptive LASSO estimate. The plot statement provides the relationship between the penalization factor <tt>lambda</tt> and model's Bayesian information criteria (BIC), which was the metric built in the package for variable selection.
  

# Other examples
See <tt>?dcalasso</tt>

1. Time-independent dataset, when the dataset can be saved and loaded as a whole: Fitting a Cox model for a time-independent dataset with 50 variables and 100,000 samples (detailed version of Key Example 1).

```
##### Generating a time-independent dataset #####
set.seed(1)
N = 1e5; p.x = 50; K = 100; n = N/K;  cor = 0.2;
bb = c(rep(0.4,4),rep(0.2,4),rep(0.1,4),rep(0.05,4))
beta0 = c(1, bb, rep(0, p.x - length(bb)))
dat.mat0 = as.data.frame(SIM.FUN(N, p.x = p.x, cor = cor, family='Cox',beta0 = beta0))
dat.mat0[,'strat'] = rep(1:20, each = N/20)

## A Cox model, without stratification of baseline hazard: Surv(u,delta)~V3+V4+.....+V52
# Example option 1: Using 1 core for computation, dividing the dataset to 10 chunks, using 2 iterations of one-step estimator for update
mod = dcalasso(as.formula(paste0('Surv(u,delta)~',paste(paste0('V',3:52),collapse='+'))),
               family = 'cox.ph',data = dat.mat0,
               K = 10, iter.os = 2)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)

# Performing model prediction
pred.link = predict(mod, newdata = dat.mat0)
pred.term = predict(mod, newdata = dat.mat0, type = 'terms')
pred.response = predict(mod, newdata = dat.mat0, type = 'response')


# Example option 2: Same model as above: Using 2 cores for parallel computation, dividing the dataset to 10 chunks, using 4 iterations of one-step estimator for update
modp = dcalasso(as.formula(paste0('Surv(u,delta)~',paste(paste0('V',3:52),collapse='+'))),
                family = 'cox.ph',data = dat.mat0,
                K = 10, iter.os = 4, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Compare the divide and conquer method unpenalized model fit against a standard Cox fit using the whole dataset
std = coxph(as.formula(paste0('Surv(u,delta)~',paste(paste0('V',3:52),collapse='+'))),
            data = dat.mat0)

plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, std$coefficients)


## A Cox model, with stratification of baseline hazard by "strat": Surv(u,delta)~strata(strat)+V3+V4+.....+V52
# Example option 1: Using 1 core for computation, dividing the dataset to 10 chunks, using 2 iterations of one-step estimator for update
mod = dcalasso(as.formula(paste0('Surv(u,delta)~strata(strat)+',paste(paste0('V',3:52),collapse='+'))),
               family = 'cox.ph',data = dat.mat0,
               K = 10, iter.os = 2)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)

# Example option 2: Same model as above: Using 2 cores for parallel computation, dividing the dataset to 10 chunks, using 4 iterations of one-step estimator for update
modp = dcalasso(as.formula(paste0('Surv(u,delta)~strata(strat)+',paste(paste0('V',3:52),collapse='+'))),
                family = 'cox.ph',data = dat.mat0,
                K = 10, iter.os = 2, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Compare the divide and conquer method unpenalized model fit against a standard Cox fit using the whole dataset
std = coxph(as.formula(paste0('Surv(u,delta)~strata(strat)+',paste(paste0('V',3:52),collapse='+'))),
            data = dat.mat0)


plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, std$coefficients)
```

2. Time-independent dataset, when the dataset are saved in multiple files: Fitting a Cox model for a time-independent dataset with 50 variables and 100,000 samples.
```
##### Generating a time-independent dataset, saving the data into 10 separate files #####
set.seed(1)
N = 1e5; p.x = 50; K = 100; n = N/K;  cor = 0.2;
bb = c(rep(0.4,4),rep(0.2,4),rep(0.1,4),rep(0.05,4))
beta0 = c(1, bb, rep(0, p.x - length(bb)))
dat.mat0 = as.data.frame(SIM.FUN(N, p.x = p.x, cor = cor, family='Cox',beta0 = beta0))
dat.mat0[,'strat'] = rep(1:20, each = N/20)
dir = "C:/"
ll = split(1:N, factor(1:10))
for (kk in 1: 10){
  df = dat.mat0[ll[[kk]],]
  saveRDS(df, file = paste0(dir,'dataTI',kk,'.rds'))
}

## A Cox model, without stratification of baseline hazard: Surv(u,delta)~V3+V4+.....+V52
# Example option 1: Using 1 core for computation, data.rds specifies that the data are saved into 10 files [each time the RAM will only contain data from 1 file], using 2 iterations of one-step estimator for update
mod = dcalasso(as.formula(paste0('Surv(u,delta)~',paste(paste0('V',3:52),collapse='+'))),
               family = 'cox.ph',
               data.rds = paste0(dir,'dataTI',1:10,'.rds'), iter.os = 2)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)

# Example option 2: Using 2 core for parallel computation, data.rds specifies that the data are saved into 10 files [each time the RAM will only contain data from 1 file], using 2 iterations of one-step estimator for update
modp = dcalasso(as.formula(paste0('Surv(u,delta)~',paste(paste0('V',3:52),collapse='+'))),
                family = 'cox.ph',
                data.rds = paste0(dir,'dataTI',1:10,'.rds'), iter.os = 2, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Compare the divide and conquer method unpenalized model fit against a standard Cox fit using the whole dataset
std = coxph(as.formula(paste0('Surv(u,delta)~',paste(paste0('V',3:52),collapse='+'))),
            data = dat.mat0)

plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, std$coefficients)


## A Cox model, with stratification of baseline hazard by "strat": Surv(u,delta)~strata(strat)+V3+V4+.....+V52
# Example option 1: Using 1 core for computation, loading data from 10 separate files (imply K=10) [each time the RAM will only contain data from 1 file], using 2 iterations of one-step estimator for update
mod = dcalasso(as.formula(paste0('Surv(u,delta)~strata(strat)+',paste(paste0('V',3:52),collapse='+'))),
               family = 'cox.ph',
               data.rds = paste0(dir,'dataTI',1:10,'.rds'), iter.os = 2)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)

# Example option 2: Using 2 core for parallel computation, data.rds specifies that the data are saved into 10 files [each time the RAM will only contain data from 1 file], using 2 iterations of one-step estimator for update
modp = dcalasso(as.formula(paste0('Surv(u,delta)~strata(strat)+',paste(paste0('V',3:52),collapse='+'))),
                family = 'cox.ph',
                data.rds = paste0(dir,'dataTI',1:10,'.rds'), iter.os = 2, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Compare the divide and conquer method unpenalized model fit against a standard Cox fit using the whole dataset
std = coxph(as.formula(paste0('Surv(u,delta)~strata(strat)+',paste(paste0('V',3:52),collapse='+'))),
            data = dat.mat0)

plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, std$coefficients)
```

3. Time-dependent dataset, when the dataset can be saved and loaded as a whole: Fitting a Cox model for a dataset with 50 time-dependent variables, 50 additional time-independent variables, and 100,000 samples (detailed version of Key Example 2).

```
########### Generating a time-dependent dataset ####################
set.seed(1)
n.subject = 1e5; p.ti = 50; p.tv = 50; K = 20; n = n.subject/K; cor = 0.2;  lambda.grid = 10^seq(-10,3,0.01);
beta0.ti = NULL
beta0.tv = NULL
dat.mat0 = as.data.frame(SIM.FUN.TVC(p.ti, p.tv, n.subject, cor, beta0.ti, beta0.tv))
dat.mat0[,'strat'] = dat.mat0[,dim(dat.mat0)[2]]%%(n.subject/20)
dat.mat0 = dat.mat0[,-(dim(dat.mat0)[2]-1)]

## A time-dependent Cox model, without stratification of baseline hazard: Surv(t0,t1,status)~V4+V5+...+V103
# Example option 1: Using 1 core for computation, dividing the data into 10 chunks, using 2 iterations of one-step estimator for update
mod = dcalasso(as.formula(paste0('Surv(t0,t1,status)~',paste(paste0('V',4:103),collapse='+'))),
               family = 'cox.ph',data = dat.mat0,
               K = 10, iter.os = 2)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)

# Example option 2: Using 2 cores for parallel computation, dividing the data into 10 chunks, using 2 iterations of one-step estimator for update
modp = dcalasso(as.formula(paste0('Surv(t0,t1,status)~',paste(paste0('V',4:103),collapse='+'))),
                family = 'cox.ph',data = dat.mat0,
                K = 10, iter.os = 2, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Compare the divide and conquer method unpenalized model fit against a standard Cox fit using the whole dataset
std = coxph(as.formula(paste0('Surv(t0,t1,status)~',paste(paste0('V',4:103),
                                                          collapse='+'))),
            data = dat.mat0)
plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, mod$coefficients.unpen)

## A time-dependent Cox model, with stratification of baseline hazard: Surv(t0,t1,status)~strata(strat)+V4+V5+...+V103
# Example option 1: Using 1 core for computation, dividing the data into 10 chunks, using 4 iterations of one-step estimator for update
mod = dcalasso(as.formula(paste0('Surv(t0,t1,status)~strata(strat)+',paste(paste0('V',4:103),collapse='+'))),
               family = 'cox.ph',data = dat.mat0,
               K = 10, iter.os = 4)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)

# Example option 2: Using 2 cores for parallel computation, dividing the data into 10 chunks, using 2 iterations of one-step estimator for update
modp = dcalasso(as.formula(paste0('Surv(t0,t1,status)~strata(strat)+',paste(paste0('V',4:103),collapse='+'))),
                family = 'cox.ph',data = dat.mat0,
                K = 10, iter.os = 4, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Compare the divide and conquer method unpenalized model fit against a standard Cox fit using the whole dataset
std = coxph(as.formula(paste0('Surv(t0,t1,status)~strata(strat)+',paste(paste0('V',4:103),
                                                                        collapse='+'))),
            data = dat.mat0)
plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, mod$coefficients.unpen)
```

4. Time-dependent dataset, when the dataset are saved in multiple files: Fitting a Cox model for a dataset with 50 time-dependent variables, 50 additional time-independent variables, and 100,000 samples.

```
########### Generating a time-dependent dataset, saving the dataset to 10 separate files ####################
set.seed(1)
n.subject = 1e5; p.ti = 50; p.tv = 50; K = 20; n = n.subject/K; cor = 0.2;  lambda.grid = 10^seq(-10,3,0.01);
beta0.ti = NULL
beta0.tv = NULL
dat.mat0 = as.data.frame(SIM.FUN.TVC(p.ti, p.tv, n.subject, cor, beta0.ti, beta0.tv))
dat.mat0[,'strat'] = dat.mat0[,dim(dat.mat0)[2]]%%(n.subject/20)
dat.mat0 = dat.mat0[,-(dim(dat.mat0)[2]-1)]
ll = split(1:dim(dat.mat0)[1], factor(1:10))
for (kk in 1: 10){
  df = dat.mat0[ll[[kk]],]
  saveRDS(df, file = paste0(dir,'dataTV',kk,'.rds'))
}


## A time-dependent Cox model, without stratification of baseline hazard: Surv(t0,t1,status)~V4+V5+...+V103
# Example option 1: Using 1 core for computation, loading dataset from 10 chunks [each time the RAM will only contain data from 1 file], using 2 iterations of one-step estimator for update
mod = dcalasso(as.formula(paste0('Surv(t0,t1,status)~',paste(paste0('V',4:103),collapse='+'))),
               family = 'cox.ph',
               data.rds = paste0(dir,'dataTV',1:10,'.rds'), iter.os = 2)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)

# Example option 2: Using 2 cores for parallel computation, loading dataset from 10 chunks [each time the RAM will only contain data from 1 file], using 2 iterations of one-step estimator for update
modp = dcalasso(as.formula(paste0('Surv(t0,t1,status)~',paste(paste0('V',4:103),collapse='+'))),
                family = 'cox.ph',
                data.rds = paste0(dir,'dataTV',1:10,'.rds'), iter.os = 2, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Compare the divide and conquer method unpenalized model fit against a standard Cox fit using the whole dataset
std = coxph(as.formula(paste0('Surv(t0,t1,status)~',paste(paste0('V',4:103),
                                                          collapse='+'))),
            data = dat.mat0)
plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, mod$coefficients.unpen)

## A time-dependent Cox model, with stratification of baseline hazard: Surv(t0,t1,status)~strata(strat)+V4+V5+...+V103
# Example option 1: Using 1 core for computation, loading dataset from 10 chunks [each time the RAM will only contain data from 1 file], using 4 iterations of one-step estimator for update
mod = dcalasso(as.formula(paste0('Surv(t0,t1,status)~strata(strat)+',paste(paste0('V',4:103),collapse='+'))),
               family = 'cox.ph',
               data.rds = paste0(dir,'dataTV',1:10,'.rds'), iter.os = 4)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)

# Example option 2: Using 2 cores for parallel computation, loading dataset from 10 chunks [each time the RAM will only contain data from 1 file], using 4 iterations of one-step estimator for update
modp = dcalasso(as.formula(paste0('Surv(t0,t1,status)~strata(strat)+',paste(paste0('V',4:103),collapse='+'))),
                family = 'cox.ph',
                data.rds = paste0(dir,'dataTV',1:10,'.rds'), K = 10, iter.os = 4, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Compare the divide and conquer method unpenalized model fit against a standard Cox fit using the whole dataset
std = coxph(as.formula(paste0('Surv(t0,t1,status)~strata(strat)+',paste(paste0('V',4:103),
                                                                        collapse='+'))),
            data = dat.mat0)
plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, mod$coefficients.unpen)
```

