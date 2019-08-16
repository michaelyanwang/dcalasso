# Description

  Fast divide-and-conquer Cox proportional hazards model with adaptive lasso, developed by [Yan Wang](https://www.researchgate.net/profile/Yan_Wang374) \<yaw719@mail.harvard.edu\>, [Tianxi Cai](https://www.hsph.harvard.edu/tianxi-cai/), and [Chuan Hong](https://dbmi.hms.harvard.edu/person/postdoctoral-fellows/chuan-hong)

  The dcalasso package aims to fit large adaptive lasso for Cox model when n>>p, by divide-and-conquer strategy, least square approximation, and one-step estimation (divide-and-conquer adaptive lasso estimator). To compute a divide-and-conquer adaptive lasso estimate, the function first finds the estimate for the corresponding unpenalized Cox proportional hazards model and then shrinks it with the adaptive lasso-penalized penalty. This divide-and-conquer adaptive lasso penalty has variable selection consistency and asymptotic normality property as the standard adaptive lasso estimator.
  
  The package can fit both adaptive lasso for time-independent Cox proportional hazards model and adaptive lasso for time-dependent Cox proportional hazards model.  
  
  The dcalasso is specialized in computing the divide-and-conquer adaptive lasso with fast computation, with the following novel features: (1) a divide-and-conquer strategy is applied for the computation of the initial unpenalized Cox proportional hazards model by splitting the observations into <tt>K</tt> chunks and processing each separately; (2) a fast linearization is applied in the estimation of the shrinkage step further reducing the computational burden.
  
  When the interest is to estimate the Cox model without adaptive lasso shrinkage, this package is also useful and the user can simply take the unpenalized estimate result. Note that when n>>p, the computation for the adaptive lasso estimate takes a very small amount of time, i.e. the dcalasso is still computationally advantageous over the standard approach to fitting Cox model for extremely large dataset when the penalized estimate is not of interest.

# Citation

  Wang, Yan, Chuan Hong, Nathan Palmer, Qian Di, Joel Schwartz, Isaac Kohane, and Tianxi Cai. "[A Fast Divide-and-Conquer Sparse Cox Regression.](https://arxiv.org/pdf/1804.00735.pdf)" arXiv preprint arXiv:1804.00735 (2018).

# Installation

```
require(devtools)
install_github("michaelyanwang/dcalasso")
require(dcalasso)
```
# Key function
<tt> dcalasso(formula, family=cox.ph(), data = NULL, data.rds = NULL, weights, subsets, na.action, offset, lambda = 10^seq(-10,3,0.01), gamma = 1, K = 20, iter.os = 2, ncores = 1) </tt>

where <tt>formula</tt> is a formula of a Cox model (see <tt>coxph</tt>), <tt>family</tt> specifies the family of the outcome (<tt>cox.ph</tt> currently), <tt>data</tt> specifies the dataset, if <tt>data</tt> is too large to load into the memory, <tt>data.rds</tt> specifies a vector of file locations where the data are split and saved, <tt>weights, subsets, na.action, offset</tt> are the same as the corresponding arguments in <tt>coxph</tt>, <tt>lambda</tt> is the penalization parameter for adaptive lasso (see <tt>glmnet</tt>), <tt>gamma</tt> is the exponent for the penalization of the adaptive penalty (see <tt>glmnet</tt>), <tt>K</tt> is the number of split which will be overwritten by <tt>length(data.rds)</tt> if <tt>data.rds</tt> is specified, <tt>iter.os</tt> is the number of iterations for the computation of the initial unpenalized estimator (default 2; a larger <tt>iter.os</tt> will result in longer computation but an unpenalized estimator closer to the results at convergence, and <tt>ncores</tt> is the number of cores used in the computation.


# Examples
See <tt>?dcalasso</tt>

```
##### Time-independent #####
set.seed(1)
N = 1e5; p.x = 50; K = 100; n = N/K;  cor = 0.2;
bb = c(rep(0.4,4),rep(0.2,4),rep(0.1,4),rep(0.05,4))
beta0 = c(1, bb, rep(0, p.x - length(bb)))
dat.mat0 = as.data.frame(SIM.FUN(N, p.x = p.x, cor = cor, family='Cox',beta0 = beta0))
dat.mat0[,'strat'] = rep(1:20, each = N/20)

## Without strata
# unicore
mod = dcalasso(as.formula(paste0('Surv(u,delta)~',paste(paste0('V',3:52),collapse='+'))),
               family = 'cox.ph',data = dat.mat0,
               K = 10, iter.os = 2)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)
pred.link = predict(mod, newdata = dat.mat0)
pred.term = predict(mod, newdata = dat.mat0, type = 'terms')
pred.response = predict(mod, newdata = dat.mat0, type = 'response')

# parallel
modp = dcalasso(as.formula(paste0('Surv(u,delta)~',paste(paste0('V',3:52),collapse='+'))),
                family = 'cox.ph',data = dat.mat0,
                K = 10, iter.os = 4, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Standard
std = coxph(as.formula(paste0('Surv(u,delta)~',paste(paste0('V',3:52),collapse='+'))),
            data = dat.mat0)

plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, std$coefficients)


## With strata
# unicore
mod = dcalasso(as.formula(paste0('Surv(u,delta)~strata(strat)+',paste(paste0('V',3:52),collapse='+'))),
               family = 'cox.ph',data = dat.mat0,
               K = 10, iter.os = 2)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)

# parallel
modp = dcalasso(as.formula(paste0('Surv(u,delta)~strata(strat)+',paste(paste0('V',3:52),collapse='+'))),
                family = 'cox.ph',data = dat.mat0,
                K = 10, iter.os = 2, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Standard
std = coxph(as.formula(paste0('Surv(u,delta)~strata(strat)+',paste(paste0('V',3:52),collapse='+'))),
            data = dat.mat0)


plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, std$coefficients)




##### Time-independent with separate file saving #####
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

## Without strata
# unicore
mod = dcalasso(as.formula(paste0('Surv(u,delta)~',paste(paste0('V',3:52),collapse='+'))),
               family = 'cox.ph',
               data.rds = paste0(dir,'dataTI',1:10,'.rds'), iter.os = 2)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)

# parallel
modp = dcalasso(as.formula(paste0('Surv(u,delta)~',paste(paste0('V',3:52),collapse='+'))),
                family = 'cox.ph',
                data.rds = paste0(dir,'dataTI',1:10,'.rds'), iter.os = 2, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Standard
std = coxph(as.formula(paste0('Surv(u,delta)~',paste(paste0('V',3:52),collapse='+'))),
            data = dat.mat0)

plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, std$coefficients)


## With strata
# unicore
mod = dcalasso(as.formula(paste0('Surv(u,delta)~strata(strat)+',paste(paste0('V',3:52),collapse='+'))),
               family = 'cox.ph',
               data.rds = paste0(dir,'dataTI',1:10,'.rds'), K = 10, iter.os = 2)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)

# parallel
modp = dcalasso(as.formula(paste0('Surv(u,delta)~strata(strat)+',paste(paste0('V',3:52),collapse='+'))),
                family = 'cox.ph',
                data.rds = paste0(dir,'dataTI',1:10,'.rds'), K = 10, iter.os = 2, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Standard
std = coxph(as.formula(paste0('Surv(u,delta)~strata(strat)+',paste(paste0('V',3:52),collapse='+'))),
            data = dat.mat0)

plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, std$coefficients)


########### Time-dependent loading as a whole ####################
set.seed(1)
n.subject = 1e5; p.ti = 50; p.tv = 50; K = 20; n = n.subject/K; cor = 0.2;  lambda.grid = 10^seq(-10,3,0.01);
beta0.ti = NULL
beta0.tv = NULL
dat.mat0 = as.data.frame(SIM.FUN.TVC(p.ti, p.tv, n.subject, cor, beta0.ti, beta0.tv))
dat.mat0[,'strat'] = dat.mat0[,dim(dat.mat0)[2]]%%(n.subject/20)
dat.mat0 = dat.mat0[,-(dim(dat.mat0)[2]-1)]

## Without strata
# unicore
mod = dcalasso(as.formula(paste0('Surv(t0,t1,status)~',paste(paste0('V',4:103),collapse='+'))),
               family = 'cox.ph',data = dat.mat0,
               K = 10, iter.os = 2)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)

# parallel
modp = dcalasso(as.formula(paste0('Surv(t0,t1,status)~',paste(paste0('V',4:103),collapse='+'))),
                family = 'cox.ph',data = dat.mat0,
                K = 10, iter.os = 2, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Standard
std = coxph(as.formula(paste0('Surv(t0,t1,status)~',paste(paste0('V',4:103),
                                                          collapse='+'))),
            data = dat.mat0)
plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, mod$coefficients.unpen)

# With strata
# unicore
mod = dcalasso(as.formula(paste0('Surv(t0,t1,status)~strata(strat)+',paste(paste0('V',4:103),collapse='+'))),
               family = 'cox.ph',data = dat.mat0,
               K = 10, iter.os = 4)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)

# parallel
modp = dcalasso(as.formula(paste0('Surv(t0,t1,status)~strata(strat)+',paste(paste0('V',4:103),collapse='+'))),
                family = 'cox.ph',data = dat.mat0,
                K = 10, iter.os = 4, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Standard
std = coxph(as.formula(paste0('Surv(t0,t1,status)~strata(strat)+',paste(paste0('V',4:103),
                                                                        collapse='+'))),
            data = dat.mat0)
plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, mod$coefficients.unpen)




########### Time-dependent separate file saving ####################
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


## Without strata
# unicore
mod = dcalasso(as.formula(paste0('Surv(t0,t1,status)~',paste(paste0('V',4:103),collapse='+'))),
               family = 'cox.ph',
               data.rds = paste0(dir,'dataTV',1:10,'.rds'), K = 10, iter.os = 2)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)

# parallel
modp = dcalasso(as.formula(paste0('Surv(t0,t1,status)~',paste(paste0('V',4:103),collapse='+'))),
                family = 'cox.ph',
                data.rds = paste0(dir,'dataTV',1:10,'.rds'), K = 10, iter.os = 2, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Standard
std = coxph(as.formula(paste0('Surv(t0,t1,status)~',paste(paste0('V',4:103),
                                                          collapse='+'))),
            data = dat.mat0)
plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, mod$coefficients.unpen)

# With strata
# unicore
mod = dcalasso(as.formula(paste0('Surv(t0,t1,status)~strata(strat)+',paste(paste0('V',4:103),collapse='+'))),
               family = 'cox.ph',
               data.rds = paste0(dir,'dataTV',1:10,'.rds'), K = 10, iter.os = 4)
sum.mod = summary(mod)
print(sum.mod, unpen = T)
plot(mod)

# parallel
modp = dcalasso(as.formula(paste0('Surv(t0,t1,status)~strata(strat)+',paste(paste0('V',4:103),collapse='+'))),
                family = 'cox.ph',
                data.rds = paste0(dir,'dataTV',1:10,'.rds'), K = 10, iter.os = 4, ncores = 2)
sum.modp = summary(modp)
print(sum.modp, unpen = T)
plot(modp)

# Standard
std = coxph(as.formula(paste0('Surv(t0,t1,status)~strata(strat)+',paste(paste0('V',4:103),
                                                                        collapse='+'))),
            data = dat.mat0)
plot(mod$coefficients.unpen, std$coefficients)
plot(modp$coefficients.unpen, mod$coefficients.unpen)
```


