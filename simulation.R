# This is an R script for simulation for the PPHIS paper.
# Internal size 1000.
# 4 external data
# 2 have the same covariate (one with different censor);
# 1 has different beta;
# 1 has different distribution.
# Continuous AFT model with Weibull baseline

rm(list=ls(all=TRUE))
suppressPackageStartupMessages({
library(survival)
library(rootSolve)
library(geepack)
library(MASS)
library(wgeesel)
library(bindata)
library(psych)
library(tidyverse)
library(foreach)
library(doParallel)
library(quadprog)
})
source("functions.R")

### main data parameter config
n = 1000 # 1000/500
n_ex1 = n_ex2 = n_ex3 = n_ex4 = 1000
n_pool_tot = c(n, n_ex1, n_ex2, n_ex3, n_ex4)
betaT = c(-0.5, 1, -1, 1, 0.5)
betaT_ex = c(-0.5, -1, -1, 1, 0.5) # different coefficient
beta_cenT = c(-0.5, -1, 1)
beta_cenT_ex = c(0.5, -1, 1) # different censoring coefficient
x_cen_ind = c(1, 3, 4)
p_ex = 0.2 # external covariate mean
shape = 1.5
shape_cen = 1
endtime = 2
time_fu = 1:endtime
time_min = 0.01
#time_fu = 1:10 / 10
x_corr = matrix(0.3, 4, 4)
diag(x_corr) = 1  
dist = "weibull"

data_ind_list = rbind(c(T, T, T, F, F),
                      c(T, F, F, T, F),
                      c(T, F, F, F, T))

cov_ind_list = rbind(c(T, T, T, T, F),
                     c(T, F, F, T, T),
                     c(T, T, T, F, T))

args = commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  fn = str_remove(args[1], ".R")
}else{
  fn = Sys.time()
}
print(fn)

# parallel computation
plist = c("survival", "rootSolve", "geepack", 
          "MASS", "wgeesel",
          "bindata", "psych", "tidyverse", "quadprog")
n.cores = detectCores()

my.cluster = parallel::makeCluster(
  n.cores, 
  type = "FORK"
)
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)

iteration = 2000
timestart = Sys.time()
result = foreach(tt = 1:iteration, .combine = comb_fun,
                 .packages = plist) %dopar% 
{
  #tt = 351
  set.seed(tt)                 
  source("functions.R")
  rslt = list()                
  
  # generate main data, (x, y_obs, delta, y_sur)
  x = cbind(1, mvrnorm(n, rep(0, 4), x_corr))
  x[, 4] = as.numeric(x[, 4] > 0)
  x[, 5] = as.numeric(x[, 5] > qnorm(0.7))
  x_cen = x[, x_cen_ind]
  
  # generate raw response
  scale = exp(x %*% betaT)
  scale_cen = exp(x_cen %*% beta_cenT)
  y_main = rweibull(n, shape, scale)
  censor = rweibull(n, shape_cen, scale_cen)
  censor = pmin(censor, endtime)
  y_obs = pmin(y_main, censor)
  delta = as.numeric(y_main <= censor)
  y_int = interval_event(y_obs, delta, cont = T)
  
  # generate external data 1, no heterogeneity
  x_ex1 = cbind(1, mvrnorm(n_ex1, rep(0, 4), x_corr))
  x_ex1[, 4] = as.numeric(x_ex1[, 4] > 0)
  x_ex1[, 5] = as.numeric(x_ex1[, 5] > qnorm(0.7))
  x_cen_ex1 = x_ex1[, x_cen_ind]
  
  # generate raw response
  scale_ex1 = exp(x_ex1 %*% betaT)
  scale_cen_ex1 = exp(x_cen_ex1 %*% beta_cenT)
  y_ex1 = rweibull(n_ex1, shape, scale_ex1)
  censor_ex1 = rweibull(n_ex1, shape_cen, scale_cen_ex1)
  censor_ex1 = pmin(censor_ex1, endtime)
  y_obs_ex1 = pmin(y_ex1, censor_ex1)
  delta_ex1 = as.numeric(y_ex1 <= censor_ex1)
  y_int_ex1 = interval_event(y_obs_ex1, delta_ex1, cont = T)
  
  # generate external data 2, different censor
  x_ex2 = cbind(1, mvrnorm(n_ex2, rep(0, 4), x_corr))
  x_ex2[, 4] = as.numeric(x_ex2[, 4] > 0)
  x_ex2[, 5] = as.numeric(x_ex2[, 5] > qnorm(0.7))
  x_cen_ex2 = x_ex2[, x_cen_ind]
  
  # generate raw response
  scale_ex2 = exp(x_ex2 %*% betaT)
  scale_cen_ex2 = exp(x_cen_ex2 %*% beta_cenT_ex)
  y_ex2 = rweibull(n_ex2, shape, scale_ex2)
  censor_ex2 = rweibull(n_ex2, shape_cen, scale_cen_ex2)
  censor_ex2 = pmin(censor_ex2, endtime)
  y_obs_ex2 = pmin(y_ex2, censor_ex2)
  delta_ex2 = as.numeric(y_ex2 <= censor_ex2)
  y_int_ex2 = interval_event(y_obs_ex2, delta_ex2, cont = T)
  
  # generate external data 3, different covariate distribution
  x_ex3 = cbind(1, mvrnorm(n_ex3, rep(0, 4), x_corr))
  x_ex3[, 4] = as.numeric(x_ex3[, 4] > 0)
  x_ex3[, 5] = as.numeric(x_ex3[, 5] > qnorm(p_ex))
  x_cen_ex3 = x_ex3[, x_cen_ind]
  
  # generate raw response
  scale_ex3 = exp(x_ex3 %*% betaT)
  scale_cen_ex3 = exp(x_cen_ex3 %*% beta_cenT)
  y_ex3 = rweibull(n_ex3, shape, scale_ex3)
  censor_ex3 = rweibull(n_ex3, shape_cen, scale_cen_ex3)
  censor_ex3 = pmin(censor_ex3, endtime)
  y_obs_ex3 = pmin(y_ex3, censor_ex3)
  delta_ex3 = as.numeric(y_ex3 <= censor_ex3)
  y_int_ex3 = interval_event(y_obs_ex3, delta_ex3, cont = T)
  
  # generate external data 4, different beta
  x_ex4 = cbind(1, mvrnorm(n_ex4, rep(0, 4), x_corr))
  x_ex4[, 4] = as.numeric(x_ex4[, 4] > 0)
  x_ex4[, 5] = as.numeric(x_ex4[, 5] > qnorm(0.7))
  x_cen_ex4 = x_ex4[, x_cen_ind]
  
  # generate raw response
  scale_ex4 = exp(x_ex4 %*% betaT_ex)
  scale_cen_ex4 = exp(x_cen_ex4 %*% beta_cenT)
  y_ex4 = rweibull(n_ex4, shape, scale_ex4)
  censor_ex4 = rweibull(n_ex4, shape_cen, scale_cen_ex4)
  censor_ex4 = pmin(censor_ex4, endtime)
  y_obs_ex4 = pmin(y_ex4, censor_ex4)
  delta_ex4 = as.numeric(y_ex4 <= censor_ex4)
  y_int_ex4 = interval_event(y_obs_ex4, delta_ex4, cont = T)

  # OLS Cox
  fit = coxph(Surv(y_obs, delta) ~ x - 1)
  beta_null = fit$coefficients[-1] # main model parameter
  rslt$beta_null = beta_null
  # rslt$vols = diag(vcov(fit))[-1]
  
  # data borrow part begins here
  
  para_list = list()
  prop_scores_pool = rep()
  elfun_total_list = list()
  
  for(k in 1:nrow(data_ind_list)){
    
    subject_ind = get_subject_ind(n_pool_tot, data_ind_list[k, ])
    
    x_pool = rbind(x, x_ex1, x_ex2, 
                   x_ex3, x_ex4)[subject_ind, cov_ind_list[k, ] ]
    y_pool = rbind(y_int, y_int_ex1, y_int_ex2, 
                   y_int_ex3, y_int_ex4)[subject_ind, ]
    n_pool = n_pool_tot[data_ind_list[k, ] ]
    
    for(i in 1:length(n_pool)){
      row_ind = 1:n_pool[i]
      if(i > 1){
        row_ind = row_ind + sum(n_pool[1:(i-1)])
      }
      fit = survreg(Surv(time = y_pool[row_ind, 1], 
                         time2 = y_pool[row_ind, 2],
                         type = "interval2") ~ x_pool[row_ind, ] - 1, 
                    dist = dist)
      para = c(-log(fit$scale), fit$coefficients)
      if(i == 1){
        para_init = para
        para0 = para
      }else{
        para_init = c(para_init, para - para0)
      }
    }
    
    threshold = 1e-3
    tau_list = c(0.001, 0.005,
                 seq(0.01, 0.2, by = 0.005))
    penalty_weight = c(rep(0, length(para0)),
                       rep(1, length(para_init) - length(para0)))
    para_search_list = rep()
    bic_list = rep()
    
    for(tau in tau_list){
      
      if(tau == min(tau_list) ||
         max(abs(para - para_init)) >= 2 * max(abs(para_init))){
        lambda = NULL
        para = para_init
      }
      # lambda = NULL
      # para = para_init
      
      para_fit = para_find(para, tau, x_pool, y_pool, elfun_weibull,
                           lambda_init = lambda,
                           threshold = threshold, maxit = 30,
                           penalty_weight = penalty_weight,
                           n_pool = n_pool)
      para = para_fit$para
      para_search_list = cbind(para_search_list, para)
      
      #para = para_list[, 1]
      lambda = lambda_find(para, x_pool, y_pool,
                           elfun = elfun_weibull, n_pool = n_pool)
      ef = elfun_weibull(para, x_pool, y_pool, n_pool = n_pool)$elef
      # bic = 2*(sum(R0der(lambda, ZZ))+n*sum(scad(para, tau))) +
      #   log(n)*sum(para != 0)
      bic = 2 * (sum(R0der(lambda, ef))) + log(sum(n_pool)) * sum(para != 0)
      bic_list = c(bic_list, bic)
      
      if(all(para[penalty_weight != 0] == 0) || bic > 2 * bic_list[1]){
        break
      }
      
    }
    #plot(tau_list[1:length(bic_list)], bic_list)
    #para_search_list
    #bic_list
    ind = which.min(bic_list)
    rslt$tau = tau_list[ind]
    para = para_search_list[, ind]
    #rslt$para = para
    para_list[[k]] = para
    
    lambda = lambda_find(para, x_pool, y_pool, 
                         elfun = elfun_weibull, n_pool = n_pool)
    elfun_total = elfun_weibull(para, x_pool, y_pool, n_pool = n_pool)
    elfun_total_list[[k]] = elfun_total
    prop_scores = 1 / (1 + t(elfun_total$elef) %*% lambda)
    prop_scores = as.vector(prop_scores)[1:n]
    # summary(prop_scores)
    
    prop_scores_pool = cbind(prop_scores_pool, prop_scores)
  }
  
  # data borrow part ends here
  
  ###calculate en cox using different weights
  cox_list = cox_asym(y_int, NULL, x[, -1], beta_null)
  
  # Average scheme optimizing IIB
  
  phis_fit = phis_multi(y_int, NULL, x[, -1], beta_null,
                        para_list, prop_scores_pool,
                        elfun_total_list, n_pool_tot, data_ind_list,
                        cox_list, 
                        var_prop_weight = rep(1, length(beta_null))
                        )
  rslt$beta_avg = phis_fit$beta_en
  rslt$Vavg = phis_fit$Ven
  
  # Aggregation scheme
  
  phis_fit = phis_multi(y_int, NULL, x[, -1], beta_null,
                        para_list, prop_scores_pool,
                        elfun_total_list, n_pool_tot, data_ind_list,
                        cox_list, agg = T)
  rslt$beta_agg = phis_fit$beta_en
  rslt$Vagg = phis_fit$Ven
  
  # PCA weight
  
  phis_fit = phis_multi(y_int, NULL, x[, -1], beta_null,
                        para_list, prop_scores_pool,
                        elfun_total_list, n_pool_tot, data_ind_list,
                        cox_list, pca = T)
  rslt$beta_pca = phis_fit$beta_en
  rslt$Vnull = phis_fit$Vnull
  rslt$Vpca = phis_fit$Ven
  # rslt$Vnull/rslt$Ven

  return(rslt)
}
timeend = Sys.time()
timeend - timestart

parallel::stopCluster(cl = my.cluster)

fdir = paste0("result/", fn, ".RData")
save.image(fdir)

betaT_ph = betaT[-1]*(-shape)
beta_null = unname(rowMeans(result$beta_null))
beta_avg = rowMeans(result$beta_avg)
beta_agg = rowMeans(result$beta_agg)
beta_pca = rowMeans(result$beta_pca)

emp_var_null = diag(var(t(result$beta_null)))
emp_var_avg = diag(var(t(result$beta_avg)))
emp_var_agg = diag(var(t(result$beta_agg)))
emp_var_pca = diag(var(t(result$beta_pca)))

est_var_null = rowMeans(result$Vnull)
est_var_avg = rowMeans(result$Vavg)
est_var_agg = rowMeans(result$Vagg)
est_var_pca = rowMeans(result$Vpca)

cp_avg = sapply(1:iteration, function(i){
  return(result$beta_avg[, i] >= betaT_ph - qnorm(0.975) *
           sqrt(result$Vavg[, i]) &
           result$beta_avg[, i] <= betaT_ph + qnorm(0.975) *
           sqrt(result$Vavg[, i]))
})
cp_agg = sapply(1:iteration, function(i){
  return(result$beta_agg[, i] >= betaT_ph - qnorm(0.975) *
           sqrt(result$Vagg[, i]) &
           result$beta_agg[, i] <= betaT_ph + qnorm(0.975) *
           sqrt(result$Vagg[, i]))
})
cp_pca = sapply(1:iteration, function(i){
  return(result$beta_pca[, i] >= betaT_ph - qnorm(0.975) *
           sqrt(result$Vpca[, i]) &
           result$beta_pca[, i] <= betaT_ph + qnorm(0.975) *
           sqrt(result$Vpca[, i]))
})

summarytable_avg = data.frame(Bias = (beta_avg-betaT_ph)*100,
                              SD = sqrt(emp_var_avg)*100,
                              SE = sqrt(est_var_avg)*100,
                              CP = rowMeans(cp_avg)*100,
                              RE = emp_var_null/emp_var_avg,
                              IIB = (1 - mean(emp_var_avg / emp_var_null)) * 100
)

summarytable_agg = data.frame(Bias = (beta_agg-betaT_ph)*100,
                              SD = sqrt(emp_var_agg)*100,
                              SE = sqrt(est_var_agg)*100,
                              CP = rowMeans(cp_agg)*100,
                              RE = emp_var_null/emp_var_agg,
                              IIB = (1 - mean(emp_var_agg / emp_var_null)) * 100
)

summarytable_pca = data.frame(Bias = (beta_pca-betaT_ph)*100,
                              SD = sqrt(emp_var_pca)*100,
                              SE = sqrt(est_var_pca)*100,
                              CP = rowMeans(cp_pca)*100,
                              RE = emp_var_null/emp_var_pca,
                              IIB = (1 - mean(emp_var_pca / emp_var_null)) * 100
)

rownames(summarytable_avg) = paste0("beta", 1:4)
round(summarytable_avg, 2)

rownames(summarytable_agg) = paste0("beta", 1:4)
round(summarytable_agg, 2)

rownames(summarytable_pca) = paste0("beta", 1:4)
round(summarytable_pca, 2)

fdir = paste0("result/", fn, ".RData")
save.image(fdir)
