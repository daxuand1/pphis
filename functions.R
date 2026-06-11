## log EL
R0der = function(lambda, ef)
{
  eps = 1/ncol(ef)
  lambda_h = 1 + as.vector(t(ef) %*% lambda)
  r0der = lambda_h
  ind = !is.na(lambda_h) & lambda_h >= eps
  r0der[ind] = log(lambda_h[ind])
  r0der[!ind] = log(eps)-1.5+2*lambda_h[!ind]/eps-lambda_h[!ind]^2/eps^2/2
  return(r0der)
}

##first derivative of log EL
R1der = function(lambda, ef)
{
  eps = 1/ncol(ef)
  lambda_h = 1 + as.vector(t(ef) %*% lambda)
  r1der = ifelse(lambda_h >= eps, 1/lambda_h,
                 2/eps-lambda_h/eps^2)
  return(r1der)
}

##second derivative of log EL
R2der = function(lambda, ef)
{
  eps = 1/ncol(ef)
  lambda_h = 1 + as.vector(t(ef) %*% lambda)
  r2der = ifelse(lambda_h >= eps, -1/lambda_h^2, -1/eps^2)
  return(r2der)
}

#function to find lambda, given theta, by Han(2019)
lambda_find = function(para, x, y, elfun, n_pool = NULL)
{
  elef = elfun(para, x, y, n_pool)$elef
  
  lambda_old = rep(0, nrow(elef))
  k=0
  tol=1e-6
  
  repeat{
    # step 1
    rl = colSums(R1der(lambda_old, elef) * t(elef))
    rll = elef %*% diag(R2der(lambda_old, elef)) %*% t(elef) 
    gamma=0
    
    # step 2
    repeat{
      update = as.vector(2^(-gamma) * ginv(rll)%*%rl)
      lambda = as.vector(lambda_old - update)
      
      
      if(sum(is.na(update)) > 0 |
         sum(update == Inf) > 0 |
         sum(update == -Inf) > 0){
        lambda = lambda_old
        break
      }else if(max(abs(update))<tol){
        break
      }
      
      index_1 = as.vector(1 + t(elef) %*% lambda) <= 1/ncol(elef)
      index_2 = sum(R0der(lambda, elef))<=sum(R0der(lambda_old, elef))
      
      if (sum(index_1)>0 | index_2>0){
        gamma=gamma+1
      }else{
        break
      }
    }
    
    # step 3
    if(max(abs(lambda-lambda_old))<tol){
      break
    }else{
      lambda_old = lambda
    }
    
  }
  
  return(lambda)
}

# ---- functions for AFT Weibull model ----

weibull_reg = function(t, x, para){
  #In addition to the regression coefficients,
  #there is a nuisance parameter for Weibull regression;
  #for numerical purpose, we take the logarithm.
  #exp(alpha) = 1/survreg's scale = rweibull's shape
  #exp(theta * x) = exp(lp) = rweibull's scale
  #S(t) = exp(-t^(exp(alpha))*exp(-lp*exp(alpha)))
  #h(t) = exp(alpha)*t^(exp(alpha)-1)*exp(-lp*exp(alpha)))
  
  alpha = para[1]
  theta = para[-1]
  
  shape = exp(alpha)
  lp = c(x %*% theta)
  err = log(t) - lp
  
  surv = exp(-t^shape*exp(-lp*shape))
  der_log_surv = log(surv)*shape*c(err,-x)
  der2_log_surv = log(surv)*shape*
    cbind(c(err,-x)*(shape*err+1),
          rbind(-x*(shape*err+1), shape*x%*%t(x)))
  
  hazard = shape * t^(shape-1) * exp(-shape*lp)
  der_log_hazard = c(1+shape*err, -shape*x)
  der2_log_hazard = shape*cbind(c(err, -x),
                           rbind(-t(x), 
                                 matrix(0, length(theta),
                                        length(theta))))
  return(list(surv = surv, 
              der_log_surv = der_log_surv,
              der2_log_surv = der2_log_surv,
              hazard = hazard,
              der_log_hazard = der_log_hazard,
              der2_log_hazard = der2_log_hazard))
}

survreg_weibull = function(para, x, y){

  ef = rep()
  def = rep()
  
  for (i in 1:nrow(x))
  {
    if(is.na(y[i, 2])){
      # L=S(t1)
      s = weibull_reg(y[i, 1], x[i, ], para)
      ef_i = s$der_log_surv
      def_i = s$der2_log_surv
    }else if(y[i, 1] == y[i, 2]){
      # L=f(t)
      s = weibull_reg(y[i, 1], x[i, ], para)
      ef_i = s$der_log_hazard + s$der_log_surv
      def_i = s$der2_log_hazard + s$der2_log_surv
    }else{
      # L=S(t1)-S(t2)
      s1 = weibull_reg(y[i, 1], x[i, ], para)
      s2 = weibull_reg(y[i, 2], x[i, ], para)
      ef_i = (s1$surv * s1$der_log_surv -
                s2$surv * s2$der_log_surv) /
        (s1$surv - s2$surv)
      def_i =
        ((s1$surv * s1$der_log_surv %*% t(s1$der_log_surv) +
            s1$surv * s1$der2_log_surv) -
           (s2$surv * s2$der_log_surv %*% t(s2$der_log_surv) +
              s2$surv * s2$der2_log_surv)) /
        (s1$surv - s2$surv) - ef_i %*% t(ef_i)
    }
    
    ef_i[is.na(ef_i) | is.infinite(ef_i)] = 0
    def_i[is.na(def_i)| is.infinite(def_i)] = 0

    ef = cbind(ef, ef_i)
    def = cbind(def, def_i)
  }
  
  return(list(ef = ef, def = def))
}

para_find_weibull = function(para, x, y){
  return(rowMeans(survreg_weibull(para, x, y)$ef))
}

elfun_weibull = function(para, x, y, n_pool = NULL)
{
  if(is.null(n_pool)){
    n_pool = nrow(x)
  }
  r = length(para)
  p = length(para) / length(n_pool)
  
  #full el ef
  elef = rep()
  delef = rep()
  # sum_delef = 0
  
  for (m in 1:length(n_pool)){
    n_m = n_pool[m]
    s = (m - 1) * p
    
    if(m == 1){
      res = survreg_weibull(para[1:p], x[1:n_m, ], y[1:n_m, ])
    }else{
      res = survreg_weibull(para[1:p] + para[1:p + s], 
                            x[1:n_m + sum(n_pool[1:(m - 1)]), ], 
                            y[1:n_m + sum(n_pool[1:(m - 1)]), ])
    }
    
    elef_m = matrix(0, r, n_m)
    elef_m[1:p + s, ] = res$ef
    delef_m = matrix(0, r, n_m * r)
    delef_m[1:p + s, get_ind(n_m * r, p, r, 0)] = 
      delef_m[1:p + s, get_ind(n_m * r, p, r, s)] = res$def
    
    elef = cbind(elef, elef_m)
    delef = cbind(delef, delef_m)
    # sum_delef = sum_delef + delef_m
  }
  
  return(list(elef = elef,
              delef = delef))
}

# --- functions for AFT exponential model ---

exp_reg = function(t, x, para){
  lp = c(x %*% para)
  surv = exp(-t * exp(-lp)) # S(t)
  der_log_surv = -c(t * exp(-lp) * (-x)) # derivative of log(S(t))
  der2_log_surv = -t * exp(-lp) * (-x) %*% t(-x) # second d
  
  hazard = exp(-lp) # h(t)
  der_log_hazard = -x # derivative of log(h(t))
  der2_log_hazard = matrix(0, nrow = length(para),
                            ncol = length(para)) # second d
  return(list(surv = surv, 
              der_log_surv = der_log_surv,
              der2_log_surv = der2_log_surv,
              hazard = hazard,
              der_log_hazard = der_log_hazard,
              der2_log_hazard = der2_log_hazard))
}

survreg_exp = function(para, x, y){
  
  ef = rep()
  def = rep()
  
  for (i in 1:nrow(x))
  {
    if(is.na(y[i, 2])){
      # L=S(t1)
      s = exp_reg(y[i, 1], x[i, ], para)
      ef_i = s$der_log_surv
      def_i = s$der2_log_surv
    }else if(y[i, 1] == y[i, 2]){
      # L=f(t)
      s = exp_reg(y[i, 1], x[i, ], para)
      ef_i = s$der_log_hazard + s$der_log_surv
      def_i = s$der2_log_hazard + s$der2_log_surv
    }else{
      # L=S(t1)-S(t2)
      s1 = exp_reg(y[i, 1], x[i, ], para)
      s2 = exp_reg(y[i, 2], x[i, ], para)
      ef_i = (s1$surv * s1$der_log_surv - 
                s2$surv * s2$der_log_surv) / 
        (s1$surv - s2$surv)
      def_i = ((s1$surv * s1$der_log_surv %*% t(s1$der_log_surv) + 
                  s1$surv * s1$der2_log_surv) - 
                 (s2$surv * s2$der_log_surv %*% t(s2$der_log_surv) + 
                    s2$surv * s2$der2_log_surv)) / 
        (s1$surv - s2$surv) - ef_i %*% t(ef_i)
    }
    
    ef_i[is.na(ef_i) | is.infinite(ef_i)] = 0
    def_i[is.na(def_i)| is.infinite(def_i)] = 0
    
    ef = cbind(ef, ef_i)
    def = cbind(def, def_i)
  }
  
  return(list(ef = ef, def = def))
}

para_find_exp = function(para, x=x_sur_pool, y=y_sur_pool){
  return(rowMeans(survreg_weibull(para, x, y)$ef))
}

elfun_exp = function(para, x, y, ...)
{
  r = length(para)
  p = length(para) / length(n_pool)
  
  #full el ef
  elef = rep()
  delef = rep()
  # sum_delef = 0
  
  for (m in 1:length(n_pool)){
    n_m = n_pool[m]
    s = (m - 1) * p
    
    if(m == 1){
      res = survreg_exp(para[1:p], x[1:n_m, ], y[1:n_m, ])
    }else{
      res = survreg_exp(para[1:p] + para[1:p + s], 
                            x[1:n_m + sum(n_pool[1:(m - 1)]), ], 
                            y[1:n_m + sum(n_pool[1:(m - 1)]), ])
    }
    
    elef_m = matrix(0, r, n_m)
    elef_m[1:p + s, ] = res$ef
    delef_m = matrix(0, r, n_m * r)
    delef_m[1:p + s, get_ind(n_m * r, p, r, 0)] = 
      delef_m[1:p + s, get_ind(n_m * r, p, r, s)] = res$def
    
    elef = cbind(elef, elef_m)
    delef = cbind(delef, delef_m)
    # sum_delef = sum_delef + delef_m
  }
  
  return(list(elef = elef,
              delef = delef))
}

# --- functions for AFT loglogistic model ---

loglogistic_reg = function(t, x, para){
  # In addition to the regression coefficients,
  # there is a nuisance parameter for log-logistic regression;
  # for numerical purpose, we take the logarithm.
  # exp(alpha) = scale = survreg's scale = rlogis's scale
  # lp = x * beta = linear predictor
  # z = (log(t) - lp) / scale
  # S(t) = 1 / (1 + exp(z))
  # h(t) = (1 / (scale * t)) * exp(z) / (1 + exp(z))
  
  p = length(para) 
  
  alpha = para[1]
  theta = para[-1]
  
  scale = exp(alpha)
  lp = sum(x * theta)
  err = log(t) - lp
  
  z  = err / scale
  ez = exp(z)
  g  = 1 / (1 + exp(-z))   # logistic(z)
  
  # Survival and hazard
  surv = 1 / (1 + ez)
  hazard = (1 / (scale * t)) * (ez / (1 + ez))
  
  # Useful vectors/matrices
  v = c(z, x / scale)      # (1+p) vector
  
  M = matrix(0, p, p)
  M[1, 1] = z
  M[1, -1] = x / scale
  M[-1, 1] = x / scale
  # M[-1,-1] already 0
  
  # --------------------------------------------------
  # 1) LOG SURVIVAL
  # --------------------------------------------------
  
  log_surv = -log(1 + exp(z))
  
  # Gradient
  der_log_surv = g * v
  
  # Hessian
  der2_log_surv = - g * (1 - g) * (v %*% t(v)) + g * M
  
  # --------------------------------------------------
  # 2) LOG HAZARD
  # --------------------------------------------------
  
  log_hazard = -alpha - log(t) + z - log(1 + exp(z))
  
  # Gradient
  der_log_hazard = c(-1, rep(0, p - 1)) - (1 - g) * v
  
  # Hessian
  der2_log_hazard = - g * (1 - g) * (v %*% t(v)) - (1 - g) * M
  
  return(list(surv = surv, 
              der_log_surv = der_log_surv,
              der2_log_surv = der2_log_surv,
              hazard = hazard,
              der_log_hazard = der_log_hazard,
              der2_log_hazard = der2_log_hazard))
}

survreg_loglogistic = function(para, x, y){
  
  ef = rep()
  def = rep()
  
  for (i in 1:nrow(x))
  {
    if(is.na(y[i, 2])){
      # L=S(t1)
      s = loglogistic_reg(y[i, 1], x[i, ], para)
      ef_i = s$der_log_surv
      def_i = s$der2_log_surv
    }else if(y[i, 1] == y[i, 2]){
      # L=f(t)
      s = loglogistic_reg(y[i, 1], x[i, ], para)
      ef_i = s$der_log_hazard + s$der_log_surv
      def_i = s$der2_log_hazard + s$der2_log_surv
    }else{
      # L=S(t1)-S(t2)
      s1 = loglogistic_reg(y[i, 1], x[i, ], para)
      s2 = loglogistic_reg(y[i, 2], x[i, ], para)
      ef_i = (s1$surv * s1$der_log_surv -
                s2$surv * s2$der_log_surv) /
        (s1$surv - s2$surv)
      def_i =
        ((s1$surv * s1$der_log_surv %*% t(s1$der_log_surv) +
            s1$surv * s1$der2_log_surv) -
           (s2$surv * s2$der_log_surv %*% t(s2$der_log_surv) +
              s2$surv * s2$der2_log_surv)) /
        (s1$surv - s2$surv) - ef_i %*% t(ef_i)
    }
    
    ef_i[is.na(ef_i) | is.infinite(ef_i)] = 0
    def_i[is.na(def_i)| is.infinite(def_i)] = 0
    
    ef = cbind(ef, ef_i)
    def = cbind(def, def_i)
  }
  
  return(list(ef = ef, def = def))
}

para_find_loglogistic = function(para, x, y){
  return(rowMeans(survreg_loglogistic(para, x, y)$ef))
}

elfun_loglogistic = function(para, x, y, n_pool = NULL)
{
  if(is.null(n_pool)){
    n_pool = nrow(x)
  }
  r = length(para)
  p = length(para) / length(n_pool)
  
  #full el ef
  elef = rep()
  delef = rep()
  # sum_delef = 0
  
  for (m in 1:length(n_pool)){
    n_m = n_pool[m]
    s = (m - 1) * p
    
    if(m == 1){
      res = survreg_loglogistic(para[1:p], x[1:n_m, ], y[1:n_m, ])
    }else{
      res = survreg_loglogistic(para[1:p] + para[1:p + s], 
                                x[1:n_m + sum(n_pool[1:(m - 1)]), ], 
                                y[1:n_m + sum(n_pool[1:(m - 1)]), ])
    }
    
    elef_m = matrix(0, r, n_m)
    elef_m[1:p + s, ] = res$ef
    delef_m = matrix(0, r, n_m * r)
    delef_m[1:p + s, get_ind(n_m * r, p, r, 0)] = 
      delef_m[1:p + s, get_ind(n_m * r, p, r, s)] = res$def
    
    elef = cbind(elef, elef_m)
    delef = cbind(delef, delef_m)
    # sum_delef = sum_delef + delef_m
  }
  
  return(list(elef = elef,
              delef = delef))
}

# --- penalty functions ---

scad = function(para, t = tau, a = 3.7){
  ifelse(abs(para)<=t, t*abs(para),
         ifelse(abs(para)<=a*t, (-para^2+2*a*t*abs(para)-t^2)/2/(a-1),
                (a+1)*t^2/2))
}

scad_der = function(para, t = tau, a = 3.7){
  t * (as.numeric(abs(para) <= t) + 
         pmax(a * t - abs(para), 0) * 
         as.numeric(abs(para) > t) / (a - 1) / t)
}

para_find = function(para_init, tau, x, y, elfun, n_pool = NULL,
                     lambda_init = NULL,
                     threshold = 1e-3, maxit = 30,
                     bound = 30, penalty_weight = NULL, ...){
  #lambda_init=NULL;tau = 0.05;x=x_pool1;y=y_pool1;
  #elfun=elfun_weibull;maxit = 30;bound=30;n_pool = n_pool1
  
  para_old = para_init
  if(is.null(lambda_init)){
    lambda = lambda_find(para_old, x, y, elfun, n_pool)
  }else{
    lambda = lambda_init
  }
  k = 0
  n = nrow(x)
  if(is.null(penalty_weight)){
    penalty_weight = rep(1, length(para_init))
  }
  
  repeat{
    
    total = elfun(para_old, x, y, n_pool = n_pool)
    ef = total$elef
    def = total$delef
    m_old = sum(R0der(lambda, ef)) + 
      sum(scad(para_old, tau) * penalty_weight) * n
    
    n = ncol(ef)
    p = nrow(ef)    
    non_zero_ind = para_old != 0 # only update non-zero coefficients
    def[, rep(!non_zero_ind, n)] = 0
    
    penalty_factor = penalty_weight * 
      ifelse(para_old == 0, 0, scad_der(para_old, tau)/abs(para_old))
    
    # calculate gradient and hessian following Han(2019)
    # h - a pxp matrix
    scaler = R1der(lambda, ef)
    h = def
    dim(h) = c(p^2, n)
    h = t(t(h) * scaler)
    dim(h) = c(p, p, n)
    h1 = apply(h, 1:2, sum)
    m_der = t(h1) %*% lambda + penalty_factor * para_old * n
    
    m_der2 = - t(h1) %*% 
      ginv(ef %*% diag(R2der(lambda,ef)) %*% t(ef)) %*% (h1) + 
      diag(penalty_factor) * n
    
    if(sum(is.na(m_der2)) > 0|
       sum(m_der2 == Inf) > 0|
       sum(m_der2 == -Inf) > 0){
      direction = rep(0, length(m_der))
    }else{
      direction = ginv(m_der2) %*% m_der
    }
    
    for(sigma in 0:5){
      step_length = 2^(-sigma)
      
      para_temp = as.vector(para_old - step_length * direction)
      ef_temp = elfun(para_temp, x, y, n_pool = n_pool)$elef
      lambda = lambda_find(para_temp, x, y, elfun, n_pool = n_pool)
      m_temp = sum(R0der(lambda, ef_temp)) +
        sum(scad(para_old, tau) + n *
              penalty_factor * (para_temp^2 - para_old^2)/2)
      if(!is.na(m_temp) & m_temp <= m_old){
        para = para_temp
        para = ifelse(abs(para) < threshold & penalty_weight != 0, 0, para)
        break
      }
    }
    
    if(max(abs(para - para_old)) < threshold){
      convergence = T
      break
    }else if( k > maxit | max(abs(para)) > bound){
      convergence = F
      break
    }else{
      para_old = para
      k = k+1
      lambda = lambda_find(para_old, x, y, elfun, n_pool = n_pool)
    }
    
  }
  
  #para = ifelse(abs(para) < 1e-3, 0, para)
  return(list(para = para,
              converge = convergence)
        )
  #return(para = para)
}

# ---- Cox model related functions ----
plee_en = function(beta, x, y, delta = NULL, prop_scores)
{
  n = nrow(x)
  if(all(x[, 1] == 1)){
    x = x[, -1]
  }
  if(is.null(delta)){
    delta = !is.na(y[, 2])
    y = y[, 1]
  }
  sf = 0
  
  for(i in 1:n){
    ind = y >= y[i]
    
    s1 = prop_scores[ind] * x[ind, ] * c(exp(x[ind, ] %*% beta))
    
    if(sum(ind) > 1){
      s1 = c(colSums(s1))
    }
    
    s2 = sum(c(exp(x[ind, ] %*% beta)) * prop_scores[ind])
    
    sf = sf + prop_scores[i] * delta[i] * (x[i, ] - s1/s2)
  }
  return(sf)
}

cum_base_hazard = function(y, delta = NULL, x, beta, time){
  if(is.null(delta)){
    delta = !is.na(y[, 2])
    y = y[, 1]
  }
  n = nrow(x)
  
  Lambdahat = 0
  for (i in 1:n)
  {
    s_n1 = sum(exp(x %*% beta)[y >= y[i]])
    Lambdahat_i = delta[i] * ifelse(y[i] <=  time, 1, 0) / s_n1
    Lambdahat = Lambdahat + Lambdahat_i
  }
  Lambdahat
}


#calculate asymptotics for naive cox
cox_asym = function(y, delta = NULL, x, beta){
  #y = y_obs; x = x[, -1]; beta = beta_en;
  if(all(x[, 1] == 1)){
    x = x[, -1]
  }
  if(is.null(delta)){
    delta = !is.na(y[, 2])
    y = y[, 1]
  }
  
  n = nrow(x)
  p = length(beta)
  
  # reorder data based on observed times
  match_id = match(1:n, order(y))
  
  delta_sort = delta[order(y)]
  x_sort = x[order(y), ]
  y_sort = y[order(y)]
  
  # calculate score function and derivatives, named G and D in Breslow (2015).
  # when correctly specified, not need to calculate D, because var(G) = D,
  # then the variance is var(G)^(-1), as in Aim 2;
  # however, here we want to be robust.
  
  time_old = 0
  cum_base_hard_seq_i_old = 0
  time_used_old = 0
  cum_sum_term_i_old = 0
  cox_ef = rep() # score
  der = 0 # derivative
  
  for (i in 1:n){
    time = y_sort[i]
    
    s_n1_seq_i = sum(exp(x_sort %*% beta)[y_sort >=  time]) / n
    
    s_nz_seq_i = 
      (x_sort * as.vector(exp(x_sort %*% beta)))[y_sort >= time, ]
    if(sum(y_sort >= time) > 1){
      s_nz_seq_i = colSums(s_nz_seq_i)
    }
    s_nz_seq_i = s_nz_seq_i / n
    
    s_nz2_seq_i = apply(
      matrix(x_sort[y_sort >= time, ], ncol = ncol(x_sort)), 
      1, function(x){
        x %*% t(x) * c(exp(x %*% beta))
      })
    dim(s_nz2_seq_i) = c(p, p, sum(y_sort >= time))
    s_nz2_seq_i = apply(s_nz2_seq_i, 1:2, sum) / n
    
    # first part of score
    cox_ef_1_i = delta_sort[i]*
      (x_sort[i, ] - s_nz_seq_i/s_n1_seq_i)
    
    # second part of score
    cum_base_hard_seq_i = 
      cum_base_hazard(y_sort, delta_sort, x_sort, beta, time)
    cum_base_hard_seq_diff = 
      cum_base_hard_seq_i - cum_base_hard_seq_i_old
    
    cum_sum_term_i = cum_sum_term_i_old + 
      (-s_nz_seq_i / s_n1_seq_i) * cum_base_hard_seq_diff
    
    cox_ef_2_i = as.vector(exp(x_sort[i, ] %*% beta))*
      (cum_sum_term_i + x_sort[i, ] * cum_base_hard_seq_i)
    
    # score
    cox_ef_i = cox_ef_1_i - cox_ef_2_i
    cox_ef = cbind(cox_ef, cox_ef_i)
    
    # derivative
    der_i = (s_nz2_seq_i - s_nz_seq_i %*% t(s_nz_seq_i) / s_n1_seq_i) *
      cum_base_hard_seq_diff
    
    #update each old term
    time_old = time
    cum_base_hard_seq_i_old = cum_base_hard_seq_i
    cum_sum_term_i_old = cum_sum_term_i
    der = der + der_i
  }
  
  #recover their order
  cox_ef = cox_ef[, match_id]
  
  # calculate statistics in V
  Sigma = cox_ef %*% t(cox_ef) / n
  # der = Sigma
  Vnull = ginv(der) %*% Sigma %*% t(ginv(der))
  
  return(list(Gamma = der, # derivative
              Sigma = Sigma, # second moment
              ef = cox_ef, # estimating function
              Vnull = Vnull / n # null variance matrix
  ))
}

#calculate asymptotics for en cox
phis_asym = function(y, delta = NULL, x, beta, 
                     para, elfun_total, n_pool){
  #y = y_obs; x = x[, -1]; beta = beta_en;
  if(is.null(delta)){
    delta = !is.na(y[, 2])
    y = y[, 1]
  }
  
  n = n_pool[1]
  p = length(beta)
  
  # reorder data based on observed times
  match_id = match(1:n, order(y))
  
  delta_sort = delta[order(y)]
  x_sort = x[order(y), ]
  y_sort = y[order(y)]
  
  # calculate score function and derivatives, named G and D in Breslow (2015).
  # when correctly specified, not need to calculate D, because var(G) = D,
  # then the variance is var(G)^(-1), as in Aim 2;
  # however, here we want to be robust.
  
  time_old = 0
  cum_base_hard_seq_i_old = 0
  time_used_old = 0
  cum_sum_term_i_old = 0
  cox_ef = rep() # score
  der = 0 # derivative
  
  for (i in 1:n){
    time = y_sort[i]
    
    s_n1_seq_i = sum(exp(x_sort %*% beta)[y_sort >=  time]) / n
    
    s_nz_seq_i = 
      (x_sort * as.vector(exp(x_sort %*% beta)))[y_sort >= time, ]
    if(sum(y_sort >= time) > 1){
      s_nz_seq_i = colSums(s_nz_seq_i)
    }
    s_nz_seq_i = s_nz_seq_i / n
    
    s_nz2_seq_i = apply(
      matrix(x_sort[y_sort >= time, ], ncol = ncol(x_sort)), 
      1, function(x){
        x %*% t(x) * c(exp(x %*% beta))
      })
    dim(s_nz2_seq_i) = c(p, p, sum(y_sort >= time))
    s_nz2_seq_i = apply(s_nz2_seq_i, 1:2, sum) / n
    
    # first part of score
    cox_ef_1_i = delta_sort[i]*
      (x_sort[i, ] - s_nz_seq_i/s_n1_seq_i)
    
    # second part of score
    cum_base_hard_seq_i = 
      cum_base_hazard(y_sort, delta_sort, x_sort, beta, time)
    cum_base_hard_seq_diff = 
      cum_base_hard_seq_i - cum_base_hard_seq_i_old
    
    cum_sum_term_i = cum_sum_term_i_old + 
      (-s_nz_seq_i / s_n1_seq_i) * cum_base_hard_seq_diff
    
    cox_ef_2_i = as.vector(exp(x_sort[i, ] %*% beta))*
      (cum_sum_term_i + x_sort[i, ] * cum_base_hard_seq_i)
    
    # score
    cox_ef_i = cox_ef_1_i - cox_ef_2_i
    cox_ef = cbind(cox_ef, cox_ef_i)
    
    # derivative
    der_i = (s_nz2_seq_i - s_nz_seq_i %*% t(s_nz_seq_i) / s_n1_seq_i) *
      cum_base_hard_seq_diff
    
    #update each old term
    time_old = time
    cum_base_hard_seq_i_old = cum_base_hard_seq_i
    cum_sum_term_i_old = cum_sum_term_i
    der = der + der_i
  }
  
  #recover their order
  cox_ef = cox_ef[, match_id]
  
  # calculate statistics in V
  Sigma = cox_ef %*% t(cox_ef) / n
  # der = Sigma
  Vnull = ginv(der) %*% Sigma %*% t(ginv(der))
  
  # calculate improved variance part
  elef = elfun_total$elef
  delef = elfun_total$delef
  
  ntot = ncol(elef) # total size of working model
  r = nrow(elef) # dim of working model
  m = length(n_pool) # number of dataset
  q = r / m # dim of working model parameter
  
  if(sum(para == 0) == 0){
    return(list(beta = beta,
                Vnull = Vnull / n,
                Ven = Vnull / n))
  }
  
  H = diag(r)
  H = matrix(H[para == 0, ], nrow = sum(para == 0))
  
  A = elef %*% t(elef) / ntot
  B = delef
  dim(B) = c(r, r, ntot)
  B = apply(B, 1:2, mean)
  C = ginv(t(B) %*% ginv(A) %*% B)
  
  S = ginv(A) - ginv(A) %*% B %*% C %*% t(B) %*% ginv(A)
  P = ginv(A) %*% B %*% C %*% t(H) %*% ginv(H %*% C %*% t(H)) %*%
    H %*% C %*% t(B) %*% ginv(A)
  
  Lambda = cox_ef %*% t(elef[, 1:n]) / n
  # Lambda = cbind(Lambda, matrix(0, p, r ))
  pi = n / ntot
  Ven = ginv(der) %*% 
    (Sigma - pi * Lambda %*% (S + P) %*% t(Lambda)) %*% t(ginv(der))
  
  return(list(beta = beta,
              Vnull = Vnull / n,
              Ven = Ven / n))
}


phis_multi = function(y, delta = NULL, x, beta,
                      para_list, prop_scores_pool,
                      elfun_total_list, n_pool_tot, data_ind_list,
                      cox_list = NULL, agg = F,
                      fixed_weight = NULL, linear_weight = NULL,
                      var_prop_weight = NULL, prediction = F, pca = F){
  if(is.null(delta)){
    delta = !is.na(y[, 2])
    y = y[, 1]
  }
  
  k = length(para_list) # number of models
  p = length(beta) # length of cox parameter
  n = nrow(x) # internal sample size
  
  if(is.null(cox_list)){
    cox_list = cox_asym(y, delta, x, beta)
  }
  
  # extract statistics from cox_list
  cox_ef = cox_list$ef
  Sigma = cox_list$Sigma
  Gamma = cox_list$Gamma
  Vnull = cox_list$Vnull * n
  
  # calculate statistics for each working model
  matrix_array = array(dim = c(p, p, k, k))
  Lambda_list = list()
  P_list = list()
  
  for(i in 1:k){
    para = para_list[[i]]
    elfun_total = elfun_total_list[[i]]
    n_pool = n_pool_tot[data_ind_list[i, ]]
    
    # calculate improved variance part
    elef = elfun_total$elef
    delef = elfun_total$delef
    
    ntot = ncol(elef) # total size of working model
    r = nrow(elef) # dim of working model
    
    H = diag(r)
    H = matrix(H[para == 0, ], nrow = sum(para == 0))
    
    A = elef %*% t(elef) / ntot
    B = delef
    dim(B) = c(r, r, ntot)
    B = apply(B, 1:2, mean)
    C = ginv(t(B) %*% ginv(A) %*% B)
    
    # S = ginv(A) - ginv(A) %*% B %*% C %*% t(B) %*% ginv(A)
    P = ginv(A) %*% B %*% C %*% t(H) %*% ginv(H %*% C %*% t(H)) %*%
      H %*% C %*% t(B) %*% ginv(A)
    
    Lambda = cox_ef %*% t(elef[, 1:n]) / n
    
    P_list[[i]] = P
    Lambda_list[[i]] = Lambda
  }
  
  for(i in 1:k){
    n_pool_i = n_pool_tot[data_ind_list[i, ]]
    
    matrix_array[, , i, i] = 
      n / sum(n_pool_i) * ginv(Gamma) %*% Lambda_list[[i]] %*% 
      P_list[[i]] %*% t(Lambda_list[[i]]) %*% t(ginv(Gamma))
    
    if(i < k){
      for(j in (i+1):k){
        joint_ind_ij = data_ind_list[i, ] & data_ind_list[j, ]
        
        if(all(joint_ind_ij == F)){
          Phi_ij = matrix(0, length(para_list[[i]]), length(para_list[[j]]))
        }else{
          n_pool_j = n_pool_tot[data_ind_list[j, ]]
          joint_ind_i = joint_ind_ij[data_ind_list[i, ]]
          joint_ind_j = joint_ind_ij[data_ind_list[j, ]]
          
          ef_i = (elfun_total_list[[i]])$elef
          ef_i = ef_i[, get_subject_ind(n_pool_i, joint_ind_i)]
          
          ef_j = (elfun_total_list[[j]])$elef
          ef_j = ef_j[, get_subject_ind(n_pool_j, joint_ind_j)]
          
          Phi_ij = ef_i %*% t(ef_j)
        }
        
        
        matrix_array[, , i, j] = n / sum(n_pool_i) / sum(n_pool_j) * 
          ginv(Gamma) %*% Lambda_list[[i]] %*% P_list[[i]] %*% Phi_ij %*%
          P_list[[j]] %*% t(Lambda_list[[j]]) %*% t(ginv(Gamma))
        
        matrix_array[, , j, i] = t(matrix_array[, , i, j])
      }
    }
  }
  
  # aggregation scheme
  if(agg == T){
    prop_scores = apply(prop_scores_pool, 1, prod)
    
    beta_en = multiroot(f = plee_en,
                        start = as.vector(beta),
                        x = x, y = y, delta = delta,
                        prop_scores = prop_scores)$root
    Ven = Vnull
    
    for(i in 1:k){
      for(j in 1:k){
        if(i == j){
          Ven = Ven - matrix_array[, , i, i]
        }else{
          Ven = Ven +  matrix_array[, , i, j]
        }
      }
    }
  }
  
  # average scheme with fixed weight
  if(!is.null(fixed_weight)){
    prop_scores = as.vector(prop_scores_pool %*% fixed_weight)
    
    beta_en = multiroot(f = plee_en,
                        start = as.vector(beta),
                        x = x, y = y, delta = delta,
                        prop_scores = prop_scores)$root
    Ven = Vnull
    
    for(i in 1:k){
      for(j in 1:k){
        if(i == j){
          Ven = Ven + 
            (fixed_weight[i]^2 - 2 * fixed_weight[i]) * matrix_array[, , i, i]
        }else{
          Ven = Ven + fixed_weight[i] * fixed_weight[j] * matrix_array[, , i, j]
        }
      }
    }
  }
  
  # average scheme optimizing linear combination
  if(!is.null(linear_weight)){
    
    # apply quadprog
    Dmat = matrix(0, k, k)
    dvec = rep(0, k)
    Amat = matrix(1, k, 1)
    bvec = c(1)
    
    for(i in 1:k){
      for(j in i:k){
        if(i == j){
          Dmat[i, i] = dvec[i] = 
            t(linear_weight) %*% matrix_array[, , i, i] %*% linear_weight
        }else{
          Dmat[i, j] = 
            t(linear_weight) %*% matrix_array[, , i, j] %*% linear_weight
        }
      }
    }
    
    weight = solve.QP(Dmat, dvec, Amat, bvec, meq = 1)$solution
    
    prop_scores = as.vector(prop_scores_pool %*% weight)
    
    beta_en = multiroot(f = plee_en,
                        start = as.vector(beta),
                        x = x, y = y, delta = delta,
                        prop_scores = prop_scores)$root
    Ven = Vnull
    
    for(i in 1:k){
      for(j in 1:k){
        if(i == j){
          Ven = Ven + 
            (fixed_weight[i]^2 - 2 * weight[i]) * matrix_array[, , i, i]
        }else{
          Ven = Ven + fixed_weight[i] * weight[j] * matrix_array[, , i, j]
        }
      }
    }
    
    
    
  }
  
  # average scheme optimizing sum of variance reduction proportion
  if(!is.null(var_prop_weight)){
    
    # apply quadprog
    Dmat = matrix(0, k, k)
    dvec = rep(0, k)
    Amat = matrix(1, k, 1)
    bvec = c(1)
    
    for(i in 1:k){
      for(j in i:k){
        if(i == j){
          Dmat[i, i] = dvec[i] = 
            sum(diag(matrix_array[, , i, i] / Vnull) * var_prop_weight)
        }else{
          Dmat[i, j] = 
            sum(diag(matrix_array[, , i, j] / Vnull) * var_prop_weight)
        }
      }
    }
    
    weight = solve.QP(Dmat, dvec, Amat, bvec, meq = 1)$solution
    
    prop_scores = as.vector(prop_scores_pool %*% weight)
    
    beta_en = multiroot(f = plee_en,
                        start = as.vector(beta),
                        x = x, y = y, delta = delta,
                        prop_scores = prop_scores)$root
    Ven = Vnull
    
    for(i in 1:k){
      for(j in 1:k){
        if(i == j){
          Ven = Ven + 
            (weight[i]^2 - 2 * weight[i]) * matrix_array[, , i, i]
        }else{
          Ven = Ven + weight[i] * weight[j] * matrix_array[, , i, j]
        }
      }
    }
    
  }
  
  # average scheme optimizing prediction variance
  if(prediction == T){
    
    # apply quadprog
    Dmat = matrix(0, k, k)
    dvec = rep(0, k)
    Amat = matrix(1, k, 1)
    bvec = c(1)
    
    for(i in 1:k){
      for(j in i:k){
        if(i == j){
          Dmat[i, i] = dvec[i] = 
            mean(apply(x, 1, function(x){
              t(x) %*% matrix_array[, , i, i] %*% x }))
        }else{
          Dmat[i, j] = 
            mean(apply(x, 1, function(x){
              t(x) %*% matrix_array[, , i, j] %*% x }))
        }
      }
    }
    
    weight = solve.QP(Dmat, dvec, Amat, bvec, meq = 1)$solution
    
    prop_scores = as.vector(prop_scores_pool %*% weight)
    
    beta_en = multiroot(f = plee_en,
                        start = as.vector(beta),
                        x = x, y = y, delta = delta,
                        prop_scores = prop_scores)$root
    Ven = Vnull
    
    for(i in 1:k){
      for(j in 1:k){
        if(i == j){
          Ven = Ven + 
            (weight[i]^2 - 2 * weight[i]) * matrix_array[, , i, i]
        }else{
          Ven = Ven + weight[i] * weight[j] * matrix_array[, , i, j]
        }
      }
    }
    
  }
  
  # PCA weights
  if(pca == T){
    
    ef = rep()
    q_pool = sapply(para_list, length) # length of working parameter
    R = matrix(0, sum(q_pool), sum(q_pool))
    
    for(i in 1:k){
      
      para = para_list[[i]]
      elfun_total = elfun_total_list[[i]]
      n_pool = n_pool_tot[data_ind_list[i, ]]
      
      # calculate improved variance part
      ef_i = elfun_total$elef
      def_i = elfun_total$delef
      
      ef_i = enlarge_ef(ef_i, n_pool_tot, data_ind_list[i, ]) 
      ef = rbind(ef, ef_i)
      
      ntot = ncol(ef_i) # total size of working model
      r = nrow(ef_i) # dim of working model
      
      H = diag(r)
      H = matrix(H[para == 0, ], nrow = sum(para == 0))
      
      A = ef_i %*% t(ef_i) / ntot
      B = def_i
      dim(B) = c(r, r, sum(n_pool))
      B = apply(B, 1:2, sum) / ntot
      C = ginv(B) %*% A  %*% ginv(t(B))
      
      P = ginv(A) %*% B %*% C %*% t(H) %*% ginv(H %*% C %*% t(H)) %*%
        H %*% C %*% t(B) %*% ginv(A)
      
      working_ind = rep(F, k)
      working_ind[i] = T
      working_ind = get_subject_ind(q_pool, working_ind)
      R[working_ind, working_ind] = A %*% P_list[[i]]
    }
    
    pca_fit = summary(prcomp(t(R %*% ef)))
    num_pc = sum(sapply(para_list, function(x){sum(x == 0)}))
    
    ef = t(as.matrix(pca_fit$rotation[, 1:num_pc])) %*% R %*% ef
    if(num_pc == 1){
      lambda = (pca_fit$sdev[1])^(-2) * rowMeans(ef)
    }else{
      lambda = diag((pca_fit$sdev[1:num_pc])^(-2)) %*% rowMeans(ef)
    }
    
    prop_scores = 1 - as.vector(t(ef[, 1:n]) %*% lambda)
    
    beta_en = multiroot(f = plee_en,
                        start = as.vector(beta),
                        x = x, y = y, delta = delta,
                        prop_scores = prop_scores)$root
    
    Psi = cox_ef %*% t(ef[, 1:n]) / ntot
    
    if(num_pc == 1){
      Ven = Vnull - ntot / n * ginv(Gamma) %*% 
        ( (pca_fit$sdev[1])^(-2) * Psi %*% t(Psi) ) %*% 
        t(ginv(Gamma))
    }else{
      Ven = Vnull - ntot / n * ginv(Gamma) %*% 
        ( Psi %*% diag((pca_fit$sdev[1:num_pc])^(-2)) %*% t(Psi) ) %*% 
        t(ginv(Gamma))
    }
    
    
  }
  
  if(agg == T | !is.null(fixed_weight) | pca == T){
    return(list(beta_en = beta_en,
                Vnull = diag(Vnull) / n,
                Ven = diag(Ven) / n))
  }else{
    return(list(beta_en = beta_en,
                Vnull = diag(Vnull) / n,
                Ven = diag(Ven) / n,
                weight = weight))
  }
  
}


# --- other functions ---

interval_event = function(y, delta, time_list = NULL, cont = F){
  # return a time interval; censored at the beginning.
  
  if(cont == T){
    res = cbind(y, ifelse(delta == 1, y, NA))
  }else{
    res = rep()
    
    for(i in 1:length(y)){
      if(y[i] <= min(time_list)){
        
        if(delta[i] == 1){
          
          time_int = c(time_min, min(time_list))
          
        }else{
          
          time_int = c(time_min, NA)
          
        }
        
      }else if(y[i] > max(time_list)){
        
        time_int = c(max(time_list), NA)
        
      }else if(delta[i] == 1){
        
        time_int = c(time_list[max(which(y[i]>time_list))],
                     time_list[min(which(y[i]<=time_list))])
        
      }else{
        
        time_int = c(time_list[max(which(y[i]>time_list))],
                     NA)
        
      }
      res = rbind(res, time_int)
    } 
  }
  
  return(res)
}

comb_fun = function(list1, list2){
  mapply(cbind, list1, list2, SIMPLIFY=F)
}

get_ind = function(n, p, q, s = 0) {
  # generate the indices for p elements from each group of q with shift s
  ind = sapply(seq(1, n, q) + s, 
                function(x) {
                  seq(x, x + p - 1)
                  }
                )
  as.vector(ind)
}

get_subject_ind = function(n_pool, ind_pool){
  # return an index vector to identify the subject index
  ind = rep(F, sum(n_pool))
  for(i in 1:length(ind_pool)){
    if(ind_pool[i] == T){
      if(i == 1){
        ind[1:(n_pool[1])] = T
      }else{
        ind[1:(n_pool[i]) + sum(n_pool[1:(i-1)])] = T
      }
    }
  }
  return(ind)
}

enlarge_ef = function(ef, n_pool, ind_pool){
  # return an enlarged estimating function pool, with 0 for subjects not in
  # this working model.
  ef_full = matrix(0, nrow(ef), sum(n_pool))
  ef_full[, get_subject_ind(n_pool, ind_pool)] = ef
  return(ef_full)
}


