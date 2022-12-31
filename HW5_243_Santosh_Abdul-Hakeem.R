

#####################################################################################
##################              HW 5 Code             ##############################
##################                                    ##############################
#####################################################################################



# Installing Pacman package for loading other packages

install.packages("pacman")
library(pacman)
p_load(tidyverse, psych, MASS, GridExtra)

# clear all variables from the environment before code starts
rm(list = ls())



#####################################################################
#####################################################################
#####################################################################
#############                                    ####################
#############        Simulation Functions        ####################
#############                                    ####################
#####################################################################
#####################################################################
#####################################################################


#######################################################################

# Minimizes and produces the minimum lambda

####################################################################

lambda_min_fn = function(fn1, y){
  
  lambda_scores = sapply(X = sq, FUN = fn1,y = y)
  
  min_lambda = sq[which.min(lambda_scores)]
  
  return(min_lambda)
  
}


#######################################################################

# Cross Validation Function

####################################################################

CV_fn = function(lambda, y){ 
  
  H = X %*% ginv( (t(X) %*% X) + (lambda * D_mat)) %*% t(X) 
  
  y_hat = H %*% y
  
  hii = diag(H)
  
  CV_score =  mean( ( (y - y_hat) / (1 - hii) )^2 )
  
  return(CV_score)
  
}



#######################################################################

# Generalized Cross Validation Function

####################################################################

GCV_fn = function(lambda, y){  
  
  H = X %*% ginv( (t(X) %*% X) + (lambda * D_mat)) %*% t(X) 
  
  y_hat = H %*% y
  
  n = length(y)
  
  GCV_score =  mean( ( (y - y_hat)^2 / (1 - (1/n) * tr(H))^2 ) )
  
  return(GCV_score)
}



#######################################################################

# Akaike Informationn Criterion Function

####################################################################

AIC_fn = function(lambda, y){  
  
  H = X %*% ginv( (t(X) %*% X) + (lambda * D_mat)) %*% t(X) 
  
  y_hat = H %*% y
  
  n = length(y)
  
  AIC_score = log(sum( (y-y_hat)^2 )) + ( (2 * (tr(H) + 1) ) / (n - tr(H) - 2 ) )
  
  return(AIC_score)
  
}



#######################################################################

# Risk Minization Function

####################################################################

sigma_hat_fn = function(y){
  
  
  sigma_lambda = lambda_min_fn(CV_fn, y)
  
  H = X %*% ginv( (t(X) %*% X) + (sigma_lambda * D_mat) ) %*% t(X) 
  
  y_hat = H %*% y
  
  n = length(y)
  
  I = diag(diag(n))
  
  sigma_hat = sum( (y - y_hat)^2 ) / tr( I - H )
  
  return(sigma_hat)
}



RM_fn = function(lambda, y){  
  
  s_hat = sigma_hat_fn(y = y)
  
  H = X %*% ginv( (t(X) %*% X) + (lambda * D_mat) ) %*% t(X) 
  
  y_hat = H %*% y
  
  n = length(y)  
  
  RM_score = sum( (y-y_hat)^2 ) + s_hat * ( 2 * tr(H)  - n )
  
  return(RM_score)
  
}



#######################################################################

# Residual Sum of Squares Function

####################################################################

RSS_fn = function(lambda, f, y){
  
  H = X %*% ginv( (t(X) %*% X) + (lambda * D_mat) ) %*% t(X) 
  
  y_hat = H %*% y
  
  numer_r = sum( (f - y_hat)^2 )
  
  return(numer_r)
  
}



#######################################################################

# Cubic Spline Design Matrix Function

####################################################################

X_mat = function(x_is){
  
  rhs <- function(x,t){ifelse(x>t,x-t,0)}
  
  design <- function(x,t){rhs(x,t)^3}
  
  knots <- seq(min(x_is), max(x_is), length.out = 30)
  
  dm <- outer(x_is,knots,design)
  
  Design_mat <- cbind(rep(1,length(x_is)),x_is, x_is^2, x_is^3,dm)
  
  return(Design_mat)
  
}



#######################################################################

# Compute Ratio of Residuals R Function

####################################################################

compute_res_min_df = function(ii, j, fn1){
  
  res_yf = fn1(x_is, j)
  
  y = res_yf$y
  f = res_yf$f
  
  rss_vec = sapply(X = sq, FUN = RSS_fn, f = f, y = y)
  
  min_rss = min(rss_vec)
  
  CV_min = lambda_min_fn(CV_fn, y)
  GCV_min = lambda_min_fn(GCV_fn, y)
  AIC_min = lambda_min_fn(AIC_fn, y)
  RM_min = lambda_min_fn(RM_fn, y)
  
  method_mins_vec =  c(CV_min, GCV_min, AIC_min, RM_min)
  
  rss_method_vec = sapply(X = method_mins_vec, FUN = RSS_fn, f = f, y = y)
  
  r_vec = rss_method_vec / min_rss
  
  r_vec_with_indx = c(ii, r_vec)
  
  return( r_vec_with_indx )
  
  
}



#######################################################################

# Vectorization for compute_res_min_df Function

####################################################################

compute_for_all_j = function(j, x_is_j, X_j, fn1){
  
  x_is <<- x_is_j
  X <<- X_j
  
  res_min_df = 1:sim_length %>% map(~compute_res_min_df(., j, fn1)) %>% do.call(rbind, .) %>% as.data.frame()
  res_min_df[,2:5] = res_min_df[,2:5] %>% apply(., 2, log) 
  colnames(res_min_df) = c("id", "CV", "GCV", "AIC", "RM") 
  return(res_min_df)
  
}


#####################################################################################
#####################################################################################
#####################################################################################
##################                                    ##############################
##################         SIMULATIONS                ##############################
##################                                    ##############################
#####################################################################################
#####################################################################################
#####################################################################################


#######################################################################

# Noise Level

####################################################################

noise_fn = function(x,j){
  
  f = { 1.5 * dnorm({x - 0.35}/0.15) } - dnorm({x-0.8}/0.04)  
  sig = 0.02 + 0.04 *(j - 1)^2
  eps = rnorm(length(x))
  y = f + sig * eps 
  res = data.frame(y = y, f = f)
  
  return(res)
  
}

x_is = {1:200 - 0.5}/200

x_is_list = 1:6 %>% lapply(function(x){x_is})

x_is = NA

X_list = 1:6 %>% map(~X_mat(x_is_list[[.]]))

X = NA


D_mat = diag(c(rep(0,4), rep(1,30)))

sim_length = 200

sim_noise_list = list()

x_is = {1:200 - 0.5}/200

x_is_list = 1:6 %>% lapply(function(x){x_is})

x_is = NA

X_list = 1:6 %>% map(~X_mat(x_is_list[[.]]))

X = NA

j = 1
lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_noise_list[[j]] = j %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], noise_fn))

print("1 out of 6")

j = 2
lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))
sim_noise_list[[j]] = j %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], noise_fn))

print("2 out of 6")

j = 3
lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))
sim_noise_list[[j]] = j %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], noise_fn))

print("3 out of 6")

j = 4
lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))
sim_noise_list[[j]] = j %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], noise_fn))

print("4 out of 6")

j = 5
lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))
sim_noise_list[[j]] = j %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], noise_fn))

print("5 out of 6")

j = 6
lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))
sim_noise_list[[j]] = j %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], noise_fn))

#######################################################################

# Noise Level Test

####################################################################

# x_is = {1:200 - 0.5}/200

# x_is_list = 1:6 %>% lapply(function(x){x_is})

# x_is = NA

# X_list = 1:6 %>% map(~X_mat(x_is_list[[.]]))

# X = NA

# j = 6
# lambda.low = 10^-12
# lambda.upp = 5
# sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))
# sim_noise_list[[j]] = j %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], noise_fn))

#######################################################################

# Design Density

####################################################################

design_den_fn = function(x,j){
  
  f = { 1.5 * dnorm({x - 0.35}/0.15) } - dnorm({x-0.8}/0.04)
  sig = 0.1
  eps = rnorm(200)
  y = f + sig * eps
  res = data.frame(y = y, f = f)
  
  return(res)
  
}

rbeta_fn = function(j){
  res = rbeta(200, {j + 4}/5, {11 - j}/5)
  return(res)
}


x_is_list = 1:6 %>% map(~rbeta_fn(.))

x_is = NA

X_list = 1:6 %>% map(~X_mat(x_is_list[[.]]))

X = NA

D_mat = diag(c(rep(0,4), rep(1,30)))

sim_design_den_list = list()

sim_length = 200

x_is_list = 1:6 %>% map(~rbeta_fn(.))

x_is = NA

X_list = 1:6 %>% map(~X_mat(x_is_list[[.]]))

X = NA

lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_design_den_list[[1]] = 1 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], design_den_fn))

print("1 out of 6")

lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_design_den_list[[2]] = 2 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], design_den_fn))

print("2 out of 6")

lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_design_den_list[[3]] = 3 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], design_den_fn))

print("3 out of 6")

lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_design_den_list[[4]] = 4 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], design_den_fn))

print("4 out of 6")

lambda.low = 10^-12
lambda.upp = 2
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_design_den_list[[5]] = 5 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], design_den_fn))

print("5 out of 6")

lambda.low = 10^-20
lambda.upp = 10^-4
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_design_den_list[[6]] = 1:6 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], design_den_fn))

#######################################################################

# Design Density Test

####################################################################

# x_is_list = 1:6 %>% map(~rbeta_fn(.))
# 
# x_is = NA
# 
# X_list = 1:6 %>% map(~X_mat(x_is_list[[.]]))
# 
# X = NA
# 
# lambda.low = 10^-15
# lambda.upp = 10^-6
# sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))
# 
# sim_design_den_list[[3]] = 3 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], design_den_fn))
# 
# print("3 out of 6")
# 
# # lambda.low = 10^-12
# # lambda.upp = 5
# sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))
# 
# sim_design_den_list[[4]] = 4 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], design_den_fn))
# 
# print("4 out of 6")
# 
# # lambda.low = 10^-12
# # lambda.upp = 2
# sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))
# 
# sim_design_den_list[[5]] = 5 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], design_den_fn))
# 
# plot_fn(3,sim_design_den_list)


#######################################################################

# Spatial Variance

####################################################################

spatial_var_fn = function(x,j){
  
  exp_j = {9 - 4*j}/5
  f_j = sqrt(x * (1-x)) * sin({2 * pi * {1 + 2^{exp_j}}} / {x + 2^exp_j} )
  sig = 0.02 
  eps = rnorm(length(x))
  y = f_j + sig * eps
  res = data.frame(y = y, f = f_j)
  
  return(res)
  
}

x_is = {1:200 - 0.5}/200

x_is_list = 1:6 %>% lapply(function(x){x_is})

x_is = NA

X_list = 1:6 %>% map(~X_mat(x_is_list[[.]]))

X = NA


sim_spatial_var_list = list()

D_mat = diag(c(rep(0,4), rep(1,30)))



x_is = {1:200 - 0.5}/200

x_is_list = 1:6 %>% lapply(function(x){x_is})

x_is = NA

X_list = 1:6 %>% map(~X_mat(x_is_list[[.]]))

X = NA



lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_spatial_var_list[[1]] = 1 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], spatial_var_fn))

print("1 out of 6")

lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_spatial_var_list[[2]] = 2 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], spatial_var_fn))

print("2 out of 6")

lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_spatial_var_list[[3]] = 3 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], spatial_var_fn))

print("3 out of 6")

lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_spatial_var_list[[4]] = 4 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], spatial_var_fn))

print("4 out of 6")

lambda.low = 10^-12
lambda.upp = 2
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_spatial_var_list[[5]] = 5 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], spatial_var_fn))

print("5 out of 6")

lambda.low = 10^-20
lambda.upp = 10^-4
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_spatial_var_list[[6]] = 1:6 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], spatial_var_fn))

#######################################################################

# Spatial Variance Test

####################################################################

# x_is = {1:200 - 0.5}/200

# x_is_list = 1:6 %>% lapply(function(x){x_is})

# x_is = NA

# X_list = 1:6 %>% map(~X_mat(x_is_list[[.]]))

# X = NA

# lambda.low = 10^-20
# lambda.upp = 10^-4
# sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

# sim_var_list[[6]] = 1:6 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], var_fn))

#######################################################################

# Variance

####################################################################

var_fn = function(x,j){
  
  v_j = ( 0.15 *{1 + 0.4 *{2 * j - 7} * {x - 0.5} } )^2
  f = { 1.5 * dnorm({x - 0.35}/0.15) } - dnorm({x-0.8}/0.04)
  eps = rnorm(200)
  y = f + sqrt(v_j) * eps
  res = data.frame(y = y, f = f)
  
  return(res)
  
}

x_is = {1:200 - 0.5}/200

x_is_list = 1:6 %>% lapply(function(x){x_is})

x_is = NA

X_list = 1:6 %>% map(~X_mat(x_is_list[[.]]))

X = NA

sim_var_list = list()

D_mat = diag(c(rep(0,4), rep(1,30)))

x_is = {1:200 - 0.5}/200

x_is_list = 1:6 %>% lapply(function(x){x_is})

x_is = NA

X_list = 1:6 %>% map(~X_mat(x_is_list[[.]]))

X = NA


lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_var_list[[1]] = 1 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], var_fn))

print("1 out of 6")

lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_var_list[[2]] = 2 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], var_fn))

print("2 out of 6")

lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_var_list[[3]] = 3 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], var_fn))

print("3 out of 6")

lambda.low = 10^-12
lambda.upp = 5
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_var_list[[4]] = 4 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], var_fn))

print("4 out of 6")

lambda.low = 10^-12
lambda.upp = 2
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_var_list[[5]] = 5 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], var_fn))

print("5 out of 6")

lambda.low = 10^-20
lambda.upp = 10^-4
sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

sim_var_list[[6]] = 1:6 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], var_fn))

#######################################################################

# Variance Test

####################################################################

# x_is = {1:200 - 0.5}/200

# x_is_list = 1:6 %>% lapply(function(x){x_is})

# x_is = NA

# X_list = 1:6 %>% map(~X_mat(x_is_list[[.]]))

# X = NA

# lambda.low = 10^-20
# lambda.upp = 10^-4
# sq = 10^(seq(log10(lambda.low), log10(lambda.upp), length.out = 50))

# sim_var_list[[6]] = 1:6 %>% map(~compute_for_all_j(., x_is_list[[.]], X_list[[.]], var_fn))






#####################################################################################
#####################################################################################
#####################################################################################
##################                                    ##############################
##################           Plot Functions           ##############################
##################                                    ##############################
#####################################################################################
#####################################################################################
#####################################################################################



#######################################################################

# Compute True Regression Function

####################################################################

true_reg_fn = function(x_is, X, fn1, fn22, sq){
  
  D_mat = diag(c(rep(0,4), rep(1,30)))
  
  res1 =  1:6 %>% map(~fn1(j = ., x = x_is[[.]])) 
  
  trial_fn = function(lista, X, D_mat, fn22, sq){
    
    y = as.vector(lista[[1]])
    res = lambda_min_fn(fn22, y)
    
    return(res)
  }
  
  y_hat_fn = function(X, lambda, D_mat, y){
    y = y[[1]]
    hat_mat = X %*% ginv({t(X) %*% X} + {lambda * D_mat}) %*% t(X) 
    y_hat = hat_mat %*% y
    return(y_hat)
    
  }
  
  if(is.numeric(fn22)){
    
    y_hat_list = 1:6 %>% imap(~y_hat_fn(X = X[[.y]], lambda = fn22, D_mat = D_mat, y = res1[[.y]]))
    
  }else{
    
    res2 = res1 %>% imap(~trial_fn(lista = .x, X = X[[.y]], D_mat = D_mat, fn22 = fn22, sq = sq))
    y_hat_list = 1:6 %>% imap(~y_hat_fn(X = X[[.y]], lambda = res2[[.y]], D_mat = D_mat, y = res1[[.y]]))
    
  }
  
  y_list = res1 %>% map(function(x){x[[1]]})
  
  x_df = as.data.frame(x_is)
  colnames(x_df) = c("x1", "x2", "x3","x4","x5","x6")
  
  y_df = as.data.frame(y_list)
  colnames(y_df) = c("y1", "y2", "y3","y4","y5","y6")
  
  y_hat_df = as.data.frame(y_hat_list)
  colnames(y_hat_df) = c("z1", "z2", "z3","z4","z5","z6")
  
  df = cbind(x_df, y_df, y_hat_df)
  
  long_df = df %>% pivot_longer(everything(),
                                names_to = c(".value", "j"),
                                names_pattern = "(.)(.)"
  )
  
  return(long_df)
  
}

#######################################################################

# Plot True Regression FUnction

####################################################################


factor_plot_fn =  function(sim_df, true_df, y_max_list = rep(NA, 6)){
  
  p1 = plot_fn(jth_plot = 1, df_list = sim_df, y_max = y_max_list[1])
  p11 = plot_true_reg_fn(df = true_df, j_num = 1)
  p2 = plot_fn(jth_plot = 2, df_list = sim_df, y_max = y_max_list[2])
  p21 = plot_true_reg_fn(df = true_df, j_num = 2)
  p3 = plot_fn(jth_plot = 3, df_list = sim_df, y_max = y_max_list[3])
  p31 = plot_true_reg_fn(df = true_df, j_num = 3)  
  p4 = plot_fn(jth_plot = 4, df_list = sim_df, y_max = y_max_list[4])
  p41 = plot_true_reg_fn(df = true_df, j_num = 4)  
  p5 = plot_fn(jth_plot = 5, df_list = sim_df, y_max = y_max_list[5])
  p51 = plot_true_reg_fn(df = true_df, j_num = 5)  
  p6 = plot_fn(jth_plot = 6, df_list = sim_df, y_max = y_max_list[6])
  p61 = plot_true_reg_fn(df = true_df, j_num = 6)  
  grid.arrange(p11, p1, p21, p2, p31, p3, p41, p4, p51, p5, p61, p6, ncol = 4)
  
}


#######################################################################

# Plot Simulation Results FUnction

####################################################################


plot_fn = function(jth_plot, df_list, y_max = NA){
  
  long_df = pivot_longer(as.data.frame(df_list[[jth_plot]]), cols = 2:5, names_to = "Method", values_to = "logr")
  long_df$Method = factor(long_df$Method, levels = c("CV", "GCV", "AIC", "RM"))
  
  p = long_df %>% ggplot(aes(x = Method, y = logr)) +
    geom_boxplot(aes(fill = Method)) +
    theme_minimal(base_size = 18)+
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.text=element_text(size=12 ,face="bold"), legend.title=element_text(size=12 ,face="bold"), axis.text=element_text(size=12 ,face="bold"),
          axis.title=element_text(size=14,face="bold")) +
    labs( title = paste0("J = ", jth_plot), x = "Method", y = "ln(r)", fill = "Methods") +
    theme(plot.title = element_text(hjust = 0.5)) + ylim(0,y_max)
  
  p
  
}


#######################################################################

# Plot True Regression FUnction

####################################################################

plot_true_reg_fn = function(df, j_num){
  
  b = c("chartreuse1", "chartreuse", "chartreuse3", "chartreuse4", "green1",  "grey")
  
  p = df %>% filter(j == j_num) %>% ggplot(aes(x = x, y = y)) +
    geom_point(color = b[j_num]) +
    geom_line(aes(y = z), color = "red")  +
    theme_minimal(base_size = 18)+
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.text=element_text(size=12 ,face="bold"), legend.title=element_text(size=12 ,face="bold"), axis.text=element_text(size=12 ,face="bold"),
          axis.title=element_text(size=14,face="bold")) +
    labs( title = paste0("J = ", j_num), x = "x", y = "y") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}



#####################################################################################
#####################################################################################
#####################################################################################
##################                                    ##############################
##################      Generate Final Results        ##############################
##################                                    ##############################
#####################################################################################
#####################################################################################
#####################################################################################

1:6 %>% map(~plot_true_reg_fn(DD_true_df, .))

x_is = {1:200 - 0.5}/200

x_is = 1:6 %>% lapply(function(x){x_is})

X = 1:6 %>% map(~X_mat(x_is[[.]]))

NN_true_df = true_reg_fn(x_is, X, noise_fn, 10^-6, sq)
SV_true_df = true_reg_fn(x_is, X, spatial_var_fn, 10^-6, sq)
VF_true_df = true_reg_fn(x_is, X, var_fn, 10^-6, sq)

x_is = 1:6 %>% map(~rbeta_fn(.))

X = 1:6 %>% map(~X_mat(x_is[[.]]))

DD_true_df = true_reg_fn(x_is, X, design_den_fn, 10^-6, sq)

1:6 %>% map(~plot_true_reg_fn(DD_true_df, .))

factor_plot_fn(sim_design_den_list, DD_true_df)