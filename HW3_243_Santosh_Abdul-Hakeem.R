# Installing Pacman package for loading other packages

# install.packages("pacman")
library(pacman)

# Loading the needed packages

p_load(tidyverse, rbenchmark)


# clear all variables from the environment before code starts

rm(list = ls())


#####################################################################################
#####################################################################################
#####################################################################################

# Question 4

#####################################################################################
#####################################################################################
#####################################################################################



#######################################################################

# Inverse Transform Algorithm function

####################################################################


inverse_transform_fn = function(g){
  
  u = runif(1)
  res = g(u)
  return(res)
  
}


#######################################################################

# Inverse of f CDF function

####################################################################

f_cdf_inverse_fn = function(x){
  res = -log(1 - x * {1 - exp(-2)})
  return(res)
}



#######################################################################

# Generate 5000 samples

####################################################################


df_4 = 1:5000 %>% map(~inverse_transform_fn(f_cdf_inverse_fn))
df_4 = data.frame(id = 1:5000, x = unlist(df_4))


#####################################################################################

# Plot the density of sample

#####################################################################################

ggplot(df_4, aes(x = as.numeric(x))) + 
  geom_density(fill = "blue", alpha = 0.25, color = "red") +
  labs(x = "", y = "Sample")  +
  ggtitle("Q4. Density plot of f(x)") +
  theme(plot.title = element_text(hjust = 0.5))




#####################################################################################
#####################################################################################
#####################################################################################

# Question 5

#####################################################################################
#####################################################################################
#####################################################################################




#####################################################################################

# Acceptance - Rejection Algorithm

#####################################################################################



accp_rej_fn_5 = function(g, g_cdf_inv, alpha = 1, N = 200){
  
  
  df = data.frame(id = 1:N, draw = NA, decision = NA )
  
  f = function(x){
    res =  exp(-x) / { 1 + x^2}
    return(res)
  }
  
  
  f_g_ratio_fn = function(x){
    res = f(x)/g(x)
    return(res)
  }
  
  
  # compute supremum of ratio
  
  count = 1
  loop_count = 1
  cat("\n ")
  
  while(count <= N){
    loop_count = loop_count + 1
    if(loop_count == 7000){
      print("Max reached")
      break
    }
    # obtain y from the inverse transform of g
    y = inverse_transform_fn(g_cdf_inv)
    
    # sample uniform u
    u = runif(1)
    if(loop_count %% 50 == 0){
      cat("\n ", loop_count, " out of ", N, "\n")
    }
    
    if(u < {f_g_ratio_fn(y)/alpha}  ){
      df[count,] = c(count, y, "accepted")
      count = count + 1
    }
    
    if(count == N){
      cat("\n\n\nHurrah!!! Completed")
    }
  }
  
  return(df)
}


#######################################################################

# G_{1} function

####################################################################


g_1_fn = function(x){
  res = exp(-x)
  return(res)
}


#######################################################################

# G_{2} function

####################################################################


g_2_fn = function(x){
  res = 2 / { pi * (1 + x^2) }
  return(res)
}


#######################################################################

# Inverse of G_{1} CDF function

####################################################################

g_1_cdf_inverse_fn = function(x){
  res = -log(1 - x)
  return(res)
}


#######################################################################

# Inverse of G_{2} CDF function

####################################################################


g_2_cdf_inverse_fn = function(x){
  res = tan({pi / 2} * x)
  return(res)
}



#######################################################################

# 5 (a)

# draw a sample of 5000 random observations
# plot the estimated density function for 0 < x < 5.

####################################################################


df1 = accp_rej_fn_5(g_1_fn, g_1_cdf_inverse_fn, alpha = 1, N = 5000)
df2 = accp_rej_fn_5(g_2_fn, g_2_cdf_inverse_fn, alpha = pi/2, N = 5000)



df1 %>% filter(0 < draw & draw < 5) %>% 
  ggplot( aes(x = as.numeric(draw))) + 
  geom_density(fill = "blue", alpha = 0.25, color = "green") +
  labs(x = "", y = "sample") +
  ggtitle("Q5. Density plot of f using g_1") +
  theme(plot.title = element_text(hjust = 0.5))

df2 %>% filter(0 < draw & draw < 5) %>% 
  ggplot( aes(x = as.numeric(draw))) + 
  geom_density(fill = "red", alpha = 0.25, color = "green") +
  labs(x = "", y = "sample") +
  ggtitle("Q5. Density plot of f using g_2") +
  theme(plot.title = element_text(hjust = 0.5))


#######################################################################

# 5 (b).

# Benchmarking acceptance - rejection algorithm via G_{1} and G_{2} functions

####################################################################

benchmark_df = benchmark("g_1" = {
  accp_rej_fn_5(g_1_fn, g_1_cdf_inverse_fn, alpha = 1, N = 200)
},
"g_2" = {
  accp_rej_fn_5(g_2_fn, g_2_cdf_inverse_fn, alpha = pi/2, N = 200)
},
replications = 1000,
columns = c("test", "replications", "elapsed",
            "relative"))
cat("\n\n")

benchmark_df
 # test      replications elapsed  relative
 # 1  g_1         1000    90.96     1.000
 # 2  g_2         1000    117.69    1.294



#####################################################################################
#####################################################################################
#####################################################################################

# Question 6

#####################################################################################
#####################################################################################
#####################################################################################



#######################################################################

# sampler for from function g

####################################################################




g_mix_sample_fn <- function(theta){
  
  # numerator of weights
  numer_weight = c(2 * gamma(theta), gamma(theta + 0.5))
  # weights
  weight = numer_weight/sum(numer_weight)
  # sample
  i = sample(c(1,2), 1, prob = weight)
  
  if(i == 1){
    res = rgamma(n = 1, shape = theta, rate = 1)
  }else{
    res = rgamma(n = 1, shape = theta + 0.5, rate = 1)
  }
  
  return(res)
  
}




#####################################################################################

# Acceptance - Rejection Algorithm

#####################################################################################



accp_rej_fn_6 = function(theta, N = 200, type){
  
  theta = theta
  
  df = data.frame(id = 1:N, draw = 1:N, type = NA)
  
  f = function(x){
    res = sqrt(4 + x) * x^{theta -1} * exp(-x)
    return(res)
  }
  
  g = function(x){
    res = { 2 * x^{theta - 1} + x^{theta - 0.5} }  * exp(-x)
    return(res)
  }
  
  f_g_ratio_fn = function(x){
    res = f(x)/g(x)
    return(res)
  }
  
  
  # compute supremum of ratio
  alpha = 1
  count = 1
  loop_count = 1
  cat("\n ")
  
  while(count <= N){
    loop_count = loop_count + 1
    if(loop_count == 7000){
      print("Max reached")
      break
    }
    # obtain y from the mixture g
    y = g_mix_sample_fn(theta)
    
    # sample uniform u
    u = runif(1)
    
    
    if(loop_count %% 50 == 0){
      cat("\n ", loop_count, " out of ", N, "\n")
    }
    
    if(u < {f_g_ratio_fn(y) /alpha} ){
      df[count,] = c(count, y, type)
      count = count + 1
    }
    
    if(count == N){
      cat("\n\n\nHurrah!!! Completed")
    }
  }
  
  return(df)
}




#####################################################################################

# Function call for the acceptance - rejection algo

#####################################################################################



theta_1 = accp_rej_fn_6(0.5, 5000, 0.5)
theta_2 = accp_rej_fn_6(1, 5000, 1)
theta_3 = accp_rej_fn_6(1.5, 5000, 1.5)

theta_df = rbind(theta_1, theta_2, theta_3)
theta_df$type = factor(theta_df$type)


#####################################################################################

# Plot the sample using theta = 0.5

#####################################################################################

theta_df %>% filter(type == 0.5) %>% 
  ggplot( aes(x = as.numeric(draw))) + 
  geom_density(fill = "green", alpha = 0.25, color = "yellow") +
  labs(x = "", y = "sample") +
  ggtitle("Q6. Density plot of f using theta = 0.5") +
  theme(plot.title = element_text(hjust = 0.5))

#####################################################################################

# Plot the sample using theta = 1.0

#####################################################################################

theta_df %>% filter(type == 1) %>% 
  ggplot( aes(x = as.numeric(draw))) + 
  geom_density(fill = "blue", alpha = 0.25, color = "yellow") +
  labs(x = "", y = "sample") +
  ggtitle("Q6. Density plot of f using theta = 1") +
  theme(plot.title = element_text(hjust = 0.5)) 


#####################################################################################

# Plot the sample using theta = 1.5

#####################################################################################

theta_df %>% filter(type == 1.5) %>%  
  ggplot( aes(x = as.numeric(draw))) + 
  geom_density(fill = "red", alpha = 0.25, color = "yellow") +
  labs(x = "", y = "sample") +
  ggtitle("Q6. Density plot of f using theta = 1.5") +
  theme(plot.title = element_text(hjust = 0.5))


#####################################################################################

# Plot the three samples using theta = {0.5, 1.0, 1.5}

#####################################################################################


theta_df %>%
ggplot( aes(x = as.numeric(draw), colour = type, fill = type)) + 
  geom_density( alpha = 0.25, color = "green") +
  labs(x = "", y = "sample") +
  ggtitle("Q6. Density plot of f using theta = {0.5, 1.0, 1.5}") +
  theme(plot.title = element_text(hjust = 0.5))

