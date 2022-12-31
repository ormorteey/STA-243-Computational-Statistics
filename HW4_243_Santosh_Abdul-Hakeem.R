# Installing Pacman package for loading other packages

# install.packages("pacman")
library(pacman)

# Loading the needed packages

p_load(tidyverse)


# clear all variables from the environment before code starts

rm(list = ls())




#####################################################################################
#####################################################################################
#####################################################################################

# Question 1

#####################################################################################
#####################################################################################
#####################################################################################




#######################################################################

# 1(a)

####################################################################


set.seed(1)
U1 <- runif(10000)
V1 <- sapply(U1, function(x){x^2}) 
Integral1 <- mean(V1)
#0.3348666



#######################################################################

# 1(b)

####################################################################

U2 <- runif(10000, -1, 1)
V2 <- runif(10000, -1, 1)
points <- cbind(U2,V2)
Integrand2 <- function(x){
  if (x[1]^2+x[2]^2 <= 1)
    {value <- 1.0}
  else {value <- 0.0}
  return(value)
}
W <- apply(points,1,Integrand2)
Integral2 <- 4*mean(W)


#######################################################################

# 1(c)

####################################################################


#######################################################################
# Sign function
####################################################################

sgn <-function(x){ if (x>0)
{sgn <- 1}
  else {
    {sgn <- -1}
  }
  return(sgn)
} 

#######################################################################
# Cube Root function
####################################################################

cube_root <- function(x){sgn(x)*(abs(x))^(1/3)} 


# X be a random variable with distribution 3/4x^2e^{-x^3/4}, x>0
# Sampling $X$:


inv_cdf <- function(x){cube_root(-4*log(1-x))}
U3 <- runif(10000)
V3 <- (sapply(U3, inv_cdf))^2
Integral3 <- mean(V3)

# Method 2 for (c)
f1 <- function(x){(3/4)*x^4*(exp(x-(.25*x^3)))}
W1 <- rexp(100000)
S2 <- sapply(W1, f1)
Int1 <- mean(S2)



#####################################################################################
#####################################################################################
#####################################################################################

# Question 2

#####################################################################################
#####################################################################################
#####################################################################################



#######################################################################

# Weight function

####################################################################



f <- function(x){(1/sqrt(2*pi))*(exp(-x^2/2)*(x>1)*(x<2))
}


U3 <- rnorm(100,1.5,0.1)
U4 <- rnorm(10000,1.5,1)
U5 <- rnorm(10000,1.5,10)
 


#######################################################################
# Integral using g1 i.e. nu =0.1
####################################################################

 
g1 <- function(x){(1/(sqrt(2*pi)*.1))*exp((-0.5/(.1)^2)*(x-1.5)^2)}

weights1 <- sapply(U3, f)/sapply(U3, g1)

hist(weights1, main = 'histogram of g1-weights', xlab = 'weights')

Integral4 <- mean(weights1)

hist(weights1)

#######################################################################
# Integral using  g2 i.e. nu = 1
####################################################################


g2 <- function(x){(1/sqrt(2*pi))*exp(-0.5*(x-1.5)^2)}

weights2 <- sapply(U4, f)/sapply(U4, g2)

hist(weights2)

Integral5 <- mean(weights2) # best one 

#######################################################################
# Integral using  g3 i.e. nu = 10
####################################################################


g3 <- function(x){(1/(sqrt(2*pi)*10))*exp(-(0.5/10)*(x-1.5)^2)}

weights3 <- sapply(U5, f)/sapply(U5, g3)

hist(weights3)

Integral6 <- mean(weights3)

integrate(f,1,2)





#####################################################################################
#####################################################################################
#####################################################################################

# Question 3

#####################################################################################
#####################################################################################
#####################################################################################





#######################################################################

# 3(a)

####################################################################


U7 <- runif(1500)
h <- function(x){1/(1+x)}
weights4 <- sapply(U7, h)
Integral7 <- mean(weights4)





#######################################################################

# 3(b)

####################################################################

n <- 1500 
U7 <- runif(n)
h <- function(x){1/(1+x)}
c1 <- function(x){1+x}
weights5 <- sapply(U7, h)
weights6 <- sapply(U7, c1)
muMC <- mean(weights5)
thetaMC <- mean(weights6)
theta <- integrate(c1,0,1)
mu_CV <- function(b){muMC - b*(thetaMC-theta$value)}
#  optimal value of b 
b <- cov(weights5,weights6)/var(weights6)
# Estimate of the integral 
mu_CV(b)

df <-  data.frame(cbind(weights5,weights6)) # Checking correlation 
pairs(df)




#######################################################################

# 3(c)

####################################################################

Var_MC <- var(weights5)/n
Var_CV <- (var(weights5)-(cov(weights6,weights5)^2/(var(weights6))))/n

#  (d)  

c2 <- function(x){h(x)^2} # sin(x)
weights7 <- sapply(U7, c2)
New_Var_CV <- (var(weights5)-(cov(weights7,weights5)^2/(var(weights7))))/n








#####################################################################################
#####################################################################################
#####################################################################################

# Question 5

#####################################################################################
#####################################################################################
#####################################################################################




#######################################################################

# 5(a)

####################################################################

p = 0.3
lambda = 2
num_gen = 100 
r = rbinom(num_gen,1,p)
x = rep(0,num_gen)
for (i in 1:num_gen){if(r[i] == 1)
{
  x[i] = rpois(1,lambda * r[i])
}
  else { x[i]==0
    
  }
}




#######################################################################

# 5(c)

####################################################################




#######################################################################
# Gibbs Sampler function
####################################################################

gibbs_sampler_fn = function(a = 2, b = 2){
  
  p = 0.3
  lambda = 2
  
  num_gen = 100
  r = rbinom(num_gen,1,p)
  x = c()
  for(i in 1:num_gen){
    
    if(r[i] == 1){
      x = c(x, rpois(1, lambda * r[i]))
    }else{
      x = c(x,0)
    }
  }
  
  len_range = 10^4
  
  #initial values
  p_init = 0.5
  
  lambda_init = 0.5
  
  p_vec = c(p_init)
  lambda_vec = c(lambda_init)
  
  for(i in 1: len_range){
    r_vec = rep(NA, num_gen)
    
    for( j in 1:num_gen){
      
      
      if(x[j] == 0){
        r_vec[j] = { p_vec[i] * exp(-lambda_vec[i]) } / { p_vec[i] * exp(-lambda_vec[i]) + {1 - p_vec[i]} }
        
      }else{
        r_vec[j] = 1
      }
      
      
    }
    
    
    lambda_vec = c(lambda_vec, rgamma(1, a + sum(x), b + sum(r_vec)))
    
    p_vec = c(p_vec, rbeta(1, 1 + sum(r_vec), num_gen - sum(r_vec) + 1))
    
    
  }
  
  
  res = c(quantile(p_vec, 0.025),
          quantile(p_vec, 0.975),
          quantile(lambda_vec, 0.025),
          quantile(lambda_vec, 0.975))
  
  return(res)
  
}



#######################################################################
# Creating grid
####################################################################

sq = seq(1,5, length.out = 100)
sq_grid = expand.grid(sq, sq)
dim(sq_grid)



#######################################################################
# Computing 95% Bayesian Credible Interval
####################################################################


df5 = data.frame(i = 1:dim(sq_grid)[1], a = NA, b = NA, 
                "lwr(p_vec)" = NA,
                "upr(p_vec)" = NA,
                "lwr(lambda_vec) "= NA,
                "upr(lambda_vec)" = NA)

for(i in 1:dim(sq_grid)[1]){
  res = gibbs_sampler_fn(sq_grid[i,1],sq_grid[i,2])
  
  df5[i,] = c(i,sq_grid[i,1],sq_grid[i,2], res)
  
}


colnames(df5) = c("i", "a", "b", "lwr_p", "upr_p", "lwr_lambda", "upr_lambda")


#######################################################################
# Percentages of credible interval containing true value
####################################################################


p_per = dim( df %>% filter( lwr_p < 0.3 & 0.3 <	upr_p))[1]

lambda_per = dim( df %>%   filter( lwr_lambda < 2 & 2 <	upr_lambda) )[1]

p_lambda_per = dim( df %>% filter( lwr_p < 0.3 & 0.3 <	upr_p)   %>%  filter( lwr_lambda < 2 & 2 <	upr_lambda) )[1]




#####################################################################################
#####################################################################################
#####################################################################################

# Question 6

#####################################################################################
#####################################################################################
#####################################################################################


#######################################################################
# Target density function
####################################################################

target_dens <- function(z){ 
  theta1 = 1.5
  theta2 = 2
  if (z <= 0){value = 0} 
  else
  {value = (z>0)*z^(-3/2)*exp(theta1*z-(theta2/z)+2*sqrt(theta1*theta2)+log(sqrt(2*theta2)))}
  return(value)
}




#######################################################################
# Gamma Density Sampler Function: Sampling y independent from x using a gamma distribution that has the same mean as target
####################################################################

gamma_calc = function(alpha, beta){
  
  N <- 1000
  x <- numeric(N)
  
  x[1] <- rgamma(1, shape = alpha, scale = beta)
  k <- 0
  u <- runif(N)
  for (i in 2:N){ xt <- x[i-1]
  y <- rgamma(1, alpha, beta) #
  num <- target_dens(y)*dgamma(xt, shape = alpha, scale = beta)
  denom <- target_dens(xt)*dgamma(y,shape = alpha, scale = beta)
  r  <- min(c(num/denom,1))
  if (u[i] <= r) { x[i] <- y}
  else {x[i] <- xt}
  k <- k+1
  }
  m_x = mean(x)
  m_1_x = mean(1/x)
  
  return(c(m_x, m_1_x)) # estimate for x and 1/x
  
}


#######################################################################
# Exact True Values of the expectations
####################################################################
m_x_exact  = 1.154701 #estimate for x 
m_1_x_exact = 1.116025   #estimate for 1/x

# creating grid

sq = seq(1,5, length.out = 100)
sq_grid = expand.grid(sq, sq)

# running through grid to find best alpha and beta

df = data.frame(id = 1:10000, alpha = NA, beta = NA, m_X = NA, m_1_x = NA)
tol = 0.001
for( i in 1:10000){
  
  res = gamma_calc(sq_grid[i,1],sq_grid[i,2])
  df[i,] = c(i,sq_grid[i,1],sq_grid[i,2], res)
  
  if(abs(res[1] - m_x_exact) < tol  & abs(res[2] - m_1_x_exact) < tol){
    break
  }
  
}



#######################################################################
# Filtering for closest  estimates
####################################################################

# 

df_filter = df %>% filter( 1.13 < m_X & m_X < 1.16  )  %>% filter( 1.112 < m_1_x & m_1_x < 1.12  )



