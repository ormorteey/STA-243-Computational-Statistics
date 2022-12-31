# install.packages("pacman")
library(pacman)



# Loading the needed packages
p_load(tidyverse)


# clear all variables from the environment before code starts
rm(list = ls())

#################################################################
# Q1 
###################################################################
####################################################
# Data 
##################################################

Y <- c(576, 635, 558, 578, 666, 580, 555, 
       661, 651, 605, 653, 575, 545, 572, 594)
Z <- c(3.39, 3.30, 2.81, 3.03, 3.44, 3.07, 3.00, 3.43, 3.36, 
       3.13, 3.12, 2.74, 2.76, 2.88, 3.96)
X <- as.data.frame(cbind(Y,Z))


# (a) Sample estimate
rho <- cor(Y,Z)
print(rho)

# Generating bootsrap samples and computing statistic 
N <- nrow(X) 
B <- 200
cor.values <- c()
for (i in 1:B) {                             
  T <- X[sample(1:N, replace = TRUE),] 
  cor.values[i] = cor(T$Y,T$Z) 
}

#(b) 
# Bootstrap SE 
t <- mean(cor.values)
se_cor <- sqrt(sum((cor.values-t)^2)/(B-1)) 
se_b <- sd(cor.values)
hist(cor.values, main = '', xlab = 'estimates of correlation coefficient')

# A more direct computation: constructing a function that computes bootstrap se 
boot_se <- function(data,B){
  R = numeric(B)
  N = nrow(data)
  #X <- data[,1]
  #Y <- data[,2]
  for (b in 1:B) {
    i <- sample(1:N, size = N, replace = TRUE)
    X <- data[,1][i]
    Y <- data[,2][i]
    R[b] <- cor(X,Y)
  }
  se <- sd(R)
  return(se)
}

boot_se(X,200)


# Jacknife SE 

JN_cor <- c()
for (i in 1:N){
  T <- X[-i,]
  JN_cor[i] <- cor(T$Y, T$Z)
}

m <- mean(JN_cor) 
se_JN <- sqrt({N-1}*mean((JN_cor-m)^2))
print(se_JN)

#(c) 

# (i) Bootstrap CI via normal theory
alpha =c(.025,.975)
print(rho + qnorm(alpha) * se_cor) #0.1768888 0.9149490 


# studentized t-confidence interval
B = 200 
r_star <- numeric(B) # the set theta_star
se_r_star <- numeric(B) # set of standard errors of each members of theta_star
for ( b in 1:B){
  data = X
  N = nrow(X)
  i <- sample(1:N, size = N, replace = TRUE)
  Y <- data[,1][i]
  Z <- data[,2][i]
  r_star[b] <- cor(Y,Z)
  M <- numeric(200) # computing se(theta*(b))
  for (j in 1:200){
    data <- as.data.frame(cbind(Y,Z)) 
    i <- sample(1:N, size = N, replace = TRUE)
    U <- data[,1][i]
    V <- data[,2][i]
    M[j] <- cor(U,V)
  }
  se_r_star[b] <- sd(M)
}

t1 <- r_star-rho  
z_star <- t1/se_r_star  # set of pivotal 
level = .95 
alpha <- 1 - level
L <- as.numeric(quantile(z_star, alpha))
U <- as.numeric(quantile(z_star, 1-alpha))
se <- boot_se(X,200)

upp <- rho - (se*L) 
upp #0.7844813
low <- rho - (se*U)
low #-0.7099513 




#############################################################
#Q2
##############################################################


# (c) +(e) Parametric bootstrap 
X <- runif(50,0,3) # drawing sample from P_theta 
x_max <- max(X) # computation of thetahat 
n = 50 
theta = 3


par_theta_star <- function(B){
  R = numeric(B)
  for (b in 1:B){
    Y <- runif(n,0,x_max) # sample from P_thetahat
    R[b] <- max(Y)
  }
  
  return(R)
}

par_boot_theta_star <- par_theta_star(5000)
hist(par_boot_theta_star, xlab = 'theta_star', main ='')


var(par_boot_theta_star)
(n*(theta^2))/((n+2)*((n+1)^2))


# They are very close to each other 

# (d) +(e) Non parametric bootstrap 

# Function that computes se of maximum of a univariate data 
non_par_theta_star <- function(X,B){
  R = numeric(B)
  N = length(X)
  for (b in 1:B){
    i <- sample(1:N, size = N, replace = TRUE)
    X <- X[i]
    R[b] <- max(X)
  }
  return(R)
}

non_par_boot_theta_star <- non_par_theta_star(X,5000)
hist(non_par_boot_theta_star, main = '', xlab = 'nonparametric theta_star')
var(non_par_boot_theta_star)




# (f) 
Inv_cdf <- function(x){3*(x^(1/50))}
Inv_cdf(.001)
U <- runif(50,0,1)
Y <- sapply(U, Inv_cdf)
hist(Y)


# graph of $g$ 

g <- function(x){(50*(x^(49))/(3^(50)))}
grid <- seq(0,3,length.out = 300)
gvalue<- sapply(grid, g)
plot(grid,gvalue, type = 'l', xlab = 'x', ylab = 'g(x)')


