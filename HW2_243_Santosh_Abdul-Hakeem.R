# Installing Pacman package for loading other packages

# install.packages("pacman")
library(pacman)

# Loading the needed packages

p_load(tidyverse, magrittr)


# clear all variables from the environment before code starts

rm(list = ls())

#####################################################################################
#####################################################################################
#####################################################################################

# Question 1. Simulate Annealing Algorithm

#####################################################################################
#####################################################################################
#####################################################################################




#####################################################################################

# Create Data Matrix

#####################################################################################

row_1 = c(0 ,1 ,2 ,4 ,9 ,8 ,3 ,2 ,1 ,5 ,7 ,1 ,2 ,9 ,3)
row_2  = c(1 ,0 ,5 ,3 ,7 ,2 ,5 ,1 ,3 ,4 ,6 ,6 ,6 ,1 ,9)
row_3  = c(2 ,5 ,0 ,6 ,1 ,4 ,7 ,7 ,1 ,6 ,5 ,9 ,1 ,3 ,4)
row_4  = c(4 ,3 ,6 ,0 ,5 ,2 ,1 ,6 ,5 ,4 ,2 ,1 ,2 ,1 ,3)
row_5  = c(9 ,7 ,1 ,5 ,0 ,9 ,1 ,1 ,2 ,1 ,3 ,6 ,8 ,2 ,5)
row_6  = c(8 ,2 ,4 ,2 ,9 ,0 ,3 ,5 ,4 ,7 ,8 ,3 ,1 ,2 ,5)
row_7  = c(3 ,5 ,7 ,1 ,1 ,3 ,0 ,2 ,6 ,1 ,7 ,9 ,5 ,1 ,4)
row_8  = c(2 ,1 ,7 ,6 ,1 ,5 ,2 ,0 ,9 ,4 ,2 ,1 ,1 ,7 ,8)
row_9  = c(1 ,3 ,1 ,5 ,2 ,4 ,6 ,9 ,0 ,3 ,3 ,5 ,1 ,6 ,4)
row_10  = c(5 ,4 ,6 ,4 ,1 ,7 ,1 ,4 ,3 ,0 ,9 ,1 ,8 ,5 ,2)
row_11  = c(7 ,6 ,5 ,2 ,3 ,8 ,7 ,2 ,3 ,9 ,0 ,2 ,1 ,8 ,1)
row_12  = c(1 ,6 ,9 ,1 ,6 ,3 ,9 ,1 ,5 ,1 ,2 ,0 ,5 ,4 ,3)
row_13  = c(2 ,6 ,1 ,2 ,8 ,1 ,5 ,1 ,1 ,8 ,1 ,5 ,0 ,9 ,6)
row_14  = c(9 ,1 ,3 ,1 ,2 ,2 ,1 ,7 ,6 ,5 ,8 ,4 ,9 ,0 ,7)
row_15  = c(3 ,9 ,4 ,3 ,5 ,5 ,4 ,8 ,4 ,2 ,1 ,3 ,6 ,7 ,0)



matrix_df = matrix(rbind(row_1,row_2,row_3,row_4,row_5,row_6,row_7,row_8,row_9,row_10,row_11,row_12,row_13,row_14,row_15), nrow = 15, ncol = 15)
matrix_df =t(matrix_df)

matrix_df
distance_matrix <- matrix_df

#######################################################################

# Objective function

####################################################################

Total_dist <- function(x){n = length(x)
distance <- distance_matrix
ynext <- 0
 for (i in 2:n){a <- x[i-1]
 b <- x[i]
 ynext <- ynext + distance[a,b]
   }
   a <- x[1]
   b <- x[n]
   ynext <- ynext +distance[a,b]
   result <- ynext
   return(result)
}

#####################################################################################

# Simulated annealing function

#####################################################################################

Sim_Ann_TSP_Modified_v1 <- function(temp, temp_min= .01, rate=.999){
   n = 15
  xbest <- xcurr <- xnext <- c(1:n)
  ynext <- Total_dist(xnext)
  ybest <- ycurr <- ynext
  while (temp>temp_min){
      for (m in 1:100){ # looping 100 times before temperature update
        i<- ceiling(runif(1,1,n))
        xnext <- xcurr
        temporary <- xnext[i]
        xnext[i] <- xnext[i-1]
        xnext[i-1] <- temporary
        ynext <- Total_dist(xnext)
        Delta <- ynext - ycurr
        accept <- exp(-Delta/temp)
        if (Delta < 0 || runif(1) < accept) {
          xcurr <- xnext
          ycurr <- ynext
       }
       if (ynext <- ybest){
        xbest <- xcurr
        ybest <- ycurr
      }
    }
    temp<-temp*rate  # updating temperature
  }
  result =  list()
  result[1] = list(xbest)
  result[2] = ybest
  return(result)
}

#####################################################################################

# Result that finds an optimal path of length 17

#####################################################################################

set.seed(0)
for (i in 1:20){print(Sim_Ann_TSP_Modified_v1(900))
}
#  3  5  7 14  4  6 13 11 15 10 12  8  2  1  9, length=17

#####################################################################################

# Checking effect of high temperature

#####################################################################################

set.seed(2)
for (i in 1:20){print(Sim_Ann_TSP_Modified_v1(1500))
}

#####################################################################################

# Checking effect of low temperature

#####################################################################################

set.seed(3)
for (i in 1:20){print(Sim_Ann_TSP_Modified_v1(500))
}









#####################################################################################
#####################################################################################
#####################################################################################

# Question 2. Genetic Algorithm

#####################################################################################
#####################################################################################
#####################################################################################



#####################################################################################

# Generate Test data

#####################################################################################


truefunction = function(x){
  t = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
  h = c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
  temp = 0
  for(i in 1:11) {
    temp = temp + h[i]/2 * (1 + sign(x - t[i]))
  }
  return(temp)
}

generate_test_data = function(n = 512){
  x = (0:(n-1))/n
  f = truefunction(x)
  set.seed(243)
  y = f+rnorm(f)/3
  
  data_result_gen = data.frame(x = x, y = y, f = f)
  return(data_result_gen)
  
}


#####################################################################################

# Randomly generate an initial population of chromosomes of size S.

#####################################################################################

generate_pop_unif =  function(size, df){
  
  population <- lapply(1:size, function(population) round(runif(nrow(df),0,1),0))
  
  return(population)
}

#####################################################################################

# Fitness function function

#####################################################################################

fitness_function = function(s,x,y, n_pop){
  
  # S - is chromosome
  # npop is the number of population
  
  df = data.frame(x = x, y = y)
  ind_ones = which(s == 1)
  b_s = ind_ones / length(s)
  df_result = data.frame( id = 1: {length(b_s) + 1}, f_j = NA, f_j_var = NA,   n_j = NA, B = NA, n_pop = NA )
  
  for(i in 1: {length(b_s) + 1}){
    
    if(i == 1){
      a = 0; b = b_s[i]
    }
    else if(i == {length(b_s) + 1}){
      a = b_s[i - 1]; b = 1
    } 
    else{
      a = b_s[i-1] ; b = b_s[i]
    }
    
    res = df[which(a <= x & x < b), 2]
    
    f_j = mean(res)
    n_j = length(res)
    f_j_var = sum({res - f_j}^2)
    B = length(ind_ones) + 1
    df_result[i, 2:6] = c(f_j, f_j_var, n_j, B, n_pop)
    # df_result[i, 2:5] = c( f_j, f_j_var, n_j, B)    
  }
  return(df_result)
}

#####################################################################################

# Compute fitness score using MDL

#####################################################################################

MDL_compute =  function(df){
  n = sum( df[,4])
  B = df[1,5] 
  result = B * log(n) + { 0.5 * sum ( log( df[, 4]) )} + { 0.5 * n * log( {1/n} * sum(df[, 3])) }
  return(result)
}

#####################################################################################

# compute fitness score using AIC

#####################################################################################

AIC_compute = function(df){
  n = sum( df[,4])
  B = df[1,5] 
  result =  { n * log( {1/n} * sum(df[, 3])) } + {2 * B * log(n)}
  return(result)
}

#####################################################################################

# Perform crossover or mutate

#####################################################################################

crossover_mutate_fn = function(pop_chr, df, p_cross, p_c){
  prob_check = runif(1, min = 0, max = 1)
  pop_chr_len = length(pop_chr)
  df$proba = df$rank / sum(df$rank)
  
  if(prob_check < p_cross){
    chr_ind = sample(1:pop_chr_len, 2 , prob = df$proba)
    result = crossover_assign_fn(pop_chr[[chr_ind[1]]], pop_chr[[chr_ind[2]]])
    return(result)
  }
  else{
    chr_ind = sample(1:pop_chr_len, 1 , prob = df$proba)
    chr = pop_chr[[chr_ind]]
    result = mutate_fn(chr,p_c)
    return(result)
  }
  
}

#####################################################################################

# Perform uniform crossover of two parents

#####################################################################################

crossover_assign_fn = function(chr1, chr2){
  
  chr_ch = c()
  for( i in 1:length(chr1)){
    chr_ch[i] = ifelse(chr1[i] == chr2[i], chr1[i], rbinom(1, 1, 0.5))
  }
  df_chr = data.frame(chr1 = chr1, chr2 = chr2, new_chr = chr_ch)
  return(df_chr)
}

#####################################################################################

# Perform mutation of a child

#####################################################################################

mutate_fn = function(chr, p_c = 0.05){
  df_chr = data.frame(id = 1:length(chr), chr = chr, new_chr = NA)
  n = length(chr)
  prob_vec = runif(n, min = 0, max = 1)
  for (ii in 1:n){
    if(prob_vec[ii] < p_c){
      chr[ii] = ifelse(chr[ii] == 0, 1, 0)
    }
  }
  df_chr$new_chr = chr
  return(df_chr)
}

#####################################################################################

# Copy best three and another 10% of the population into the next generation

#####################################################################################

copied_over_fn = function(pop_chr, df){
  
  pop_chr_len = length(pop_chr)
  df$proba = df$rank / sum(df$rank)
  best_set = df$ix[{pop_chr_len - 2}:pop_chr_len]
  chr_ind = sample(1:pop_chr_len, round(0.1 * pop_chr_len) , prob = df$proba)
  
  chr_ind = c(chr_ind, best_set )
  result = pop_chr[chr_ind]
  return(result)
  
}

#####################################################################################

# Compute the data for piecewise plot

#####################################################################################

piecewise_compute_function = function(s,x,y, n_pop, minimizer_fitness_df){
  
  # S - is chromosome
  # npop is number of population
  df = data.frame(x = x, y = y)
  df_piece = data.frame(x = NA, y = NA)
  ind_ones = which(s == 1)
  b_s = ind_ones / length(s)
  
  df_result = data.frame(id = 1: length(x))
  
  for(i in 1: {length(b_s) + 1}){
    
    if(i == 1){
      a = 0; b = b_s[i]
    }
    else if(i == {length(b_s) + 1}){
      a = b_s[i - 1]; b = 1
    } 
    else{
      a = b_s[i-1] ; b = b_s[i]
    }
    
    res_minimizer_df = df[which(a <= x & x < b), ]
    f_j_minimizer = minimizer_fitness_df$f_j[i]
    res_minimizer_df$y = rep(f_j_minimizer, nrow(res_minimizer_df))
    df_piece = rbind(df_piece, res_minimizer_df)
    
  }
  df_piece  = df_piece[2:nrow(df_piece),]
  df_result = cbind(df_result, df_piece)
  return(df_result)
}


#####################################################################################

# Genetic Algorithm

#####################################################################################

genetic_algorithm_fn = function(noisy_data, mtd = 1){
  
  # set population size for chromosome
  population_size = 300
  
  if( {"x" %in% names(noisy_data)} & {"y" %in% names(noisy_data)}){
    x = noisy_data$x
    y = noisy_data$y
  }else{
    x = noisy_data[,1]
    y = noisy_data[,2]
  }
  
  # S_pop = generate_pop(population_size, noisy_data)
  S_pop = generate_pop_unif(population_size, noisy_data)
  stop_criterion_external = 1
  
  # Set while condition 
  while_cond = TRUE
  
  # N_same initialization
  N_same = c()
  min_mbl_new = -1
  
  mult_count = 4
  cat("\n Number of iterations out of ")
  cat(mult_count * population_size)
  cat("\n")
  
  
  while(while_cond){
    
    #setting new min as old
    min_mbl_old = min_mbl_new
    
    mbl = c()
    f_val = data.frame(id = NA, f_j = NA, f_j_var = NA,   n_j = NA, B = NA, n_pop = NA )
    
    # compute MBL or AIC for each chromosome on the population
    
    compute_mtd = ifelse(mtd == 1, MDL_compute, AIC_compute)
    
    for(i in 1:length(S_pop)){
      fit_df = fitness_function(S_pop[[i]],x,y,i)
      fit_df = fit_df %>% filter(f_j != "NaN")
      f_val = rbind(f_val, fit_df)
      mbl[i] = compute_mtd(fit_df)
    }
    
    mbl_rank_df = sort(mbl, decreasing = TRUE, index.return = TRUE) %>% as.data.frame()
    temp_df = data.frame(rank = 1:nrow(mbl_rank_df))
    mbl_rank_df %<>% cbind(temp_df,.)
    best_mbl = mbl_rank_df[population_size,]
    best_mbl
    
    min_mbl_new = best_mbl$x
    # round to a value for degree of convergence
    min_mbl_new = round(min_mbl_new, 1)
    
    
    
    
    # Inital conditions
    
    p_cross = 0.9
    p_c = 0.05
    S_pop_new = copied_over_fn(S_pop, mbl_rank_df)
    stop_criterion_internal = length(S_pop_new) + 1
    
    
    while(length(S_pop_new) < length(S_pop)){
      new_chr_df = crossover_mutate_fn(S_pop, mbl_rank_df, p_cross, p_c)
      S_pop_new[length(S_pop_new) + 1] = list(new_chr_df$new_chr)
      S_pop_new = unique(S_pop_new)
      stop_criterion_internal = stop_criterion_internal + 1
      if(stop_criterion_internal >  {length(S_pop) + 20}){
        cat("\n\n Problem with generating population \n\n")
        break
      }
    }
    
    
    
    if( min_mbl_old == min_mbl_new ){
      N_same = c(N_same,min_mbl_new)
    }
    
    
    # Num of successive iterations to be equal for valid convergence as conv_n_same
    conv_n_same = 20
    if(length(N_same) >= {conv_n_same + 3}){
      N_same_unique = N_same[{length(N_same) - conv_n_same}:length(N_same)]
      N_same_unique = unique(N_same_unique)
      if(length(N_same_unique) == 1){
        cat("\n")
        cat("\n")
        cat("Hurrah!! Convergence Attained")
        cat("\n")
        cat("\n")
        break
      }
    }
    
    if(stop_criterion_external ==  3 * population_size){
      cat("\nExternal break for 3 * Pop\n")
      break
    }
    S_pop = S_pop_new
    
    cat(stop_criterion_external)
    cat(" out of ")
    cat(mult_count * population_size)
    cat("\n")
    stop_criterion_external = stop_criterion_external + 1
    
  }
  
  mbl_rank_df[population_size,]
  minimizer_id = mbl_rank_df$ix[population_size]
  minimizer_chrm = S_pop_new[minimizer_id]
  minimizer_fitness_df = fitness_function(S_pop_new[[minimizer_id]],x,y,population_size)
  minimizer_fitness_df = minimizer_fitness_df %>% filter(f_j != "NaN")
  plot_df = piecewise_compute_function(S_pop_new[[minimizer_id]],x,y,population_size, minimizer_fitness_df)
  
  # Plot the noisy data and the piecewise approximation
  
  if(mtd == 1){
    plot(x,y, main = "Plot of Noisy data and MBL fit")
  }else{
    plot(x,y, main = "Plot of Noisy data and AIC fit")
  }
  
  lines(plot_df$x, plot_df$y)
  
}

#####################################################################################

# End of Genetic Algorithm

#####################################################################################


#####################################################################################

# Test Genetic Algorithm for AIC and MDL

#####################################################################################

trydata = generate_test_data(300)

# Three plots should be generated
# genetic_algorithm_fn(noisy data, method indicator)
# Method indicator is either 1 for MDL or 2 for AIC

genetic_algorithm_fn(trydata, 1) # MDL
genetic_algorithm_fn(trydata, 2) # AIC
trydata %$% plot(x,y, main = " Plot of test data and denoised data ") # Noisy data
trydata %$% lines(x,f) # Original Denoised Data
