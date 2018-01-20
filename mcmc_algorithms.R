# Metropolis Algorithm --------------------------------------------------------------

flips <- 41; heads <- 13
a <- 2; b <- 2

metropolis_algorithm <- function(samples, theta_seed, sd){
  
  theta_curr <- theta_seed
  posterior_thetas <- rep(NA, samples) # Create vector to store sampled parameters
  
  for (i in 1:samples){
    
    theta_prop <- rnorm(n = 1, mean = theta_curr, sd = sd) # Proposal distribution
    
    # If the proposed parameter is outside its range then set it equal to its current
    # value. Otherwise keep the proposed value
    theta_prop <- ifelse((theta_prop < 0 | theta_prop > 1), theta_curr, theta_prop)
    
    # Bayes' numerators
    posterior_prop <- dbeta(theta_prop, a, b) * dbinom(heads, flips, theta_prop)
    posterior_curr <- dbeta(theta_curr, a, b) * dbinom(heads, flips, theta_curr)
    
    # Calculate probability of accepting
    p_accept_theta_prop <- min(posterior_prop/posterior_curr, 1.0)
    
    rand_unif <- runif(n = 1)
    
    # Probabilistically accept proposed theta
    theta_select <- ifelse(p_accept_theta_prop > rand_unif, theta_prop, theta_curr)
    
    posterior_thetas[i] <- theta_select
    
    # Reset theta_curr for the next iteration of the loop
    theta_curr <- theta_select
  }
  return(posterior_thetas)
}

set.seed(555)
posterior_thetas <- metropolis_algorithm(samples = 10000, theta_seed = 0.9, sd = 0.05)

opar <- par()
par(mar=c(2.5,3.5,3,2.1), mgp = c(1.7, 1, 0))

d <- density(posterior_thetas)
plot( d
      ,main = expression(paste('Kernel Density Plot for ', theta))
      ,xlab = expression(theta)
      ,ylab = 'density'
      ,yaxt = 'n'
      ,cex.lab = 1.3
)
polygon(d, col='dodgerblue1', border='dark blue')

a <- 20; b <- 20 # Prior parameters

# Create a sequence of thetas ranging from 0.1 to 0.9 since this is where the density is
min <- 0.1; max <- 0.9
thetas <- seq(from = min, to = max, length.out = 1000)

# Pass the sequence of thetas to calculate the likelihood, prior, and  posterior
likelihood <- dbinom(heads, flips, thetas)
prior <- dbeta(thetas, a, b)
posterior <- dbeta(thetas, (a + heads), (b + flips - heads))

# Metropolis Algorithm in 2-D -------------------------------------------------------

flips_1 <- 25; heads_1 <- 17
flips_2 <- 9; heads_2 <- 1
a <- 10; b <- 10

thetas <- seq(from = 0, to = 1.0, length.out = 50)

prior_1 <- dbeta(thetas, a, b)
prior_2 <- dbeta(thetas, a, b)

likelihood_1 <- dbinom(heads_1, flips_1, thetas)
likelihood_2 <- dbinom(heads_2, flips_2, thetas)

posterior_1 <- dbeta(thetas, (a + heads_1), (b + flips_1 - heads_1))
posterior_2 <- dbeta(thetas, (a + heads_2), (b + flips_2 - heads_2))

# Matrix of densities. Multiply (*) because of independence.
joint_prior <- outer(prior_1, prior_2, FUN='*' )
joint_likelihood <- outer(likelihood_1, likelihood_2, FUN='*')
joint_posterior <- outer(posterior_1, posterior_2, FUN = '*')

library(MASS) # Needed for mvrnorm() function

flips <- c(flips_1, flips_2) # Create a vector of flips
heads <- c(heads_1, heads_2) # Create a vector of heads

metropolis_algorithm_2d <- function(samples, thetas_seed, cov_mat){# 2D metropolis 
  # algorithm function
  
  thetas_curr <- thetas_seed
  posterior_thetas <- matrix(NA, nrow = samples, ncol = 2) # Create a matrix to store 
  # results
  for (i in 1:samples){# Similar for loop as before in 1-D
    
    # 2-dimensional proposal distribution
    thetas_prop <- mvrnorm(n = 1, mu = thetas_curr, Sigma = cov_mat)
    
    # If any of the values in thetas_prop is not in [0, 1.0] then replace their value with
    # the same element in thetas_curr, otherwise keep their values
    thetas_prop <- mapply(function(x, y) ifelse((x > 1 | x < 0), y, x),
                          x = thetas_prop, 
                          y = thetas_curr)
    
    # Calculate Bayes' numerators
    posterior_prop <- dbeta(thetas_prop, a, b) * dbinom(heads, flips, thetas_prop)
    posterior_curr <- dbeta(thetas_curr, a, b) * dbinom(heads, flips, thetas_curr)
    
    # Accept/reject logic
    p_accept_theta_prop = pmin(posterior_prop/posterior_curr, 1.0)
    
    # Draw 2 random uniform values
    rand_unif <- runif(n = 2)
    
    # If probability of accept is greater than the random uniform then accept 
    # thetas_prop otherwise return thetas_curr
    thetas_select <- mapply(function(x, y, w, z) ifelse(x > y, w, z),
                            x = p_accept_theta_prop,
                            y = rand_unif,
                            w = thetas_prop,
                            z = thetas_curr)
    
    posterior_thetas[i, ] <- thetas_select
    
    # Reset theta_curr for the next iteration of the loop
    thetas_curr <- thetas_select
    
  }
  return(posterior_thetas)
}

# Call the function
set.seed(225)
posterior_thetas_2d <- metropolis_algorithm_2d(20000, thetas_seed = c(0.5, 0.5)
                                               ,cov_mat = matrix(c(0.05^2, 0, 0, 0.05^2)
                                                                 ,ncol=2))

set.seed(481) # For reproducibility
warmup<- seq(5000) # 1,2,...,5000
posterior_thetas_2d_post_warmup <- posterior_thetas_2d[-warmup, ] # Remove warmup period

n <- nrow(posterior_thetas_2d_post_warmup) # Number of samples to draw from the exact
# marginal posterior distributions. We set this
# equal to the number of samples post warmup

exact_theta_1 <- rbeta(n, (a + heads_1), (b + flips_1 - heads_1)) # Random sample from the
# first marginal
exact_theta_2 <- rbeta(n, (a + heads_2), (b + flips_2 - heads_2)) # Random sample from the
# second marginal

apprx_theta_1 <- posterior_thetas_2d_post_warmup[ ,1]
apprx_theta_2 <- posterior_thetas_2d_post_warmup[ ,2]

# Combine into a matrix so that we can apply the quantile() function to each column of the
# matrix with one function call
matrx_thetas <- cbind(exact_theta_1, apprx_theta_1, exact_theta_2, apprx_theta_2)

# Let's look at the 2.5th, 5th, 25th, median, 75th, 95th, and 97.5th quantile
quant_thetas <- apply(matrx_thetas, MARGIN = 2, FUN = quantile, c(0.025,0.05, 0.25, 0.5
                                                                  ,0.75, 0.95, 0.975))

# MCMC Diagnostics ------------------------------------------------------------------

# Trace plots
# define our metropolis sampling process
samples <- 10000
theta_seeds <- c(0.05, 0.50, 0.95)
sd <- 0.075 

flips <- 41; heads <- 13 # data

a <- 10; b <- 10 # prior beta distribution parameters

set.seed(124)
chains <- mapply(metropolis_algorithm, samples = samples, theta_seed = theta_seeds
                 ,sd = sd)

# Autocorrelation
warmup <- seq(100) # Conservatively set the warmup period. We discard these samples from
# our correlogram

# Remove warmup samples
chains_post_warmup <- chains[-warmup, ]

# Apply the acf() across each column (margin 2) of the chains_post_warmup matrix
acfs <- apply(chains_post_warmup, MARGIN = 2, FUN = acf, plot=FALSE)

# Effective sample size
# Index into the first chain
chain1_post_warmup <- chains_post_warmup[, 1]
n <- length(chain1_post_warmup)

acf_1 <- acf(chain1_post_warmup, plot=FALSE, lag.max = 1000)[[1]] # Grab the
# autocorrelations at
# each lag

# Loop through each consecutive sum of autocorrelations while the sums are non-negative. 
# The resulting t is our terminal point
acf_sums <- 999 # Initialize to be a positive number
t <- 1
while (acf_sums >= 0){
  t <- t + 1
  acf_sums <- acf_1[t - 1] + acf_1[t]
}

# Compute the denominator. The sum of the acf_1 vector is from the 2nd element to element 
# t because the second element is the first lag, k = 1
denominator <- 1 + 2*sum(acf_1[2:t]) 

# Compute an approximation of the effective number of samples, n_eff
n_eff <- n / denominator

# Thining
h = 5 # Thin our chain by keeping only 1 out of every 5 samples
n <- length(chain1_post_warmup)
keep_seq <- seq(1, n, h) # Sequence of integers that represent the indices of the chain
# we will keep
chain1_post_warmup_thin <- chain1_post_warmup[keep_seq]

thin_n <- length(chain1_post_warmup_thin)

# Extract autocorrelations at each lag
acf_1_thin <- acf(chain1_post_warmup_thin, plot=FALSE, lag.max = thin_n/5)[[1]]

acf_sums <- 999
t <- 1
while (acf_sums >= 0){
  t <- t + 1
  acf_sums <- acf_1_thin[t - 1] + acf_1_thin[t]
}

thin_denominator <- 1 + 2*sum(acf_1_thin[2:t])
thin_n_eff <- thin_n / thin_denominator

# Gelman-Rubin Convergence Diagnostic
# Helper function to split the chains in half
split_chain <- function(x){
  len <- length(x)
  split_x <- split(x, seq_along(x) <= len/2)
  return(split_x)
}

compute_Rhat <- function(chains){# Function to compute Gelman-Rubin Diagnostic
  
  chains_split <- unlist(apply(chains, MARGIN = 2, FUN = split_chain)
                         , recursive = FALSE) # Split chains in half
  
  n <- length(chains_split[[1]]) # Get number of samples in a split chain
  
  Ws <- sapply(chains_split, var) # Compute the variance of each chain. These are the 
  # components of W
  W <- mean(Ws) # Take the mean of the components of W. This is W
  
  chain_means <- sapply(chains_split, mean) # Compute the mean of each chain
  
  B <- var(chain_means) # Compute the variance of the chain means to get B
  
  V <- (n-1)/n * W + B
  
  Rhat <- sqrt(V / W)
  
  return(Rhat)
}

# Create a sequence of terminal number of samples. We will compute Gelman-Rubin
# after 10 samples, 60 samples, 110 samples,...,total_n samples
total_n <- nrow(chains_post_warmup)
terminal_n <- seq(from = 10, to = total_n, by = 50)

# Helper function that will be called in the subsequent loop to limit the chain size to
# `n` samples
truncate_chain <- function(chain, n){
  new_chain <- chain[1:n]
  return(new_chain)
}

Rhats <- rep(0, length(terminal_n)) # Create vector of 0's to store results 

# Loop through each terminal_n by truncating the chains at the value of n and then calling
# the compute_Rhat function. Store the results in the Rhats vector that we initialized 
# above
i <- 1
for (n in terminal_n){
  truncated_chains <- apply(chains_post_warmup, MARGIN = 2, FUN = truncate_chain, n)
  Rhats[i] <- compute_Rhat(truncated_chains)
  i <- i + 1
}

# Acceptance Rates
# We'll provide two ways to do this. The first might be more intuitive to some readers as 
# it better illustrates what we're actually doing. The second way is less code but perhaps
# a bit more opaque

# Intuitive way
compute_accRate_intuitive <- function(chain){
  n <- length(chain) # Number of samples
  chain_front <- chain[-n] # Create a chain with element 1,2,3,..,n-1
  chain_back <- chain[-1] # Create a chain with element 2,3,4,...n
  accepts <- sum(chain_front != chain_back) # Count how many successive elements are 
  # different from the next element
  acceptance_rate = accepts / (n - 1) # We subtract 1 because we could only make (n-1) 
  # comparisons
  return(acceptance_rate)
}

# Opaque way. Employs R's unique() function which keeps only unique elements in a vector
compute_accRate_opaque <- function(chain){
  n <- length(chain)
  acceptance_rate = (length(unique(chain)) - 1) / (n - 1) # Subtract 1 from numerator
  # because the first theta will
  # always get counted as a unique
  # value
  return (acceptance_rate)
}

# Apply both functions to our matrix of chains. Results should be equivalent
acceptance_rate_intuitive <- apply(chains_post_warmup, MARGIN = 2, 
                                   FUN = compute_accRate_intuitive)

acceptance_rate_opaque <- apply(chains_post_warmup, MARGIN = 2, 
                                FUN = compute_accRate_opaque)
