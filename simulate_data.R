library(tidyverse)

# Choose distribution and parameters
skew <- FALSE
df <- 5L # degrees of freedom for t distribution; ignored if skew=TRUE
a <- 4L # first component of alpha in skew-normal distribution; ignored if skew=FALSE
noise_multiplier <- 1/6 # values used: 1/6, 1/3, 2/3, 4/3

# Set some constants
K_true <- 75L # clusters
P <- 10L # markers
M <- 10L # samples
N <- 3e5L # events per sample
set.seed(1258) # reproducibility

# Generate phenotypes. Optional, tabulate number of positive markers
comp <- generate_cluster_phenos(K_true, P)
# count_positive <- apply(comp, 1, function(v) length(which(v == 2)))
# table(count_positive)

# Generate mixture proportions
pro_baseline <- rbeta(K_true, shape1 = 2, shape2 = 8) # baseline proportions
pro_by_sample <- seq(M) %>%
  lapply(function(x) pro_baseline * runif(K_true, min = 0.7, max = 1.3)) %>%
  lapply(function(x) x / sum(x)) %>%
  do.call(what = rbind) # add some random variation across samples

# Generate means
means_baseline <- cbind(rep(1.1,P) , rep(3,P))
sd_marginal <- get_sd_marginal(P) # legacy function; should be replaced with runif(2*P,0.05,0.2)

means_by_sample <- array(dim = c(M, K_true,P))
for (m in seq(M)) {
  for (p in seq(P)) {
    means_by_sample[m,,p] <- vapply(comp[,p], function(x) rnorm(1, mean = means_baseline[p,x], 
                                                                sd = noise_multiplier * sd_marginal[p,x]),
                               numeric(1))
  }
}

# Generate variances: for simplicity, the same for all samples
var_baseline <- array(dim = c(P,P,K_true))
for (k in seq_len(K_true)) {
  c <- comp[k,]
  d <- vapply(c(1:P), function(x) sd_marginal[x,c[x]], numeric(1))
  var_baseline[,,k] <- rWishart(1, df = P, Sigma = diag(d) / P)
}

# Generate data
data <- matrix(nrow = N * M, ncol = P)
names <- vapply(c(1:P), function(x) paste("M", x, sep = ""), character(1))
colnames(data) <- names
colnames(comp) <- names
true_labels <- integer(N * M)

start <- 1
for (m in seq(M)) {
  for (k in seq(K_true)) {
    # How many events drawn from this cluster?
    this_N <- ifelse(k == K_true, m*N - start + 1, floor(pro_by_sample[m,k] * N))
    stop <- start + this_N - 1
    
    if (skew) {
      alpha <- c(rep(a,1),rep(0,p-1))
      data[c(start:stop),] <- sn::rmsn(this_N, alpha = alpha,
                                       xi = means_by_sample[m,k,], 
                                       Omega = var_baseline[,,k])
    }
    else if (is.infinite(df)) {
      data[c(start:stop),] <- mvtnorm::rmvnorm(this_N, mean = means_by_sample[m,k,], 
                                              sigma = var_baseline[,,k])
    }
    else {
      scale <- 1- 2/ df # Keep the same variance regardless of df
      data[c(start:stop),] <- mvtnorm::rmvt(this_N, df = df, delta = means_by_sample[m,k,], 
                                            sigma = scale* var_baseline[,,k])
    }
    true_labels[c(start:stop)] <- k
    
    start <- stop + 1
  }
}

# Inspect data
suppressWarnings( plot_kdes_with_means(data, means_by_sample) )


###########
s <- apply(means_by_sample, c(2,3), sd)
noise_frac <- mean(s) / 1.9



