generate_cluster_phenos <- function(K_true, P) {
  alpha <- rep(0.5, P)
  mat <- DirichletReg::rdirichlet(2 * K_true, alpha)
  comp <- matrix(vapply(mat, function(v) ifelse(v > 1.6/P, 2L, 1L), integer(1)), nrow = 2 * K_true)
  comp <- choose_unique_phenos(comp, K_true - 1)
  comp <- rbind(rep(1L,P) , comp)
  return(comp)
}

choose_unique_phenos <- function(comp, n_unique) {
  phenos <- apply(comp, 1, function(x) paste(x, collapse = ""))
  unique_pheno_list <- split(seq_along(phenos), phenos)
  unique_indices <- as.integer(vapply(unique_pheno_list, function(x) x[1], integer(1)))
  # sel <- sample(unique_indices,n_unique)
  sel <- unique_indices[seq(n_unique)]
  
  return(comp[sel,])
}

count_pheno <- function(comp) {
  phenos <- apply(comp, 1, function(x) paste(x, collapse = ""))
  return(length(unique(phenos)))
}

get_sd_marginal <- function(P) {
  col <- runif(P, min = 0.6, max = 0.95)
  pros_marginal <- cbind(col, 1-col)
  sd_marginal <- 1/sqrt(pros_marginal) * mean(sqrt(pros_marginal)) * 0.5
  
  return(sd_marginal)
}


plot_kdes_with_means <- function(data, means) {
  plot_list <- list()
  P <- ncol(data)
  
  for (p in seq(P)) {
    # kde <- KernSmooth::bkde(data[,p])
    sel <- which(data[,p]>-1 & data[,p]<5)
    kde <- KernSmooth::bkde(data[sel,p])
    
    df_kde <- data.frame(x = kde$x, y = kde$y)
    df_means <- data.frame(m = as.vector(means[,,p]))
    
    g <- ggplot() +
      geom_line(data = df_kde, mapping = aes(x=x, y=y)) +
      geom_point(data = df_means, mapping = aes(x=m, y=0), size = 2, alpha = 0.3) +
      xlim(-1,5) +
      labs(title = names[p], x = "", y = "") +
      theme_bw()
    
    plot_list[[p]] <- g
  }
  
  return(gridExtra::grid.arrange(grobs = plot_list, ncol = 5))
}
