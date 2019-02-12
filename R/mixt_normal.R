#' random sampling of mixtured normal distribution
#' 
#' @param n sample size
#' @param mu population means
#' @param sigma population sd
#' @param pi mixtured rate
#' 
#' @importFrom dplyr %>%
#' @importFrom purrr pmap_dbl
#' @importFrom stats rnorm
#' @export
#' 
random_mixt_normal <- function(n, mu, sigma, pi){
  # elementary statistics
  n_clusters <- length(mu)
  pi <- pi / sum(pi)
  # sanity check (equality of dimensions)
  if (length(sigma) != n_clusters | length(pi) != n_clusters){
    stop("Error : dimension of mu, sigma and pi must be the same.")
  }
  # cluster that each point are sampled
  cluster_idxs <- sample(x = 1:n_clusters, size = n, prob = pi, replace = TRUE)
  # sampling
  obj <- list(cluster_idx = cluster_idxs) %>%
    pmap_dbl(.f = function(cluster_idx){
      rnorm(1, mean = mu[cluster_idx], sd = sigma[cluster_idx])
      })
  # output
  return(obj)
}


#' probability density function of mixtured normal distribution
#' 
#' @param x sample point
#' @param mu population means
#' @param sigma population sd
#' @param pi mixtured rate
#'
#' @importFrom dplyr %>%
#' @importFrom purrr pmap_dbl
#' @importFrom purrr invoke_map_dbl
#' @importFrom stats dnorm
#' @export
#' 
density_mixt_normal <- function(x, mu, sigma, pi){
  n_clusters <- length(mu)
  pi <- pi / sum(pi)
  # sanity check
  if (length(sigma) != n_clusters | length(pi) != n_clusters){
    stop("Error : dimension of mu, sigma and pi must be the same.")
  }
  # calculation of probability density
  # for each sample points
  density_each <- function(x){
    list(cluster_id = 1:n_clusters) %>%
    pmap_dbl(.f = function(cluster_id){
      pi[cluster_id] * dnorm(x, mean = mu[cluster_id], sd = sigma[cluster_id])
      }) %>%
    sum()
  }
  # For all sample points
  obj <- invoke_map_dbl(.f = density_each, .x = x)
  # output
  return(obj)
}


#' LL function : log likelihood of mixtured normal distribution
#' 
#' @param x sample point
#' @param mu population means
#' @param sigma population sd
#' @param pi mixtured rate
#' 
#' @importFrom dplyr %>%
#' @importFrom purrr pmap_dbl
#' @importFrom purrr invoke_map_dbl
#' @importFrom stats dnorm
#' @export
#' 
LL_mixt_normal <- function(x, mu, sigma, pi){
  # elementary statistics
  n_clusters <- length(mu)
  pi <- pi / sum(pi)
  # sanity check
  if (length(sigma) != n_clusters | length(pi) != n_clusters){
    stop("Error : dimension of mu, sigma and pi must be the same.")
  }
  # calculation of log likelihood
  # For each sample points
  LL_each <- function(x){
    list(cluster_id = 1:n_clusters) %>%
      pmap_dbl(.f = function(cluster_id){
        pi[cluster_id] * dnorm(x, mean = mu[cluster_id], sd = sigma[cluster_id])
        }) %>%
      sum() %>%
      log()
  }
  # For all sample points
  obj <- invoke_map_dbl(.f = LL_each, .x = x) %>% sum()
  # output
  return(obj)
}


#' EM function : EM algorithm for mixtured normal distribution
#' 
#' @param x sample
#' @param max_iter maximum number of iterations
#' @param tol tolerance
#' @param init_mu initial values of population means
#' @param init_sigma initial values of population sds
#' @param init_pi initial values of population rate
#' 
#' @importFrom dplyr data_frame
#' @importFrom purrr map
#' @importFrom purrr map_chr
#' @importFrom stats dnorm
#' @export
#' 
EM_mixt_normal <- function(x, max_iter, tol, init_mu, init_sigma, init_pi){
  # elementary statistics
  sample_size <- length(x)    # sample size
  n_clusters <- length(init_mu)
  init_pi <- init_pi / sum(init_pi)
  # sanity check
  if (length(init_sigma) != n_clusters | length(init_pi) != n_clusters){
    stop("Error : dimension of mu, sigma and pi must be the same.")
  }
  # settings of initial values
  mu <- init_mu; sigma <- init_sigma; pi <- init_pi
  # history
  LL_history <- LL_mixt_normal(x = x, mu = mu, sigma = sigma, pi = pi)
  # EM algorithm
  for (iter in 1:max_iter) {
    gamma <- x %>%
      map(.f = function(x){
        pi * dnorm(x = x, mu, sigma) / density_mixt_normal(x = x, mu, sigma, pi)
        }) %>%
      unlist() %>%
      matrix(byrow = TRUE, ncol = n_clusters)
    colnames(gamma) <- as.character(1:n_clusters) %>%
      map_chr(.f = ~paste0("gamma", .))
    sum_gamma <- colSums(gamma)
    # MLE of pi
    pi <- sum_gamma / sample_size
    # MLE of mu 
    mu <- as.vector(x %*% gamma) / sum_gamma
    # MLE of sigma
    ncol <- length(x); nrow <- n_clusters
    x_matrix <- matrix(x, byrow = TRUE, nrow = nrow, ncol = ncol)
    mu_matrix <- matrix(mu, nrow = nrow, ncol = ncol)
    sigma <- sqrt(diag((x_matrix-mu_matrix)%*%(gamma*t(x_matrix-mu_matrix))) / sum_gamma)
    # recoding history...
    LL_history <- append(LL_history, LL_mixt_normal(x = x, mu = mu, sigma = sigma, pi = pi))
    # coverage
    if(abs(LL_history[iter+1]-LL_history[iter]) < tol){
      break    
    }
  }
  LL_history_df <- data_frame(iter = 0:iter, LL_history = LL_history)
  estimated_cluster <- max.col(gamma)
  obj <- list(mu = mu, sigma = sigma, pi = pi, LL = LL_history_df, estimated_cluster = estimated_cluster, n_iter = iter)
  return(obj)
}