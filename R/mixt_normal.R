## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))


#' Random Sampling of 1-dim. Mixtures of Normal Distribution
#' 
#' Set sample size n and parameters of MixtNormal(mu,sigma,pi).
#' By definition, dimension of mu, sigma and pi must be the same.
#' In the function, parameter `pi` is normalized.
#' @param n <int scalar> : sample size
#' @param mu <double vector> : population means
#' @param sigma <double vector> : population sds
#' @param pi <double vector> : mixtured rate
#' @return The output will be <double vector> of size n that you input.
#' 
#' @importFrom dplyr %>%
#' @importFrom purrr pmap_dbl
#' @importFrom stats rnorm
#' @export 
#' 
random_mixt_normal <- function(n, mu, sigma, pi){
  # 1. Preliminaries and Error Corrections
  # record a number of components from dimension of the parameter `mu`
  n_components <- length(mu)
  pi <- pi / sum(pi)
  # sanity check for dimension of parameters
  if (length(sigma) != n_components | length(pi) != n_components){
    stop("Error : dimension of mu, sigma and pi must be the same.")
  }
  # 2. Main Calculation
  # decide a component that each point is sampled...
  component_idx <- sample(x = 1:n_components, size = n, prob = pi, replace = TRUE)
  # sampling
  obj <- list(component_idx = component_idx) %>%
    pmap_dbl(.f = function(component_idx){
      rnorm(1, mean = mu[component_idx], sd = sigma[component_idx])
      })
  # 3. output
  return(obj)
}





#' Probability Density Function of 1-dim. Mixtures of Normal Distribution
#' 
#' Set sample point x and parameters of MixtNormal(mu,sigma,pi).
#' By definition, dimension of mu, sigma and pi must be the same.
#' In the function, parameter `pi` is normalized.
#' @param x <double vector> : sample point
#' @param mu <double vector> : population means
#' @param sigma <double vector> : population sd
#' @param pi <double vector> : mixtured rate
#' @return The output is <data_frame> that is probability densities of each saple point at each component and mixture distribution.
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr as_tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr everything
#' @importFrom purrr map
#' @importFrom stringr str_c
#' @importFrom stats dnorm
#' @export
#' 
density_mixt_normal <- function(x, mu, sigma, pi){
  # 1. Preliminaries and Error Correction
  # number of componets calculated by dimension of parameter `mu`
  n_components <- length(mu)
  pi <- pi / sum(pi)
  # sanity check of each parameter's dimension
  if (length(sigma) != n_components | length(pi) != n_components){
    stop("Error : dimension of mu, sigma and pi must be the same.")
  }
  # 2. Main Calculation
  # calculation of probability density
  # pd_each_component : output matrix is probability densities of each component
  # pd : output vector is probability density
  pd_each_component <- as.list(x) %>%
    map(.f = ~dnorm(., mean = mu, sd = sigma)) %>%
    unlist() %>%
    matrix(byrow = TRUE, ncol = n_components)
  pd <- as.numeric(pd_each_component %*% pi)
  # 3. output
  obj <- as_tibble(pd_each_component)
  colnames(obj) <- str_c("component_", 1:n_components)
  obj <- obj %>%
    mutate(prob_density = pd, x = x) %>%
    select(x, prob_density, everything())
  return(obj)
}





#' Log Likelihood of 1-dim. Mixtures of Normal Distribution
#' 
#' @param x <double vector> : sample point
#' @param mu <double vector> : population means
#' @param sigma <double vector> : population sd
#' @param pi <double vector> : mixtured rate
#' @return The output is <double vector> that is loglikelihood of each point of x 
#' 
#' @importFrom dplyr %>%
#' @export
#' 
LL_mixt_normal <- function(x, mu, sigma, pi){
  # 1. Preliminaries and Error Correction
  # number of components calculated from dimension of mu
  n_components <- length(mu)
  # sanity check
  if (length(sigma) != n_components | length(pi) != n_components){
    stop("Error : dimension of mu, sigma and pi must be the same.")
  }
  # 2. Main Calculation
  # calculation of log likelihood
  obj <- sum(log(density_mixt_normal(x, mu = mu, sigma = sigma, pi = pi)$prob_density))
  # output
  return(obj)
}





#' EM Algorithm for 1-dim. Mixtures of Normal Distribution
#' 
#' @param x <double vector> : sample
#' @param max_iter <int scalar> : maximum number of iterations
#' @param tol <double scalar> : tolerance
#' @param init_mu <double vector> : initial values of population means
#' @param init_sigma <double vector> : initial values of population sds
#' @param init_pi <double vector> : initial values of population rate
#' @param history <logical double> : If you need history, you set the value TRUE. Defalut is FALSE.
#' @return The output is the set (MLE of params, LL, estimated component of each sample point, n_iter). If you set history == TRUE, you can get hisotry of params and LL.
#' 
#' @importFrom dplyr data_frame
#' @importFrom dplyr select
#' @importFrom dplyr bind_rows
#' @importFrom purrr map_chr
#' @export
#' 
EM_mixt_normal <- function(x, max_iter, tol, init_mu, init_sigma, init_pi, history = FALSE){
  # 1. Preliminalies and Error Correction
  # number of components and sample size
  sample_size <- length(x)    # sample size
  n_components <- length(init_mu)
  init_pi <- init_pi / sum(init_pi)
  # sanity check for dimension of parameters
  if (length(init_sigma) != n_components | length(init_pi) != n_components){
    stop("Error : dimension of mu, sigma and pi must be the same.")
  }
  # 2. Main Calculation
  # settings of initial values
  mu <- init_mu; sigma <- init_sigma; pi <- init_pi
  # history
  LL_history <- LL_mixt_normal(x = x, mu = mu, sigma = sigma, pi = pi)
  params_history <- data_frame(iter = 0, component = 1:n_components, mu = mu, sigma = sigma, pi = pi)
  # EM algorithm
  for (iter in 1:max_iter) {
    # E-step : gammma
    likelihood <- density_mixt_normal(x = x, mu = mu, sigma = sigma, pi = pi)
    gamma_each_component <- as.matrix(likelihood[ , 3:(2+n_components)] / likelihood$prob_density)
    colnames(gamma_each_component) <- as.character(1:n_components) %>%
      map_chr(.f = ~paste0("gamma", .))
    sum_gamma <- colSums(gamma_each_component)
    # M-step
    pi <- sum_gamma / sample_size    # MLE of pi
    mu <- as.vector(x %*% gamma_each_component) / sum_gamma    # MLE of mu
    ncol <- length(x); nrow <- n_components    # MLE of sigma
    x_matrix <- matrix(x, byrow = TRUE, nrow = nrow, ncol = ncol)
    mu_matrix <- matrix(mu, nrow = nrow, ncol = ncol)
    sigma <- sqrt(
      diag((x_matrix - mu_matrix) %*% (gamma_each_component * t(x_matrix - mu_matrix))) / sum_gamma
      )
    LL_history <- append(LL_history,        # recoding history...
                         LL_mixt_normal(x = x, mu = mu, sigma = sigma, pi = pi))
    params_history <- bind_rows(
      params_history,
      data_frame(iter = iter, component = 1:n_components, mu = mu, sigma = sigma, pi = pi)
    )
    # coveragement
    if(abs(LL_history[iter+1]-LL_history[iter]) < tol){
      break    
    }
  }
  # 3. output
  n_iter <- iter
  LL_history_df <- data_frame(iter = 0:n_iter, LL_history = LL_history)
  data_estimated_component <- data_frame(x = x, estimated_component = max.col(gamma_each_component))
  LL_last <- LL_history_df %>% filter(iter == n_iter) %>% .$LL_history
  params_last <- params_history %>% filter(iter == n_iter) %>% select(-iter)
  if (history == TRUE){
    obj <- list(params_history = params_history,
                LL_history = LL_history_df,
                params = params_last,
                LL = LL_last,
                estimated_component = data_estimated_component,
                n_iter = n_iter)
  } else{
    obj <- list(params = params_last, 
                LL = LL_last,
                estimated_component = data_estimated_component,
                n_iter = n_iter)
  }
  return(obj)
}
