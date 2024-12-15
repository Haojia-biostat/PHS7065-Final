# define function to generate a single dataset
gen_data <- function(
    size,       # sample size
    mu, sigma,  # for X1
    p,          # for X2
    sigma_e # for epsilon
){
  X1 <- rnorm(size, mu, sigma)
  X2 <- rbinom(size, 1, p)
  Y <- 1 + 2*X1 + 3*X2 + rnorm(size, 0, sigma_e)
  return(data.frame(
    X1 = X1,
    X2 = X2,
    Y = Y
  ))
}

# density functions
w1 <- function(x1, mu_1, mu_0, sigma_0, truncation = NULL){
  w <- exp((1/(2*sigma_0^2))*(2*x1-mu_1-mu_0)*(mu_1-mu_0))
  if(!is.null(truncation)) w[w > truncation] <- truncation
  return(w)
}
w2 <- function(x1, sigma_1, mu_0, sigma_0, truncation = NULL){
  w <- (sigma_0/sigma_1)*exp(-0.5*((1/sigma_1^2)-(1/sigma_0^2))*(x1-mu_0)^2)
  if(!is.null(truncation)) w[w > truncation] <- truncation
  return(w)
}
w3 <- function(x2, p_1, p_0, truncation = NULL){
  w <- (p_1/p_0)^x2*((1-p_1)/(1-p_0))^(1-x2)
  if(!is.null(truncation)) w[w > truncation] <- truncation
  return(w)
}


# define function to combine labeled and unlabeled data
comb_data <- function(
    n, N,               # sample sizes
    mu_0, sigma_0, p_0, # X distribution in labeled data
    mu_1, sigma_1, p_1, # X distribution in unlabeled data
    sigma_e,            # noise level
    truncation = NULL   # truncation for weights
){
  # labeled data
  labeled_data <- gen_data(n, mu_0, sigma_0, p_0, sigma_e)
  # unlabeled data
  unlabeled_data <- gen_data(N, mu_1, sigma_1, p_1, sigma_e)
  
  # fit model
  fit <- lm(Y ~ X1 + X2, data = labeled_data)
  labeled_data$Y_hat <- predict(fit, labeled_data[,1:2])
  unlabeled_data$Y_hat <- predict(fit, unlabeled_data[,1:2])
  
  labeled_data$w <- 1
  # add weight
  if(mu_0 != mu_1){
    labeled_data$w <- w1(labeled_data$X1, mu_1 = mu_1, mu_0 = mu_0, sigma_0 = sigma_0, truncation = truncation)
  }
  if(sigma_0 != sigma_1){
    labeled_data$w <- w2(labeled_data$X1, sigma_1 = sigma_1, mu_0 = mu_0, sigma_0 = sigma_0, truncation = truncation)
  }
  if(p_0 != p_1){
    labeled_data$w <- w3(labeled_data$X2, p_1 = p_1, p_0 = p_0, truncation = truncation)
  }
  
  # generate output
  return(list(
    settings = c(
      n = n, N = N, 
      mu_0 = mu_0, sigma_0 = sigma_0, p_0 = p_0, 
      mu_1 = mu_1, sigma_1 = sigma_1, p_1 = p_1, 
      sigma_e = sigma_e
    ),
    labeled = labeled_data,
    unlabeled = unlabeled_data
  ))
}

# ---- ppi++ ----
lambda_ppi <- function(n, N, Y_n, Y_hat_n, Y_hat_N, w){
  cov <- cov(w*Y_n, w*Y_hat_n)
  var_N <- var(Y_hat_N)
  var_n <- var(w*Y_hat_n)
  lam <- cov/n/(var_N/N+var_n/n)
  return(lam)
}
pointest_ppi <- function(lam, n, N, Y_n, Y_hat_n, Y_hat_N, w){
  mean(w*Y_n) + lam*mean(unlist(Y_hat_N)) - lam*mean(Y_hat_n*w)
}
std_ppi <- function(lam, n, N, Y_n, Y_hat_n, Y_hat_N, w){
  Y_hat_N_var <- var(lam*Y_hat_N)/N
  rectifier_var <- var(w*(Y_n-lam*Y_hat_n))/n
  std <- sqrt(Y_hat_N_var+rectifier_var)
  return(std)
}
ci_ppi <- function(pointest, std, alpha) {
  z <- qnorm(1-alpha/2)
  ll <- pointest - z*std
  ul <- pointest + z*std
  return(c(ll, ul))
}
gen_res <- function(n, N, Y_n, Y_hat_n, Y_N, Y_hat_N, w, lam = NULL, alpha = 0.5) {
  if(is.null(lam)) {
    lam <- lambda_ppi(n, N, Y_n, Y_hat_n, Y_hat_N, w)
  }
  pointest <- pointest_ppi(lam, n, N, Y_n, Y_hat_n, Y_hat_N, w)
  std <- std_ppi(lam, n, N, Y_n, Y_hat_n, Y_hat_N, w)
  ci <- ci_ppi(pointest, std, alpha)
  return(c(
    lam = lam,
    truemean = mean(Y_N),
    pointest = pointest,
    std = std,
    ll = ci[1],
    ul = ci[2],
    ci_width = ci[2] - ci[1],
    coverage = (ci[1] <= mean(Y_N) & mean(Y_N) <= ci[2])
  ))
}

# ---- cross ppi++ ----
lambda_ppi_cross <- function(n, N, Y_n, Y_hat_n, Y_hat_N, w, K, fold){
  cov <- sum(sapply(1:K, \(j) cov(w[fold == j]*Y_n[fold == j], w[fold == j]*Y_hat_n[fold == j])) * tapply(fold, fold, length)) / n
  var_N <- sum(apply(Y_hat_N, 2, var)) / (K^2)
  var_n <- sum(sapply(1:K, \(j) var(w[fold == j]*Y_hat_n[fold == j])) * tapply(fold, fold, length)) / n
  lam <- cov/n/(var_N/N+var_n/n)
  return(lam)
}
std_ppi_cross <- function(lam, n, N, Y_n, Y_hat_n, Y_hat_N, w, K, fold){
  Y_hat_N_var <- sum(apply(Y_hat_N, 2, var)) *lam^2 / (K^2*N)
  rectifier_var <- sum(sapply(1:K, \(j) var(w[fold == j]*(Y_n[fold == j] - lam*Y_hat_n[fold == j]))) * tapply(fold, fold, length)) / (n^2)
  std <- sqrt(Y_hat_N_var+rectifier_var)
  return(std)
}
gen_res_cross <- function(n, N, Y_n, Y_hat_n, Y_N, Y_hat_N, w, K, fold, lam = NULL, alpha = 0.5) {
  if(is.null(lam)) {
    lam <- lambda_ppi_cross(n, N, Y_n, Y_hat_n, Y_hat_N, w, K, fold)
  }
  pointest <- pointest_ppi(lam, n, N, Y_n, Y_hat_n, Y_hat_N, w)
  std <- std_ppi_cross(lam, n, N, Y_n, Y_hat_n, Y_hat_N, w, K, fold)
  ci <- ci_ppi(pointest, std, alpha)
  return(c(
    lam = lam,
    truemean = mean(Y_N),
    pointest = pointest,
    std = std,
    ll = ci[1],
    ul = ci[2],
    ci_width = ci[2] - ci[1],
    coverage = (ci[1] <= mean(Y_N) & mean(Y_N) <= ci[2])
  ))
}
