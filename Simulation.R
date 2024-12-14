# Setup ####
N <- 10000
n <- c(100, 200, 500, 1000, 2000, 5000)
sigma_e <- c(0.1, 1, 2, 5)
mu_0 <- 0
sigma_0 <- 1
p_0 <- 0.5
mu_1 <- c(0.2, 1, 2)
sigma_1 <- c(0.1, 0.5, 2)
p_1 <- c(0.55, 0.7, 0.85)
K <- 5

w1 <- function(x1, mu_1, mu_0, sigma_0){
  exp((1/(2*sigma_0^2))*(2*x1-mu_1-mu_0)*(mu_1-mu_0))
}
w2 <- function(x1, sigma_1, mu_0, sigma_0){
  (sigma_0/sigma_1)*exp(-0.5*((1/sigma_1^2)-(1/sigma_0^2))*(x1-mu_0)^2)
}
w3 <- function(x2, p_1, p_0){
  (p_1/p_0)^x2*((1-p_1)/(1-p_0))^(1-x2)
}

# data simulation ####
gen_data <- function(
    size, 
    mu, sigma, # for X1
    p, # for X2
    noise_level # vector, for Y = Y_hat + e
){
  X1 <- rnorm(size, mu, sigma)
  X2 <- rbinom(size, 1, p_0)
  Y <- 1 + 2*X1 + 3*X2 + rnorm(size, 0, noise_level)
  return(data.frame(
    X1 = X1,
    X2 = X2,
    Y = Y
  ))
}

# unlabeled data
# parameter combinations
param_comb <- rbind(
  # no covariate shift
  # c(mu = mu_0, sigma = sigma_0, p = p_0),
  # scenario 1: covariate shift in X1 with mu changed
  data.frame(mu = mu_1, sigma = sigma_0, p = p_0, scenario = "Scenario 1"),
  # scenario 2: covariate shift in X1 with sigma changed
  data.frame(mu = mu_0, sigma = sigma_1, p = p_0, scenario = "Scenario 2"),
  # scenario 3: covariate shift in X2 with p changed
  data.frame(mu = mu_0, sigma = sigma_0, p = p_1, scenario = "Scenario 3")
)

comb_data <- function(
    n, N,
    mu_0, sigma_0, p_0,
    mu_1, sigma_1, p_1,
    sigma_e
) {
  # labeled data
  labeled_data <- gen_data(n, mu_0, sigma_0, p_0, sigma_e)
  # unlabeled data
  unlabeled_data <- gen_data(N, mu_1, sigma_1, p_1, sigma_e)
  
  # fit model
  fit <- lm(Y ~ X1 + X2, data = labeled_data)
  labeled_data$Y_hat <- predict(fit, labeled_data[,1:2])
  unlabeled_data$Y_hat <- predict(fit, unlabeled_data[,1:2])
  
  # add weight
  if(mu_0 != mu_1){
    labeled_data$w <- w1(labeled_data$X1, mu_1 = mu_1, mu_0 = mu_0, sigma_0 = sigma_0)
  }
  if(sigma_0 != sigma_1){
    labeled_data$w <- w2(labeled_data$X1, sigma_1 = sigma_1, mu_0 = mu_0, sigma_0 = sigma_0)
  }
  if(p_0 != p_1){
    labeled_data$w <- w3(labeled_data$X2, p_1 = p_1, p_0 = p_0)
  }
  
  return(list(
    settings = c(n = n, N = N, mu_0 = mu_0, sigma_0 = sigma_0, p_0 = p_0, mu_1 = mu_1, sigma_1 = sigma_1, p_1 = p_1, sigma_e = sigma_e),
    labeled = labeled_data,
    unlabeled = unlabeled_data
  ))
}

settings <- merge(expand.grid(n = n, N = N, sigma_e = sigma_e), param_comb, all = T)

# ppi++ ####
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

# cross ppi++ ####
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

results <- vector(mode = "list", length = 100)
for(i in 1:100) {
  combined_data_list_no_cross <- mapply(
    comb_data,
    settings$n,
    settings$N,
    mu_0,
    sigma_0,
    p_0,
    settings$mu,
    settings$sigma,
    settings$p,
    settings$sigma_e,
    SIMPLIFY = F
  )    
  
  lamopt_ppi_res <- sapply(
    combined_data_list_no_cross,
    \(x) {
      res <- gen_res(
        n = x$settings["n"],
        N = x$settings["N"],
        Y_n = x$labeled$Y,
        Y_hat_n = x$labeled$Y_hat,
        Y_N = x$unlabeled$Y,
        Y_hat_N = x$unlabeled$Y_hat,
        w = x$labeled$w,
        alpha = 0.05
      )
      return(res)
    }
  )
  lamopt_nodensratio_ppi_res <- sapply(
    combined_data_list_no_cross,
    \(x) {
      res <- gen_res(
        n = x$settings["n"],
        N = x$settings["N"],
        Y_n = x$labeled$Y,
        Y_hat_n = x$labeled$Y_hat,
        Y_N = x$unlabeled$Y,
        Y_hat_N = x$unlabeled$Y_hat,
        w = 1,
        alpha = 0.05
      )
      return(res)
    }
  )
  lam1_ppi_res <- sapply(
    combined_data_list_no_cross,
    \(x) {
      res <- gen_res(
        n = x$settings["n"],
        N = x$settings["N"],
        Y_n = x$labeled$Y,
        Y_hat_n = x$labeled$Y_hat,
        Y_N = x$unlabeled$Y,
        Y_hat_N = x$unlabeled$Y_hat,
        w = x$labeled$w,
        lam = 1,
        alpha = 0.05
      )
      return(res)
    }
  )
  lam0_ppi_res <- sapply(
    combined_data_list_no_cross,
    \(x) {
      res <- gen_res(
        n = x$settings["n"],
        N = x$settings["N"],
        Y_n = x$labeled$Y,
        Y_hat_n = x$labeled$Y_hat,
        Y_N = x$unlabeled$Y,
        Y_hat_N = x$unlabeled$Y_hat,
        w = x$labeled$w,
        lam = 0,
        alpha = 0.05
      )
      return(res)
    }
  )
  
  # PPI++
  res_summary_lamopt <- data.frame(
    settings,
    t(lamopt_ppi_res) |> as.data.frame()
  )
  # PPI++, no density ratio
  res_summary_lamopt_nodensratio <- data.frame(
    settings,
    t(lamopt_nodensratio_ppi_res) |> as.data.frame()
  )
  # standard PPI
  res_summary_lam1 <- data.frame(
    settings,
    t(lam1_ppi_res) |> as.data.frame()
  )
  # CLT
  res_summary_lam0 <- data.frame(
    settings,
    t(lam0_ppi_res) |> as.data.frame()
  )
  
  # add K to labeled data
  combined_data_list_cross <- lapply(
    combined_data_list_no_cross,
    \(x) {
      x$labeled$fold <- sample(1:K, nrow(x$labeled), replace = TRUE)
      return(x)
    }
  )
  # predict Y_hat_N by K models
  combined_data_list_cross <- lapply(
    combined_data_list_cross,
    \(x) {
      labeled <- x$labeled
      fit_list <- lapply(
        1:K,
        \(j) lm(Y ~ X1+X2, data = labeled, subset = fold == j)
      )
      k_coef <- sapply(fit_list, unlist(coef)) |> t() |> as.data.frame()
      labeled$Y_hat <- apply(labeled, 1, \(obs) {
        coef <- k_coef[obs["fold"],] |> unlist()
        coef[1] + coef[2]*obs["X1"] + coef[3]*obs["X2"]
      })
      
      unlabeled <- x$unlabeled
      Y_hat_N <- data.frame(lapply(
        fit_list,
        \(fit) predict(fit, newdata = unlabeled[,1:2])
      ))
      
      x <- list(
        settings = x[[1]], 
        labeled = labeled, 
        unlabeled = x[[3]], 
        coef = k_coef, 
        Y_hat_N = Y_hat_N
      )
      
      return(x)
    }
  )
  
  lamopt_ppi_res_cross <- sapply(
    combined_data_list_cross,
    \(x) {
      res <- gen_res_cross(
        n = x$settings["n"],
        N = x$settings["N"],
        Y_n = x$labeled$Y,
        Y_hat_n = x$labeled$Y_hat,
        Y_N = x$unlabeled$Y,
        Y_hat_N = x$Y_hat_N,
        w = x$labeled$w,
        K = K,
        fold = x$labeled$fold,
        alpha = 0.05
      )
      return(res)
    }
  )
  
  res_summary_lamopt_cross <- data.frame(
    settings,
    t(lamopt_ppi_res_cross) |> as.data.frame()
  )
  
  colnames(res_summary_lam0) <- colnames(res_summary_lam1) <- colnames(res_summary_lamopt) <- colnames(res_summary_lamopt_nodensratio) <- colnames(res_summary_lamopt_cross) <- c("n", "N", "noise_level", "mu", "sigma", "p", "scenario", "lam", "truemean", "pointest", "std", "ll", "ul", "ci_width", "coverage")
  
  # different n, with everything else the same
  res_summary_all <- rbind(
    cbind(res_summary_lam0, model = "classical CLT"),
    cbind(res_summary_lam1, model = "standard PPI"),
    cbind(res_summary_lamopt, model = "PPI++"),
    cbind(res_summary_lamopt_nodensratio, model = "PPI++, no density ratio"),
    cbind(res_summary_lamopt_cross, model = "Cross PPI++")
  )
  
  results[[i]] <- res_summary_all
  print(i)
}

saveRDS(settings, file = "Final/simulation_settings.rds")
saveRDS(results, file = "Final/simulation_results.rds")
