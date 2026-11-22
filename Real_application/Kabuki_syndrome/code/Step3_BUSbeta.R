##############################################################
################    Kabuki syndrome data    ##################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

library(invgamma)
library(truncnorm)
library(aricode)


getMode <- function(v){
  uniqv <- unique(v)
  return(as.numeric(uniqv[which.max(tabulate(match(v, uniqv)))]))
}

logit <- function(x){ return(log(x / (1 - x))) } 

logBe <- function(y_vec, mu_vec, psi_vec){
  s1 = mu_vec * psi_vec
  s2 = psi_vec - mu_vec * psi_vec
  logBeta_rvs = dbeta(y_vec, shape1 = s1, shape2 = s2, log = TRUE)
  return(logBeta_rvs)
}

rBe <- function(mu, psi){
  s1 = mu * psi
  s2 = psi - mu * psi
  beta_rv = rbeta(1, shape1 = s1, shape2 = s2)
  return(beta_rv)
}

Dirichlet <- function(c_vec){
  n = length(c_vec)
  gamma_rvs = sapply(1:n, function(i){rgamma(1, shape = c_vec[i], rate = 1)})
  return(gamma_rvs/sum(gamma_rvs))
} 


#################### Full Conditionals ######################## 

update_tau_xi <- function(xi_t){
  return(sqrt(rinvgamma(1, shape = p_tau + (K-1)*J/2, rate = q_tau + sum(xi_t[ ,-1]^2)/2)))
}

update_pi <- function(Z_t){
  return(t(sapply(1:B, function(b){Dirichlet(sapply(1:K, function(k){sum(Z_t[[b]] == k)}) + c)})))
}

update_Z <- function(pi_t, alpha_t, xi_t, gamma_t, psi_t){
  temp = lapply(1:B, function(b){
    sapply(1:n_vec[b], function(i){
      log_prob = sapply(1:K, function(k){log(pi_t[b, k]) + sum(logBe(Y[[b]][i, ], 1/(1 + exp(-alpha_t-xi_t[, k]-gamma_t[b, ])), psi_t[b, ]))})
      Prob = exp(log_prob-max(log_prob)) / sum(exp(log_prob-max(log_prob)))
      sample(1:K, size = 1, prob = Prob)
    })
  })
  return(temp)
}

update_alpha <- function(Z_t, alpha_t, xi_t, gamma_t, psi_t){
  ## update alpha by Metropolis-Hasting step
  alpha_temp = alpha_t
  proposal = rnorm(J, mean = alpha_t, sd = phi_alpha)
  r1 = dnorm(proposal, mean = ita_alpha, sd = tau_alpha, log = TRUE) - dnorm(alpha_t, mean = ita_alpha, sd = tau_alpha, log = TRUE)
  sumBe <- function(alpha){
    rowSums(sapply(1:B, function(b){rowSums(sapply(1:n_vec[b], function(i){logBe(Y[[b]][i, ], 1/(1 + exp(-alpha-xi_t[, Z_t[[b]][i]]-gamma_t[b, ])), psi_t[b, ])}))}))
  }
  r2 = sumBe(proposal) - sumBe(alpha_t)
  prob = pmin(exp(r1 + r2), 1)
  temp = runif(J)
  alpha_temp[temp <= prob] = proposal[temp <= prob]
  return(alpha_temp)
}

update_xi <- function(tau_xi_t, Z_t, alpha_t, xi_t, gamma_t, psi_t){
  ## update xi by Metropolis-Hasting step
  xi_temp = xi_t
  for(k in 2:K){
    proposal = rnorm(J, mean = xi_t[, k], sd = phi_xi[, k])
    r1 = dnorm(proposal, mean = 0, sd = tau_xi_t, log = TRUE) - dnorm(xi_t[, k], mean = 0, sd = tau_xi_t, log = TRUE)
    sumBe <- function(xi){
      rowSums(sapply(1:B, function(b){
        I_seq = which(Z_t[[b]]==k)
        if(length(I_seq) == 0) rep(0, J)
        else rowSums(sapply(I_seq, function(i){logBe(Y[[b]][i, ], 1/(1 + exp(-alpha_t-xi-gamma_t[b, ])), psi_t[b, ])}))
      }))
    }
    r2 = sumBe(proposal) - sumBe(xi_t[, k])
    prob = pmin(exp(r1 + r2), 1)
    temp = runif(J)
    xi_temp[temp <= prob, k] = proposal[temp <= prob]
  }
  return(xi_temp)
}


update_gamma <- function(Z_t, alpha_t, xi_t, gamma_t, psi_t){
  ## update gamma by Metropolis-Hasting step
  gamma_temp = gamma_t
  for(b in 2:B){
    proposal = rnorm(J, mean = gamma_t[b, ], sd = phi_gamma[b, ])
    r1 = dnorm(proposal, mean = 0, sd = tau_gamma, log = TRUE) - dnorm(gamma_t[b, ], mean = 0, sd = tau_gamma, log = TRUE)
    sumBe <- function(gamma){
      rowSums(sapply(1:n_vec[b], function(i){logBe(Y[[b]][i, ], 1/(1 + exp(-alpha_t-xi_t[, Z_t[[b]][i]]-gamma)), psi_t[b, ])}))
    }
    r2 = sumBe(proposal) - sumBe(gamma_t[b, ])
    prob = pmin(exp(r1 + r2), 1)
    temp = runif(J)
    gamma_temp[b, temp <= prob] = proposal[temp <= prob]
  }
  return(gamma_temp)
}


update_psi <- function(Z_t, alpha_t, xi_t, gamma_t, psi_t){
  ## update psi by Metropolis-Hasting step
  psi_temp = psi_t
  for(b in 1:B){
    proposal = rtruncnorm(J, a = 0, b = Inf, mean = psi_t[b, ], sd = phi_psi[b, ])
    r1 = dlnorm(proposal, meanlog = ita_psi, sdlog = tau_psi, log = TRUE) - dlnorm(psi_t[b, ], meanlog = ita_psi, sdlog = tau_psi, log = TRUE)
    sumBe <- function(psi){
      rowSums(sapply(1:n_vec[b], function(i){logBe(Y[[b]][i, ], 1/(1 + exp(-alpha_t-xi_t[, Z_t[[b]][i]]-gamma_t[b, ])), psi)}))
    }
    r2 = sumBe(proposal) - sumBe(psi_t[b, ])
    r3 = log(dtruncnorm(psi_t[b, ], a=0, b=Inf, mean = proposal, sd = phi_psi[b, ])) - log(dtruncnorm(proposal, a=0, b=Inf, mean = psi_t[b, ], sd = phi_psi[b, ]))
    prob = pmin(exp(r1 + r2 + r3), 1)
    temp = runif(J)
    psi_temp[b, temp <= prob] = proposal[temp <= prob]
  }
  return(psi_temp)
}



## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Read data Y_preprocessed
Y_preprocessed <- read.table("../input_data/Y_preprocessed.txt", sep = "\t", header = T)
Y_preprocessed <- as.data.frame(Y_preprocessed)

B <- length(unique(Y_preprocessed $ batch))
J <- ncol(Y_preprocessed) - 2
K <- 2 ### From BIC analysis
n_vec <- as.vector(table(Y_preprocessed $ batch))


## Data Y list
Y <- list()
for(b in 1:B){
  ind = which(Y_preprocessed $ batch == b)
  Y[[b]] = as.matrix(Y_preprocessed[ind, c(-J-1, -J-2)])
}


## Read true subgroups Z_real
load("../input_data/Z_real.RData")


## Set seed for BUSbeta in application to the Kabuki syndrome dataset
seed = 11
set.seed(seed)


## Hyper-parameters
c <- 1
ita_alpha <- 0
tau_alpha <- 0.1
p_tau <- 7
q_tau <- 1
tau_gamma <- 0.05
ita_psi <- 0.65
tau_psi <- 0.8
phi_alpha <- rep(0.05, J)
phi_xi <- matrix(0.1, nrow = J, ncol = K)
phi_gamma <- matrix(0.03, nrow = B, ncol = J)
phi_psi <- matrix(1.875, nrow = B, ncol = J)


## Initiation for parameters
pi_0 <- matrix(NA, nrow = B, ncol = K)

Z_0 <- list()
for(b in 1:B){
  Z_0[[b]] <- sample(1:K, size = n_vec[b], replace = T)
  pi_0[b, ] <- as.data.frame(table(Z_0[[b]])) $ Freq / n_vec[b]
}

### xi[j,1]==0, gamma[1,j]==0
alpha_0 <- rnorm(J, mean = ita_alpha, sd = tau_alpha)
tau_xi_0 <- rinvgamma(1, shape = p_tau, rate = q_tau)
xi_0 <- matrix(c(rep(0,J), rnorm(J*(K-1), mean = 0, sd = tau_xi_0)), nrow = J, ncol = K)
gamma_0 <- matrix(c(rep(0,J), rnorm((B-1)*J, mean = 0, sd = tau_gamma)), nrow = B, ncol = J, byrow = T)
psi_0 <- matrix(rlnorm(B*J, meanlog = ita_psi, sdlog = tau_psi), nrow = B, ncol = J, byrow = T)

pi_t = pi_0
Z_t = Z_0
alpha_t = alpha_0
xi_t = xi_0
gamma_t = gamma_0
psi_t = psi_0
tau_xi_t = tau_xi_0

alpha_t_record = list()
xi_t_record = list()
gamma_t_record = list()
psi_t_record = list()
Z_t_record = list()


## Gibbs sampling
t1 = Sys.time()
cat("  running the Gibbs sampler ...\n")
iter = 3000
iterhalf = iter / 2
for(t in 1 : iter){
  if(t == 1) cat("######## GIBBS SAMPLING STARTS ########", "\n")
  ### tau_xi
  tau_xi_t = update_tau_xi(xi_t)
  ### pi
  pi_t = update_pi(Z_t)
  ### Z
  Z_t = update_Z(pi_t, alpha_t, xi_t, gamma_t, psi_t)
  ### alpha
  alpha_t = update_alpha(Z_t, alpha_t, xi_t, gamma_t, psi_t)
  ### xi
  xi_t = update_xi(tau_xi_t, Z_t, alpha_t, xi_t, gamma_t, psi_t)
  ### gamma
  gamma_t = update_gamma(Z_t, alpha_t, xi_t, gamma_t, psi_t)
  ### psi
  psi_t = update_psi(Z_t, alpha_t, xi_t, gamma_t, psi_t)
  
  alpha_t_record[[t]] = alpha_t
  xi_t_record[[t]] = xi_t
  gamma_t_record[[t]] = gamma_t
  psi_t_record[[t]] = psi_t
  Z_t_record[[t]] = Z_t
  
  if(t %% (iter / 10) == 0) cat("######## GIBBS SAMPLING IS ", t / iter * 100, "% COMPLETE ########", "\n")
} 
t2 = Sys.time()

cat(paste0("  The Gibbs sampler takes: ", round(difftime( t2, t1, units = "mins"), 3), " mins", "\n"))



## Estimated parameters
alpha_hat <- 0
xi_hat <- 0
gamma_hat <- 0
psi_hat <- 0
Z_BUSbeta <- list()
mu_hat <- list()
sigma_hat <- list()
sigma_weighted_hat <- list()


for(t in (iter-iterhalf+1):iter){
  alpha_hat <- alpha_hat + alpha_t_record[[t]]
  xi_hat <- xi_hat + xi_t_record[[t]]
  gamma_hat <- gamma_hat + gamma_t_record[[t]]
  psi_hat <- psi_hat + psi_t_record[[t]]
}
alpha_hat <- alpha_hat / iterhalf
xi_hat <- xi_hat / iterhalf
gamma_hat <- gamma_hat / iterhalf
psi_hat <- psi_hat / iterhalf


## BUSbeta clustering label
for(b in 1:B){
  Z_BUSbeta[[b]] = rep(NA, times = n_vec[b])
  for(i in 1:n_vec[b]){
    Z_BUSbeta[[b]][i] = getMode(sapply((iter-iterhalf+1):iter, function(t){Z_t_record[[t]][[b]][i]}))
  }
}

## Compute ARI
ARI_BUSbeta <- data.frame(Batch1 = round(ARI(Z_real[[1]], Z_BUSbeta[[1]]), 3), 
                          Batch2 = round(ARI(Z_real[[2]], Z_BUSbeta[[2]]), 3), 
                          Overall = round(ARI(unlist(Z_real), unlist(Z_BUSbeta)), 3))

## Save ARI
write.csv(ARI_BUSbeta, row.names = "BUSbeta ARI", file = "../result_data/ARI_BUSbeta.csv")



###################### batch correction #########################

Y_BUSbeta <- list()
Y_BUSbeta[[1]] = Y[[1]]

for(b in 1:B){
  mu_hat[[b]] = 1 / (1 + exp(-matrix(rep(alpha_hat, K), nrow =J, ncol = K) - xi_hat - matrix(rep(gamma_hat[b, ], K), nrow = J, ncol = K)))
  sigma_hat[[b]] = sqrt(mu_hat[[b]] * (1 - mu_hat[[b]]) / (psi_hat[b, ] + 1)) # list b with J * K matrix
  sigma_weighted_hat[[b]] = sqrt(apply(sapply(1:K, function(k){sigma_hat[[b]][, k]^2 * length(which(Z_BUSbeta[[b]]==k))}), MARGIN = 1, sum) / n_vec[b])
}

X1 <- list()
for(b in 1:B){
  X1[[b]] = matrix(0, nrow = n_vec[b], ncol = K)
  for(i in 1:n_vec[b]){
    X1[[b]][i, Z_BUSbeta[[b]][i]] = 1
  }
}

X2 <- matrix(0, nrow = sum(n_vec), ncol = B-1)
for(b in 2:B) X2[(sum(n_vec[1:(b-1)])+1):sum(n_vec[1:b]), b-1] = rep(1, times = n_vec[b])


Y_combined = do.call(rbind, Y)
design_X = cbind(rep(1, times = sum(n_vec)), do.call(rbind, X1)[, -1], X2) ### n * (K+B-1) matrix
design_X_col0_ind = which(colSums(design_X != 0) == 0) ### l all-zero columns
l = length(design_X_col0_ind)
design_X_cleaned0 <- design_X[, colSums(design_X != 0) > 0] ### n * (K+B-1-l) matrix


fit <- matrix(0, nrow = K+B-1-l, ncol = J) ### (K+B-1-l) * J matrix
for(j in 1:J) fit[, j] = lm.fit(design_X_cleaned0, Y_combined[, j]) $ coefficients
omega_hat = fit[1, ] ### J dimensinal vector
m_hat = fit[2:(K-l), ] ### (K-1-l) * J matrix
if(l > 0){
  temp = matrix(0, nrow = K, ncol = J)
  temp[-c(1, design_X_col0_ind), ] = m_hat
  m_hat = temp ### K * J matrix
} else m_hat = rbind(rep(0, times = J), m_hat) ### K * J matrix


eta_hat = fit[(K-l+1):(K-l+B-1), ] ### (B-1) * J matrix
eta_hat = rbind(rep(0, times = J), eta_hat) ### B * J matrix

temp_eps = 0
for(b in 1:B){
  for(i in 1:n_vec[b]){
    temp_eps = temp_eps + (Y[[b]][i, ] - omega_hat - X1[[b]][i, ] %*% m_hat - eta_hat[b, ])^2
  }
}
eps_hat = sqrt(temp_eps / sum(n_vec)) ### J dimensinal vector

for(b in 2:B){
  Y_BUSbeta[[b]] = matrix(NA, nrow = n_vec[b], ncol = J)
  for(i in 1:n_vec[b]){
    Y_BUSbeta[[b]][i, ] = ((Y[[b]][i, ] - omega_hat - X1[[b]][i, ] %*% m_hat - eta_hat[b, ]) / eps_hat + eta_hat[b, ] - mu_hat[[b]][, Z_BUSbeta[[b]][i]] + 1 / (1 + exp(-alpha_hat - xi_hat[, Z_BUSbeta[[b]][i]]))) * eps_hat / (sigma_weighted_hat[[b]] / sigma_weighted_hat[[1]]) + omega_hat + X1[[b]][i, ] %*% m_hat
    ## [0,1] truncation
    less0 = which(Y_BUSbeta[[b]][i, ] < 0)
    more1 = which(Y_BUSbeta[[b]][i, ] > 1)
    Y_BUSbeta[[b]][i, less0] = rep(0, times = length(less0))
    Y_BUSbeta[[b]][i, more1] = rep(1, times = length(more1))
  }
}


## Save Corrected data and clustering label
save(Y_BUSbeta, file = "../result_data/Y_BUSbeta.RData")
save(Z_BUSbeta, file = "../result_data/Z_BUSbeta.RData")







