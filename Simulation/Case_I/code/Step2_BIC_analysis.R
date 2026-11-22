##############################################################
#################     Simulation Case I    ###################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

library(invgamma)
library(truncnorm)
library(matrixStats)


getMode <- function(v){
  uniqv <- unique(v)
  return(as.numeric(uniqv[which.max(tabulate(match(v, uniqv)))]))
}


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

update_tau_xi_BIC <- function(K_value, p_tau, q_tau, xi_t){
  return(sqrt(rinvgamma(1, shape = p_tau + (K_value-1)*J/2, rate = q_tau + sum(xi_t[ ,-1]^2)/2)))
}

update_pi_BIC <- function(K_value, c, Z_t){
  if(K_value == 1) return(as.matrix(sapply(1:B, function(b){Dirichlet(c)})))
  else return(t(sapply(1:B, function(b){Dirichlet(sapply(1:K_value, function(k){sum(Z_t[[b]] == k)}) + c)})))
}

logBeMat <- function(y_mat, mu_mat, psi_mat){
  s1_mat = mu_mat * psi_mat
  s2_mat = psi_mat - mu_mat * psi_mat
  logBeta_mat = t(dbeta(t(y_mat), shape1 = s1_mat, shape2 = s2_mat, log = TRUE))
  return(logBeta_mat)
}

update_Z_BIC <- function(K_value, pi_t, alpha_t, xi_t, gamma_t, psi_t){
  temp = list()
  for(b in 1:B){
    log_prob = sapply(1:K_value, function(k){log(pi_t[b, k]) + rowSums(logBeMat(Y[[b]], 1/(1 + exp(-alpha_t-xi_t[, k]-gamma_t[b, ])), psi_t[b, ]))})
    Prob = exp(log_prob-rowMaxs(log_prob)) / rep(rowSums(exp(log_prob-rowMaxs(log_prob))),times=K_value)
    temp[[b]] = sapply(1:n_vec[b], function(i){sample(1:K_value, size = 1, prob = Prob[i, ])})
  }
  return(temp)
}

update_alpha_BIC <- function(eta_alpha, tau_alpha, phi_alpha, Z_t, alpha_t, xi_t, gamma_t, psi_t){
  ## update alpha by Metropolis-Hasting step
  alpha_temp = alpha_t
  proposal = rnorm(J, mean = alpha_t, sd = phi_alpha)
  r1 = dnorm(proposal, mean = eta_alpha, sd = tau_alpha, log = TRUE) - dnorm(alpha_t, mean = eta_alpha, sd = tau_alpha, log = TRUE)
  sumBe <- function(alpha){
    rowSums(sapply(1:B, function(b){rowSums(sapply(1:n_vec[b], function(i){logBe(Y[[b]][i, ], 1/(1 + exp(-alpha-xi_t[, Z_t[[b]][i]]-gamma_t[b, ])), psi_t[b, ])}))}))
  }
  r2 = sumBe(proposal) - sumBe(alpha_t)
  prob = pmin(exp(r1 + r2), 1)
  temp = runif(J)
  alpha_temp[temp <= prob] = proposal[temp <= prob]
  return(alpha_temp)
}

update_xi_BIC <- function(K_value, phi_xi, tau_xi_t, Z_t, alpha_t, xi_t, gamma_t, psi_t){
  ## update xi by Metropolis-Hasting step
  if(K_value == 1){
    return(xi_t)
  }
  xi_temp = xi_t
  for(k in 2:K_value){
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


update_gamma_BIC <- function(tau_gamma, phi_gamma, Z_t, alpha_t, xi_t, gamma_t, psi_t){
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


update_psi_BIC <- function(eta_psi, tau_psi, phi_psi, Z_t, alpha_t, xi_t, gamma_t, psi_t){
  ## update psi by Metropolis-Hasting step
  psi_temp = psi_t
  for(b in 1:B){
    proposal = rtruncnorm(J, a = 0, b = Inf, mean = psi_t[b, ], sd = phi_psi[b, ])
    r1 = dlnorm(proposal, meanlog = eta_psi, sdlog = tau_psi, log = TRUE) - dlnorm(psi_t[b, ], meanlog = eta_psi, sdlog = tau_psi, log = TRUE)
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



################# Gibbs sampling funciton ##################### 

GIBBS_BIC <- function(K_value){
  
  ## Hyper-parameters
  c <- 1
  eta_alpha <- 0
  tau_alpha <- 0.1
  p_tau <- 7
  q_tau <- 1
  tau_gamma <- 0.05
  eta_psi <- 0.65
  tau_psi <- 0.8
  phi_alpha <- rep(0.05, J)
  phi_xi <- matrix(0.1, nrow = J, ncol = K_value)
  phi_gamma <- matrix(0.03, nrow = B, ncol = J)
  phi_psi <- matrix(1.875, nrow = B, ncol = J)
  
  ## Initiation for parameters
  pi_0 <- matrix(NA, nrow = B, ncol = K_value)
  Z_0 <- list()
  for(b in 1:B){
    Z_0[[b]] <- sample(1:K_value, size = n_vec[b], replace = T)
    pi_0[b, ] <- as.data.frame(table(Z_0[[b]])) $ Freq / n_vec[b]
  }
  
  ### xi[j,1]==0, gamma[1,j]==0
  alpha_0 <- rnorm(J, mean = eta_alpha, sd = tau_alpha)
  tau_xi_0 <- rinvgamma(1, shape = p_tau, rate = q_tau)
  xi_0 <- matrix(c(rep(0,J), rnorm(J*(K_value-1), mean = 0, sd = tau_xi_0)), nrow = J, ncol = K_value)
  gamma_0 <- matrix(c(rep(0,J), rnorm((B-1)*J, mean = 0, sd = tau_gamma)), nrow = B, ncol = J, byrow = T)
  psi_0 <- matrix(rlnorm(B*J, meanlog = eta_psi, sdlog = tau_psi), nrow = B, ncol = J, byrow = T)
  
  
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
  pi_t_record = list()
  
  cat("  running the Gibbs sampler ...\n")
  iter = 1000
  iterhalf = iter / 2
  for(t in 1 : iter){
    
    if(t == 1) cat("######## GIBBS SAMPLING STARTS ########", "\n")
    tau_xi_t = update_tau_xi_BIC(K_value, p_tau, q_tau, xi_t)
    pi_t = update_pi_BIC(K_value, c, Z_t)
    Z_t = update_Z_BIC(K_value, pi_t, alpha_t, xi_t, gamma_t, psi_t)
    alpha_t = update_alpha_BIC(eta_alpha, tau_alpha, phi_alpha, Z_t, alpha_t, xi_t, gamma_t, psi_t)
    xi_t = update_xi_BIC(K_value, phi_xi, tau_xi_t, Z_t, alpha_t, xi_t, gamma_t, psi_t)
    gamma_t = update_gamma_BIC(tau_gamma, phi_gamma, Z_t, alpha_t, xi_t, gamma_t, psi_t)
    psi_t = update_psi_BIC(eta_psi, tau_psi, phi_psi, Z_t, alpha_t, xi_t, gamma_t, psi_t)
    
    alpha_t_record[[t]] = alpha_t
    xi_t_record[[t]] = xi_t
    gamma_t_record[[t]] = gamma_t
    psi_t_record[[t]] = psi_t
    pi_t_record[[t]] = pi_t
    
    if(t %% (iter / 10) ==0) cat("######## GIBBS SAMPLING for K =", K_value, "IS", t / iter * 100, "% COMPLETE ########", "\n")
  } 
  
  alpha_hat <- 0
  xi_hat <- 0
  gamma_hat <- 0
  psi_hat <- 0
  pi_hat <- 0
  
  for(t in (iter-iterhalf+1):iter){
    alpha_hat <- alpha_hat + alpha_t_record[[t]]
    xi_hat <- xi_hat + xi_t_record[[t]]
    gamma_hat <- gamma_hat + gamma_t_record[[t]]
    psi_hat <- psi_hat + psi_t_record[[t]]
    pi_hat <- pi_hat +pi_t_record[[t]]
  }
  alpha_hat <- alpha_hat / iterhalf
  xi_hat <- xi_hat / iterhalf
  gamma_hat <- gamma_hat / iterhalf
  psi_hat <- psi_hat / iterhalf
  pi_hat <- pi_hat / iterhalf
  
  BIC_value = -2 * sum(sapply(1:B, function(b){sum(sapply(1:n_vec[b], function(i){
    temp = sapply(1:K_value, function(k){
      log(pi_hat[b, k])+sum(logBe(Y[[b]][i, ], 1/(1 + exp(-alpha_hat-xi_hat[, k]-gamma_hat[b, ])), psi_hat[b, ]))
    })
    log(sum(exp(temp - max(temp))))+max(temp)
  }))})) + (K_value*J+(2*B-1)*J)*log(sum(n_vec)*J)
  
  return(BIC_value)
}



## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Read data Y
load("../input_data/Y.RData")
B <- length(Y)
J <- ncol(Y[[1]])
n_vec <- sapply(Y, nrow)
 

## Set seed for BUSbeta in simulation case I
seed = 2002
set.seed(seed)



## BIC analysis from K = 1 to 10
BIC_rec = c(NA)
for(k in 1 : 10){
  set.seed(seed)
  BIC_rec[k] = GIBBS_BIC(k)
}


# Save BIC values for K = 1 to 10
# save(BIC_rec, file = "../result_data/BIC_rec.RData")
# BIC_rec = c(31913.7, -243222.3, -349904.5, -353703.7, -425472.3, -326880.9, -398688.7, -366636.5, -355390.4, -339762.5)


## Save BIC plot
png(file="../figures/Figure1(d).png", width = 460, height = 400)
options(scipen = 999)
BIC_plot = plot(x = 1:10, y = BIC_rec, pch = 16, cex = 0.8, cex.axis = 20 / par("ps"))
lines(x = 1:10, y = BIC_rec)
dev.off()


