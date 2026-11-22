##############################################################
#################     Simulation Case I    ###################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

library(invgamma)
library(truncnorm)
library(coda)


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


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Read data Y
load("../input_data/Y.RData")
load("../input_data/Z_real.RData")
B <- length(Y)
J <- ncol(Y[[1]])
K <- 5 ### From BIC analysis
n_vec <- sapply(Y, nrow)


## Set seed for all Markov chains of BUSbeta in simulation case I
seed = 2002
set.seed(seed)

## Number of Markov chains
n_chain = 4


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
phi_xi <- matrix(0.1, nrow = J, ncol = K)
phi_gamma <- matrix(0.03, nrow = B, ncol = J)
phi_psi <- matrix(1.875, nrow = B, ncol = J)


## Perturb initial values of selected parameters to form n_chain Markov chains
eta_alpha_chain = c(
  eta_alpha,
  eta_alpha + 0.01,
  eta_alpha - 0.01,
  eta_alpha + 0.02
)

tau_alpha_chain = c(
  tau_alpha,
  tau_alpha * 1.02,
  tau_alpha * 0.98,
  tau_alpha * 1.05
)

tau_gamma_chain = c(
  tau_gamma,
  tau_gamma * 1.02,
  tau_gamma * 0.98,
  tau_gamma * 1.05
)

eta_psi_chain = c(
  eta_psi,
  eta_psi + 0.01,
  eta_psi - 0.01,
  eta_psi + 0.02
)

tau_psi_chain = c(
  tau_psi,
  tau_psi * 1.02,
  tau_psi * 0.98,
  tau_psi * 1.05
)



##############################################################
## Run all n_chain Markov chains
##############################################################

alpha_record_chain = list()
gamma_record_chain = list()
psi_record_chain = list()


for(chain in 1:n_chain){
  
  set.seed(seed)
  
  ## Initiation for parameters
  pi_0 <- matrix(NA, nrow = B, ncol = K)
  
  Z_0 <- list()
  for(b in 1:B){
    Z_0[[b]] <- sample(1:K, size = n_vec[b], replace = T)
    pi_0[b, ] <- as.data.frame(table(Z_0[[b]])) $ Freq / n_vec[b]
  }
  
  tau_xi_0 <- rinvgamma(1, shape = p_tau, rate = q_tau)
  xi_0 <- matrix(c(rep(0,J), rnorm(J*(K-1), mean = 0, sd = tau_xi_0)), nrow = J, ncol = K)
  
  alpha_0 <- rnorm(J, mean = eta_alpha_chain[chain], sd = tau_alpha_chain[chain])
  
  gamma_0 <- matrix(c(rep(0,J), rnorm((B-1)*J, mean = 0, sd = tau_gamma_chain[chain])), nrow = B, ncol = J, byrow = T)
  
  psi_0 <- matrix(rlnorm(B*J, meanlog = eta_psi_chain[chain], sdlog = tau_psi_chain[chain]), nrow = B, ncol = J, byrow = T)
  
  ## Initialize the selected parameters 
  alpha_t = alpha_0
  gamma_t = gamma_0
  psi_t = psi_0
  
  alpha_t_record = list()
  gamma_t_record = list()
  psi_t_record = list()
  
  ## Initialize other parameters and keep original initialization
  pi_t = pi_0
  Z_t = Z_0
  xi_t = xi_0
  tau_xi_t = tau_xi_0
  
  
  ## Gibbs sampling
  cat("  running the Gibbs sampler for the chain", chain, "...\n")
  iter = 1000
  iterhalf = iter / 2
  for(t in 1:iter){
    if(t == 1) cat("######## GIBBS SAMPLING STARTS ########", "\n")
    
    tau_xi_t = update_tau_xi(xi_t)
    
    pi_t = update_pi(Z_t)
    
    Z_t = update_Z(pi_t, alpha_t, xi_t, gamma_t, psi_t)
    
    alpha_t = update_alpha(Z_t, alpha_t, xi_t, gamma_t, psi_t)
    
    xi_t = update_xi(tau_xi_t, Z_t, alpha_t, xi_t, gamma_t, psi_t)
    
    gamma_t = update_gamma(Z_t,alpha_t, xi_t, gamma_t, psi_t)
    
    psi_t = update_psi(Z_t, alpha_t, xi_t, gamma_t, psi_t)
    
    
    alpha_t_record[[t]] = alpha_t
    gamma_t_record[[t]] = gamma_t
    psi_t_record[[t]] = psi_t
    
    if(t %% (iter / 10) == 0) cat("######## GIBBS SAMPLING IS ", t / iter * 100, "% COMPLETE ########", "\n")
    
  }
  
  
  alpha_record_chain[[chain]] = alpha_t_record
  gamma_record_chain[[chain]] = gamma_t_record
  psi_record_chain[[chain]] = psi_t_record
  
}



##############################################################
## Potential scale reduction factor: R_hat
##############################################################

## Function to calculate the percent quantile of R_hat
get_percent_Rhat <- function(record_chain, end_iter, type, percent){
  
  mcmc_chain = list()
  
  for(chain in 1:n_chain){
    temp = record_chain[[chain]][(iterhalf+1):end_iter]
    if(type == "alpha"){
      ## iteration * J
      temp = do.call(rbind, temp)
    }else{
      ## iteration * (B*J)
      temp = do.call(rbind, lapply(temp, c))
    }
    mcmc_chain[[chain]] = mcmc(temp)
  }
  
  ## R_hat for each parameter
  result = gelman.diag(
    mcmc.list(mcmc_chain),
    autoburnin = FALSE,
    multivariate = FALSE
  )$psrf[,1]
  
  if(type == "alpha"){
    ## Percent quantile among J alpha parameters
    return(as.numeric(quantile(result, 
                               probs = percent,
                               na.rm = TRUE)))
  }else{
    ## Reshape R_hat of B*J parameters
    result = matrix(result,
                    nrow = B,
                    ncol = J)
    ## Percent quantile within each batch
    return(apply(result,
                 1,
                 function(x) quantile(x,
                                      probs = percent,
                                      na.rm = TRUE)))
  }
}



## Calculate the 95th percentile of R_hat at different iterations
check_iter = c(50,100,250,500)

Rhat_result = matrix(
  NA,
  nrow = 6,
  ncol = 4
)

rownames(Rhat_result) = c(
  "alpha",
  "gamma_2",
  "gamma_3",
  "psi_1",
  "psi_2",
  "psi_3"
)

colnames(Rhat_result) = iterhalf + check_iter

for(i in 1:length(check_iter)){
  
  end_iter = iterhalf + check_iter[i]
  
  Rhat_result["alpha",i] = get_percent_Rhat(alpha_record_chain, end_iter, type = "alpha", 0.95)
  
  Rhat_gamma_temp = get_percent_Rhat(gamma_record_chain, end_iter, type = "gamma", 0.95)
  Rhat_result["gamma_2",i] = Rhat_gamma_temp[2]
  Rhat_result["gamma_3",i] = Rhat_gamma_temp[3]
  
  Rhat_psi_temp = get_percent_Rhat(psi_record_chain, end_iter, type = "psi", 0.95)
  Rhat_result["psi_1",i] = Rhat_psi_temp[1]
  Rhat_result["psi_2",i] = Rhat_psi_temp[2]
  Rhat_result["psi_3",i] = Rhat_psi_temp[3]
  
}

Rhat_result = as.data.frame(Rhat_result)


## Save the 95th percentile of R_hat
write.csv(round(Rhat_result, 3), file = "../result_data/R_hat_95.csv")





