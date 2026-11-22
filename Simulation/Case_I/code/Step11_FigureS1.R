##############################################################
#################     Simulation Case I    ###################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

library(invgamma)
library(truncnorm)
library(aricode)
library(ggplot2)


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

update_tau_xi_SA <- function(xi_t, p_tau, q_tau){
  return(sqrt(rinvgamma(1, shape = p_tau + (K-1)*J/2, rate = q_tau + sum(xi_t[ ,-1]^2)/2)))
}

update_pi_SA <- function(Z_t, c){
  return(t(sapply(1:B, function(b){Dirichlet(sapply(1:K, function(k){sum(Z_t[[b]] == k)}) + c)})))
}

update_Z_SA <- function(pi_t, alpha_t, xi_t, gamma_t, psi_t){
  temp = lapply(1:B, function(b){
    sapply(1:n_vec[b], function(i){
      log_prob = sapply(1:K, function(k){log(pi_t[b, k]) + sum(logBe(Y[[b]][i, ], 1/(1 + exp(-alpha_t-xi_t[, k]-gamma_t[b, ])), psi_t[b, ]))})
      Prob = exp(log_prob-max(log_prob)) / sum(exp(log_prob-max(log_prob)))
      sample(1:K, size = 1, prob = Prob)
    })
  })
  return(temp)
}

update_alpha_SA <- function(Z_t, alpha_t, xi_t, gamma_t, psi_t, eta_alpha, tau_alpha, phi_alpha){
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

update_xi_SA <- function(tau_xi_t, Z_t, alpha_t, xi_t, gamma_t, psi_t, phi_xi){
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


update_gamma_SA <- function(Z_t, alpha_t, xi_t, gamma_t, psi_t, tau_gamma, phi_gamma){
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


update_psi_SA <- function(Z_t, alpha_t, xi_t, gamma_t, psi_t, eta_psi, tau_psi, phi_psi){
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





## Function of sensitivity analysis for each hyper-parameter
hyperparameter_sensitivity_analysis <- function(hyperparam_name, value) {
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
  
  
  if (hyperparam_name == "c") {
    c <- value
  }
  else if (hyperparam_name == "eta_alpha") {
    eta_alpha <- value
  }
  else if (hyperparam_name == "tau_alpha") {
    tau_alpha <- value
  }
  else if (hyperparam_name == "p_tau") {
    p_tau <- value
  }
  else if (hyperparam_name == "q_tau") {
    q_tau <- value
  }
  else if (hyperparam_name == "tau_gamma") {
    tau_gamma <- value
  }
  else if (hyperparam_name == "eta_psi") {
    eta_psi <- value
  }
  else if (hyperparam_name == "tau_psi") {
    tau_psi <- value
  }
  else if (hyperparam_name == "phi_alpha") {
    phi_alpha <- rep(value, J)
  }
  else if (hyperparam_name == "phi_xi") {
    phi_xi <- matrix(value, nrow = J, ncol = K)
  }
  else if (hyperparam_name == "phi_gamma") {
    phi_gamma <- matrix(value, nrow = B, ncol = J)
  }
  else if (hyperparam_name == "phi_psi") {
    phi_psi <- matrix(value, nrow = B, ncol = J)
  }
  else
    return("ERROR!")
  
  
  ## Set seed for BUSbeta in simulation case I
  seed = 2002
  set.seed(seed)
  
  ## Initiation for parameters
  pi_0 <- matrix(NA, nrow = B, ncol = K)
  
  Z_0 <- list()
  for (b in 1:B) {
    Z_0[[b]] <- sample(1:K, size = n_vec[b], replace = T)
    pi_0[b, ] <- as.data.frame(table(Z_0[[b]]))$Freq / n_vec[b]
  }
  
  ### xi[j,1]==0, gamma[1,j]==0
  alpha_0 <- rnorm(J, mean = eta_alpha, sd = tau_alpha)
  tau_xi_0 <- rinvgamma(1, shape = p_tau, rate = q_tau)
  xi_0 <- matrix(c(rep(0, J), rnorm(
    J * (K - 1), mean = 0, sd = tau_xi_0
  )), nrow = J, ncol = K)
  gamma_0 <- matrix(c(rep(0, J), rnorm((B - 1) * J, mean = 0, sd = tau_gamma
  )),
  nrow = B,
  ncol = J,
  byrow = T)
  psi_0 <- matrix(
    rlnorm(B * J, meanlog = eta_psi, sdlog = tau_psi),
    nrow = B,
    ncol = J,
    byrow = T
  )
  
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
  cat("  running the Gibbs sampler for",
      hyperparam_name,
      "=",
      value,
      "...\n")
  iter = 1000
  iterhalf = iter / 2
  for (t in 1:iter) {
    ### tau_xi
    tau_xi_t = update_tau_xi_SA(xi_t, p_tau, q_tau)
    ### pi
    pi_t = update_pi_SA(Z_t, c)
    ### Z
    Z_t = update_Z_SA(pi_t, alpha_t, xi_t, gamma_t, psi_t)
    ### alpha
    alpha_t = update_alpha_SA(Z_t, alpha_t, xi_t, gamma_t, psi_t, eta_alpha, tau_alpha, phi_alpha)
    ### xi
    xi_t = update_xi_SA(tau_xi_t, Z_t, alpha_t, xi_t, gamma_t, psi_t, phi_xi)
    ### gamma
    gamma_t = update_gamma_SA(Z_t, alpha_t, xi_t, gamma_t, psi_t, tau_gamma, phi_gamma)
    ### psi
    psi_t = update_psi_SA(Z_t, alpha_t, xi_t, gamma_t, psi_t, eta_psi, tau_psi, phi_psi)
    
    alpha_t_record[[t]] = alpha_t
    xi_t_record[[t]] = xi_t
    gamma_t_record[[t]] = gamma_t
    psi_t_record[[t]] = psi_t
    Z_t_record[[t]] = Z_t
    
    if (t %% (iter / 10) == 0)
      cat("######## GIBBS SAMPLING IS ",
          t / iter * 100,
          "% COMPLETE ########",
          "\n")
  }
  
  
  ## Estimated parameters
  alpha_hat <- 0
  xi_hat <- 0
  gamma_hat <- 0
  psi_hat <- 0
  Z_BUSbeta <- list()
  mu_hat <- list()
  sigma_hat <- list()
  sigma_weighted_hat <- list()
  
  
  for (t in (iter - iterhalf + 1):iter) {
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
  for (b in 1:B) {
    Z_BUSbeta[[b]] = rep(NA, times = n_vec[b])
    for (i in 1:n_vec[b]) {
      Z_BUSbeta[[b]][i] = getMode(sapply((iter - iterhalf + 1):iter, function(t) {
        Z_t_record[[t]][[b]][i]
      }))
    }
  }
  
  ## Compute ARI
  ARI_BUSbeta <- data.frame(
    Batch1 = round(ARI(Z_real[[1]], Z_BUSbeta[[1]]), 3),
    Batch2 = round(ARI(Z_real[[2]], Z_BUSbeta[[2]]), 3),
    Batch3 = round(ARI(Z_real[[3]], Z_BUSbeta[[3]]), 3),
    Overall = round(ARI(unlist(Z_real), unlist(Z_BUSbeta)), 3)
  )
  
  return(ARI_BUSbeta)
  
}


hyperparam_names <- c(
  "c",
  "eta_alpha",
  "tau_alpha",
  "p_tau",
  "q_tau",
  "tau_gamma",
  "eta_psi",
  "tau_psi",
  "phi_alpha",
  "phi_xi",
  "phi_gamma",
  "phi_psi"
)

## Set the range for varying hyper-parameters
hyperparam_value_list <- list(
  c = seq(0.5, 1.5, by = 0.05),
  eta_alpha = seq(-0.2, 0.2, by = 0.02),
  tau_alpha = seq(0.05, 0.15, by = 0.005),
  p_tau = seq(2.5, 12.5, by = 0.5),
  q_tau = seq(0.5, 1.5, by = 0.05),
  tau_gamma = seq(0.005, 0.105, by = 0.005),
  eta_psi = seq(0.25, 1.25, by = 0.05),
  tau_psi = seq(0.4, 1.4, by = 0.05),
  phi_alpha = seq(0.01, 0.11, by = 0.005),
  phi_xi = seq(0.05, 0.15, by = 0.005),
  phi_gamma = seq(0.01, 0.05, by = 0.002),
  phi_psi = seq(0.5, 3, by = 0.125)
)




## Run sensitivity analysis and create 12 data frames
for (hyperparam_name in names(hyperparam_value_list)){
  
  hyperparam_value <- hyperparam_value_list[[hyperparam_name]]
  value_ARI_list <- list()
  for(ind in 1:21){
    new_row = hyperparameter_sensitivity_analysis(hyperparam_name, hyperparam_value[ind])
    value_ARI_list[[ind]] = new_row
  }
  
  df_name <- paste0("df_", hyperparam_name, "_ARI")
  assign(
    df_name,
    cbind(hyperparameter = hyperparam_name, value = hyperparam_value, do.call(rbind, value_ARI_list))
  )
}


# df_c_ARI = data.frame(hyperparameter = "c", 
#                      value = seq(0.5, 1.5, by = 0.05),
#                      Batch1 = 1, Batch2 = 1, Batch3 = 1, Overall = 1)

# df_eta_alpha_ARI = data.frame(hyperparameter = "eta_alpha", 
#                               value = seq(-0.2, 0.2, by = 0.02),
#                               Batch1 = 1, Batch2 = 1, Batch3 = 1, Overall = 1)

# df_tau_alpha_ARI = data.frame(hyperparameter = "tau_alpha", 
#                               value = seq(0.05, 0.15, by = 0.005),
#                               Batch1 = 1, Batch2 = 1, Batch3 = 1, Overall = 1)

# df_p_tau_ARI = data.frame(hyperparameter = "p_tau", 
#                           value = seq(2.5, 12.5, by = 0.5),
#                           Batch1 = c(0.331,0.854,0.854,0.854,0.854,rep(1,16)), 
#                           Batch2 = c(0.677,0.895,0.895,0.895,0.895,rep(1,16)), 
#                           Batch3 = c(0.538,0.839,0.839,0.839,0.839,rep(1,16)), 
#                           Overall = c(0.455,0.761,0.761,0.761,0.761,rep(1,16)))

# df_q_tau_ARI = data.frame(hyperparameter = "q_tau", 
#                           value = seq(0.5, 1.5, by = 0.05),
#                           Batch1 = c(rep(1,20),0.854), 
#                           Batch2 = c(rep(1,20),0.895), 
#                           Batch3 = c(rep(1,20),0.839), 
#                           Overall = c(rep(1,20),0.761))

# df_tau_gamma_ARI = data.frame(hyperparameter = "tau_gamma", 
#                           value = seq(0.005, 0.105, by = 0.005),
#                           Batch1 = 1, 
#                           Batch2 = c(1,1,0.981,0.993,rep(1,17)), 
#                           Batch3 = c(0.839,0.825,1,1,rep(1,17)), 
#                           Overall = c(0.733,0.735,0.989,0.996,rep(1,17)))

# df_eta_psi_ARI = data.frame(hyperparameter = "eta_psi", 
#                               value = seq(0.25, 1.25, by = 0.05),
#                               Batch1 = c(0.854,0.854,0.854,0.854,1,1,0.844,rep(1,14)), 
#                               Batch2 = c(0.895,0.895,0.895,0.895,1,1,0.895,rep(1,14)), 
#                               Batch3 = c(0.839,0.839,0.839,0.839,1,1,0.839,rep(1,14)), 
#                               Overall = c(0.761,0.761,0.761,0.761,1,1,0.758,rep(1,14)))

# df_tau_psi_ARI = data.frame(hyperparameter = "tau_psi", 
#                             value = seq(0.4, 1.4, by = 0.05),
#                             Batch1 = c(0.816,0.816,0.816,rep(1,18)), 
#                             Batch2 = c(0.895,0.895,0.895,rep(1,18)), 
#                             Batch3 = c(0.839,0.839,0.839,rep(1,18)), 
#                             Overall = c(0.781,0.781,0.781,rep(1,18)))

# df_phi_alpha_ARI = data.frame(hyperparameter = "phi_alpha", 
#                       value = seq(0.01, 0.11, by = 0.005),
#                       Batch1 = 1, Batch2 = 1, Batch3 = 1, Overall = 1)

# df_phi_xi_ARI = data.frame(hyperparameter = "phi_xi", 
#                               value = seq(0.05, 0.15, by = 0.005),
#                               Batch1 = 1, Batch2 = 1, Batch3 = 1, Overall = 1)

# df_phi_gamma_ARI = data.frame(hyperparameter = "phi_gamma", 
#                            value = seq(0.01, 0.05, by = 0.002),
#                            Batch1 = 1, Batch2 = 1, Batch3 = 1, Overall = 1)

# df_phi_psi_ARI = data.frame(hyperparameter = "phi_psi", 
#                              value = seq(0.5, 3, by = 0.125),
#                              Batch1 = 1, Batch2 = 1, Batch3 = 1, Overall = 1)


# df_SA_ARI = rbind(df_c_ARI, df_eta_alpha_ARI, df_tau_alpha_ARI, df_p_tau_ARI,
#                   df_q_tau_ARI, df_tau_gamma_ARI, df_eta_psi_ARI, df_tau_psi_ARI,
#                   df_phi_alpha_ARI, df_phi_xi_ARI, df_phi_gamma_ARI, df_phi_psi_ARI)



## Draw ARI plots for each hyper-parameter
c_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "c"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "c", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(a).png", plot = c_SY_plot, width=9, height=6)


eta_alpha_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "eta_alpha"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "eta_alpha", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(b).png", plot = eta_alpha_SY_plot, width=9, height=6)


tau_alpha_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "tau_alpha"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "tau_alpha", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(c).png", plot = tau_alpha_SY_plot, width=9, height=6)


p_tau_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "p_tau"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "p_tau", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(d).png", plot = p_tau_SY_plot, width=9, height=6)
  

q_tau_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "q_tau"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "q_tau", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(e).png", plot = q_tau_SY_plot, width=9, height=6)


tau_gamma_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "tau_gamma"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "tau_gamma", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(f).png", plot = tau_gamma_SY_plot, width=9, height=6)


eta_psi_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "eta_psi"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "eta_psi", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(g).png", plot = eta_psi_SY_plot, width=9, height=6)


tau_psi_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "tau_psi"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "tau_psi", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(h).png", plot = tau_psi_SY_plot, width=9, height=6)


phi_alpha_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "phi_alpha"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "phi_alpha", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(i).png", plot = phi_alpha_SY_plot, width=9, height=6)


phi_xi_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "phi_xi"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "phi_xi", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(j).png", plot = phi_xi_SY_plot, width=9, height=6)


phi_gamma_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "phi_gamma"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "phi_gamma", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(k).png", plot = phi_gamma_SY_plot, width=9, height=6)


phi_psi_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "phi_psi"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "phi_psi", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(l).png", plot = phi_psi_SY_plot, width=9, height=6)







