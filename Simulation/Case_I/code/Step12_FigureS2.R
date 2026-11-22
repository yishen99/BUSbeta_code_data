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


## Set seed for BUSbeta in simulation case I
seed = 2002
set.seed(seed)


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


## Initiation for parameters
pi_0 <- matrix(NA, nrow = B, ncol = K)

Z_0 <- list()
for(b in 1:B){
  Z_0[[b]] <- sample(1:K, size = n_vec[b], replace = T)
  pi_0[b, ] <- as.data.frame(table(Z_0[[b]])) $ Freq / n_vec[b]
}

### xi[j,1]==0, gamma[1,j]==0
alpha_0 <- rnorm(J, mean = eta_alpha, sd = tau_alpha)
tau_xi_0 <- rinvgamma(1, shape = p_tau, rate = q_tau)
xi_0 <- matrix(c(rep(0,J), rnorm(J*(K-1), mean = 0, sd = tau_xi_0)), nrow = J, ncol = K)
gamma_0 <- matrix(c(rep(0,J), rnorm((B-1)*J, mean = 0, sd = tau_gamma)), nrow = B, ncol = J, byrow = T)
psi_0 <- matrix(rlnorm(B*J, meanlog = eta_psi, sdlog = tau_psi), nrow = B, ncol = J, byrow = T)

pi_t = pi_0
Z_t = Z_0
alpha_t = alpha_0
xi_t = xi_0
gamma_t = gamma_0
psi_t = psi_0
tau_xi_t = tau_xi_0

pi_t_record = list()
alpha_t_record = list()
xi_t_record = list()
gamma_t_record = list()
psi_t_record = list()
Z_t_record = list()


## Gibbs sampling
t1 = Sys.time()
cat("  running the Gibbs sampler ...\n")
iter = 1000
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
  
  pi_t_record[[t]] = pi_t
  alpha_t_record[[t]] = alpha_t
  xi_t_record[[t]] = xi_t
  gamma_t_record[[t]] = gamma_t
  psi_t_record[[t]] = psi_t
  Z_t_record[[t]] = Z_t
  
  if(t %% (iter / 10) == 0) cat("######## GIBBS SAMPLING IS ", t / iter * 100, "% COMPLETE ########", "\n")
} 
t2 = Sys.time()

cat(paste0("  The Gibbs sampler takes: ", round(difftime( t2, t1, units = "mins"), 3), " mins", "\n"))




## Choose parameters randomly
set.seed(111)
### pi_1k, pi_2k, pi_3k
pi_trace_ind <- sample(1:K, 1, replace = FALSE)

set.seed(222)
### alpha_j
alpha_trace_ind <- sample(1:J, 3, replace = FALSE)

set.seed(333)
### xi_j3
xi_trace_ind <- sample(1:J, 3, replace = FALSE)

set.seed(444)
### gamma_2j
gamma_trace_ind <- sample(1:J, 3, replace = FALSE)

set.seed(555)
### psi_1j, psi_2j, psi_3j
psi_trace_ind <- sample(1:J, 1, replace = FALSE)

print(list(pi=pi_trace_ind, alpha=alpha_trace_ind, xi=xi_trace_ind, gamma=gamma_trace_ind, psi=psi_trace_ind))



## Create data frames and store trajectory values of parameters
pi_trace_1 <- c(NA)
pi_trace_2 <- c(NA)
pi_trace_3 <- c(NA)
alpha_trace_1 <- c(NA)
alpha_trace_2 <- c(NA)
alpha_trace_3 <- c(NA)
xi_trace_1 <- c(NA)
xi_trace_2 <- c(NA)
xi_trace_3 <- c(NA)
gamma_trace_1 <- c(NA)
gamma_trace_2 <- c(NA)
gamma_trace_3 <- c(NA)
psi_trace_1 <- c(NA)
psi_trace_2 <- c(NA)
psi_trace_3 <- c(NA)

for(t in 1:iter){
  pi_trace_1[t] <- pi_t_record[[t]][1, pi_trace_ind]
  pi_trace_2[t] <- pi_t_record[[t]][2, pi_trace_ind]
  pi_trace_3[t] <- pi_t_record[[t]][3, pi_trace_ind]
  alpha_trace_1[t] <- alpha_t_record[[t]][alpha_trace_ind[1]]
  alpha_trace_2[t] <- alpha_t_record[[t]][alpha_trace_ind[2]]
  alpha_trace_3[t] <- alpha_t_record[[t]][alpha_trace_ind[3]]
  xi_trace_1[t] <- xi_t_record[[t]][xi_trace_ind[1], 3]
  xi_trace_2[t] <- xi_t_record[[t]][xi_trace_ind[2], 3]
  xi_trace_3[t] <- xi_t_record[[t]][xi_trace_ind[3], 3]
  gamma_trace_1[t] <- gamma_t_record[[t]][2, gamma_trace_ind[1]]
  gamma_trace_2[t] <- gamma_t_record[[t]][2, gamma_trace_ind[2]]
  gamma_trace_3[t] <- gamma_t_record[[t]][2, gamma_trace_ind[3]]
  psi_trace_1[t] <- psi_t_record[[t]][1, psi_trace_ind]
  psi_trace_2[t] <- psi_t_record[[t]][2, psi_trace_ind]
  psi_trace_3[t] <- psi_t_record[[t]][3, psi_trace_ind]
}


df_trace <- data.frame(
  pi_trace_1, pi_trace_2, pi_trace_3,
  alpha_trace_1, alpha_trace_2, alpha_trace_3, 
  xi_trace_1, xi_trace_2, xi_trace_3,
  gamma_trace_1, gamma_trace_2, gamma_trace_3,
  psi_trace_1, psi_trace_2, psi_trace_3
)


## Draw trace plots for each chosen parameter
### pi
pi_1_trace_plot <- ggplot(df_trace[(iterhalf+1):iter, ], aes(x = (iterhalf+1):iter, y = pi_trace_1)) +
  geom_line() + 
  theme_gray() +
  labs(x = "MCMC iteration", y = NULL) + 
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS2(a).png", plot = pi_1_trace_plot, width=9, height=6)

pi_2_trace_plot <- ggplot(df_trace[(iterhalf+1):iter, ], aes(x = (iterhalf+1):iter, y = pi_trace_2)) +
  geom_line() + 
  theme_gray() +
  labs(x = "MCMC iteration", y = NULL) + 
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS2(b).png", plot = pi_2_trace_plot, width=9, height=6)

pi_3_trace_plot <- ggplot(df_trace[(iterhalf+1):iter, ], aes(x = (iterhalf+1):iter, y = pi_trace_3)) +
  geom_line() + 
  theme_gray() +
  labs(x = "MCMC iteration", y = NULL) + 
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS2(c).png", plot = pi_3_trace_plot, width=9, height=6)


### alpha
alpha_1_trace_plot <- ggplot(df_trace[(iterhalf+1):iter, ], aes(x = (iterhalf+1):iter, y = alpha_trace_1)) +
  geom_line() + 
  theme_gray() +
  labs(x = "MCMC iteration", y = NULL) + 
  coord_cartesian(ylim = c(-0.25, 0.25)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS2(d).png", plot = alpha_1_trace_plot, width=9, height=6)

alpha_2_trace_plot <- ggplot(df_trace[(iterhalf+1):iter, ], aes(x = (iterhalf+1):iter, y = alpha_trace_2)) +
  geom_line() + 
  theme_gray() +
  labs(x = "MCMC iteration", y = NULL) + 
  coord_cartesian(ylim = c(-0.25, 0.25)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS2(e).png", plot = alpha_2_trace_plot, width=9, height=6)

alpha_3_trace_plot <- ggplot(df_trace[(iterhalf+1):iter, ], aes(x = (iterhalf+1):iter, y = alpha_trace_3)) +
  geom_line() + 
  theme_gray() +
  labs(x = "MCMC iteration", y = NULL) + 
  coord_cartesian(ylim = c(-0.25, 0.25)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS2(f).png", plot = alpha_3_trace_plot, width=9, height=6)


### xi
xi_1_trace_plot <- ggplot(df_trace[(iterhalf+1):iter, ], aes(x = (iterhalf+1):iter, y = xi_trace_1)) +
  geom_line() + 
  theme_gray() +
  labs(x = "MCMC iteration", y = NULL) + 
  coord_cartesian(ylim = c(1, 1.9)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS2(g).png", plot = xi_1_trace_plot, width=9, height=6)

xi_2_trace_plot <- ggplot(df_trace[(iterhalf+1):iter, ], aes(x = (iterhalf+1):iter, y = xi_trace_2)) +
  geom_line() + 
  theme_gray() +
  labs(x = "MCMC iteration", y = NULL) + 
  coord_cartesian(ylim = c(1, 1.9)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS2(h).png", plot = xi_2_trace_plot, width=9, height=6)

xi_3_trace_plot <- ggplot(df_trace[(iterhalf+1):iter, ], aes(x = (iterhalf+1):iter, y = xi_trace_3)) +
  geom_line() + 
  theme_gray() +
  labs(x = "MCMC iteration", y = NULL) + 
  coord_cartesian(ylim = c(1, 1.9)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS2(i).png", plot = xi_3_trace_plot, width=9, height=6)

### gamma
gamma_1_trace_plot <- ggplot(df_trace[(iterhalf+1):iter, ], aes(x = (iterhalf+1):iter, y = gamma_trace_1)) +
  geom_line() + 
  theme_gray() +
  labs(x = "MCMC iteration", y = NULL) + 
  coord_cartesian(ylim = c(-0.05, 0.2)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS2(j).png", plot = gamma_1_trace_plot, width=9, height=6)

gamma_2_trace_plot <- ggplot(df_trace[(iterhalf+1):iter, ], aes(x = (iterhalf+1):iter, y = gamma_trace_2)) +
  geom_line() + 
  theme_gray() +
  labs(x = "MCMC iteration", y = NULL) + 
  coord_cartesian(ylim = c(-0.05, 0.2)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS2(k).png", plot = gamma_2_trace_plot, width=9, height=6)

gamma_3_trace_plot <- ggplot(df_trace[(iterhalf+1):iter, ], aes(x = (iterhalf+1):iter, y = gamma_trace_3)) +
  geom_line() + 
  theme_gray() +
  labs(x = "MCMC iteration", y = NULL) + 
  coord_cartesian(ylim = c(-0.05, 0.2)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS2(l).png", plot = gamma_3_trace_plot, width=9, height=6)


### psi
psi_1_trace_plot <- ggplot(df_trace[(iterhalf+1):iter, ], aes(x = (iterhalf+1):iter, y = psi_trace_1)) +
  geom_line() + 
  theme_gray() +
  labs(x = "MCMC iteration", y = NULL) + 
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS2(m).png", plot = psi_1_trace_plot, width=9, height=6)

psi_2_trace_plot <- ggplot(df_trace[(iterhalf+1):iter, ], aes(x = (iterhalf+1):iter, y = psi_trace_2)) +
  geom_line() + 
  theme_gray() +
  labs(x = "MCMC iteration", y = NULL) + 
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS2(n).png", plot = psi_2_trace_plot, width=9, height=6)

psi_3_trace_plot <- ggplot(df_trace[(iterhalf+1):iter, ], aes(x = (iterhalf+1):iter, y = psi_trace_3)) +
  geom_line() + 
  theme_gray() +
  labs(x = "MCMC iteration", y = NULL) + 
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS2(o).png", plot = psi_3_trace_plot, width=9, height=6)










