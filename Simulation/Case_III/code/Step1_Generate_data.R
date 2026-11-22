##############################################################
#################     Simulation Case III    ###################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################


rBe <- function(mu, psi){
  s1 = mu * psi
  s2 = psi - mu * psi
  beta_rv = rbeta(1, shape1 = s1, shape2 = s2)
  return(beta_rv)
}


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


B <- 3
K <- 5
J <- 1000

## Y: dim B * (n_b (b=1,...,B) * J) matrix list
## mu: dim B * (J * K) matrix list
## n_vec: dim B vector
## pi: dim B * J matrix
## Z: dim B * n_b (b=1,...,B) vector list
## alpha: dim J vector
## xi: dim J * K matrix with first column being zeros
## gamma: dim B * J matrix with first row being ones
## psi: dim B * J matrix

n_vec <- c(200, 220, 240)

pi_real <- matrix(c(0.3,0.2,0.0,0.2,0.3,
                    0.7,0.3,0.0,0.0,0.0,
                    0.0,0.0,0.4,0.4,0.2), nrow = B, ncol = K, byrow = T) 


Z_real <- list()
for(b in 1:B){
  Z_real[[b]] <- c(rep(1, pi_real[b,1]*n_vec[b]), rep(2, pi_real[b,2]*n_vec[b]), 
                   rep(3, pi_real[b,3]*n_vec[b]), rep(4, pi_real[b,4]*n_vec[b]), 
                   rep(5, pi_real[b,5]*n_vec[b]))
}

alpha_real <- rep(0, times=J)
xi_real <- matrix(c(rep(0,J), rep(1.5,J), rep(-1.49,J), rep(0.5,J), rep(-0.51,J)), nrow = J, ncol = K)
gamma_real <- matrix(c(rep(0,J), rep(0.2,J), rep(-0.2,J)), nrow = B, ncol = J, byrow = T)
psi_real <- matrix(c(rep(2,J), rep(6,J), rep(4,J)), nrow = B, ncol = J, byrow = T)


## Set seed for data generation in simulation case III
set.seed(1)


## Data generation
mu <- list()
Y <- list()
for(b in 1:B){
  mu[[b]] <- 1 / (1 + exp(-matrix(rep(alpha_real, K), nrow =J, ncol = K) - xi_real - matrix(rep(gamma_real[b, ], K), nrow = J, ncol = K)))
  Y[[b]] <- matrix(NA, nrow = n_vec[b], ncol = J)
  for(i in 1:n_vec[b]){
    Y[[b]][i, ] = sapply(1:J, function(j){rBe(mu[[b]][j, Z_real[[b]][i]], psi_real[b, j])})
  }
}


## True subgroup mean
truemean_by_subgroup <- list()
for(k in 1:K) truemean_by_subgroup[[k]] = matrix(rep(1/(1+exp(-alpha_real-xi_real[, k])), times = length(which(unlist(Z_real)==k))), ncol=J, byrow = TRUE)


## Save Raw simulated data and subgroup label
save(Y, file = "../input_data/Y.RData")
save(Z_real, file = "../input_data/Z_real.RData")

## Save true subgroup mean
save(truemean_by_subgroup, file = "../input_data/truemean_by_subgroup.RData")


