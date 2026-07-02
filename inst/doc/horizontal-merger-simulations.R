## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(mergersim)

## -----------------------------------------------------------------------------
library(BB)
library(rootSolve)
library(numDeriv)

## -----------------------------------------------------------------------------
alpha_true  <- -0.9
delta_true <- c(.81,.93,.82)
c_j <- c(.05,.31,.30)

## -----------------------------------------------------------------------------
own_pre = diag(3)
own_post <- own_pre
own_post[1,2] <- 1
own_post[2,1] <- 1

## -----------------------------------------------------------------------------
x0 <- c_j*1.1
out1 <- multiroot(f = bertrand_foc, start = x0, 
                  own = own_pre, alpha= alpha_true, 
                  delta = delta_true, cost = c_j)

p1 <- out1$root
p1
share1 <- (exp(delta_true + alpha_true*p1))/(1+sum(exp(delta_true + alpha_true*p1)))
share1

## -----------------------------------------------------------------------------
x0 <- c_j*1.1
out1_post <- multiroot(f = bertrand_foc ,start = x0, 
                  own = own_post, alpha= alpha_true, 
                  delta = delta_true, cost = c_j)

p1_post <- out1_post$root
p1_post
share1_post <- (exp(delta_true + alpha_true*p1_post))/(1+sum(exp(delta_true + alpha_true*p1_post)))
share1_post

# margins. Note that multi-product firm sets same level margins.
(p1_post - c_j)
(p1_post - c_j)/p1_post

# price effect
(p1_post - p1)/p1_post


## ----results = "hide", warning=FALSE------------------------------------------
x00 <- c(-1)
wt_matrix <- diag(c(1,1,1,1000,1000,1000))

out1 <- optim(f = bertrand_calibrate, par = x00, 
                   own = own_pre, price = p1, 
                   shares = share1, cost  = c_j,
                   weight = wt_matrix)




## -----------------------------------------------------------------------------

# note that this optimization recovers the true demand parameters
out1$par
alpha_true

delta_cal <- log(share1) - log(1-sum(share1)) - out1$par*p1
delta_cal
delta_true

## ----eval = TRUE--------------------------------------------------------------
c_obs <- c_j
c_obs[3] <- NA
c_obs

## ----results = "hide", warning=FALSE------------------------------------------
x00 <- c(-1)
wt_matrix <- diag(c(1,1,1,1000,1000,1000))

out1b <- optim(f = bertrand_calibrate, par = x00, 
                   own = own_pre, price = p1, 
                   shares = share1, cost  = c_obs,
                   weight = wt_matrix)




## -----------------------------------------------------------------------------

# note that this optimization recovers the true demand parameters
out1b$par
alpha_true

delta_cal <- log(share1) - log(1-sum(share1)) - out1b$par*p1
delta_cal
delta_true

## ----eval = TRUE--------------------------------------------------------------
share2 <- (exp(delta_true + alpha_true*c_j))/(1+sum(exp(delta_true + alpha_true*c_j)))

p2 <- c_j + log(1 - own_pre%*%share2)/(alpha_true*own_pre%*%share2)
p2_post <- c_j + log(1 - own_post%*%share2)/(alpha_true*own_post%*%share2)

(p2_post - p2)/p2

## -----------------------------------------------------------------------------
c_obs <- c_j
c_obs[3] <- NA
c_obs

## -----------------------------------------------------------------------------
wt_matrix <- diag(c(1,1))

result4 <- BBoptim(f = ssa_calibrate, par = c(-.2),
                    lower = c(-Inf), upper = c(-0.0001),
                    own=own_pre, price = p2, share = share2,
                    cost = c_obs, weight = wt_matrix)

alpha4 <- result4$par
alpha4       # recover true value

## -----------------------------------------------------------------------------
cost4 <- p2 - log(1 - own_pre%*%share2) / (alpha4*own_pre%*%share2)
cost4

## -----------------------------------------------------------------------------
c_j

## -----------------------------------------------------------------------------
p4 <- cost4 + log(1 - own_pre%*%share2)/(alpha4*own_pre%*%share2)
p4_post <- cost4 + log(1 - own_post%*%share2)/(alpha4*own_post%*%share2)

(p4_post - p4)/p4

## -----------------------------------------------------------------------------
1/(p4[1] - c_obs[1]) * log(1 - share2[1]) / share2[1]

## -----------------------------------------------------------------------------
x0 <- c_j*1.5

out3 <- multiroot(f = bargain_foc, start = x0, own = own_pre, 
                  alpha = alpha_true, delta = delta_true, cost = c_j,
                  lambda = 0.5)
p3 <- out3$root
share3 <- (exp(delta_true + alpha_true*p3))/(1+sum(exp(delta_true + alpha_true*p3)))

print(p3)
print(share3)


## ----eval = TRUE--------------------------------------------------------------
J <- length(c_j)
alpha_start <- -1.2
delta_start <- rep(1,J)
x00 <- c(alpha_start,delta_start)
wt_matrix <- diag(J*2)


bargain_calibrate(param = x00, 
                       own = own_pre, price = p3, 
                       shares = share3, cost  = c_j,
                       weight = wt_matrix, lambda = 0.5)


out3 <- BBoptim(f = bargain_calibrate, par = x00, 
                own = own_pre, price = p3, 
                shares = share3, cost  = c_j,
                weight = wt_matrix, lambda = 0.5)

# check if we recovered correct demand parameters
# finding good initial values is important.
alpha3 <- out3$par[1]
delta3 <- out3$par[2:4]
alpha_true
alpha3
delta_true
delta3

## -----------------------------------------------------------------------------
x0 <- c_j*1.5

out4 <- multiroot(f = ssbargain_foc, start = x0, own = own_pre, 
                       alpha = alpha_true, delta = delta_true, 
                  cost = c_j, lambda = 0.5)
p4 <- out4$root
share4 <- (exp(delta_true + alpha_true*c_j))/(1+sum(exp(delta_true + alpha_true*c_j)))

print(p4)
print(share4)


## -----------------------------------------------------------------------------
alpha_start <- -1.2
delta_start <- rep(1,J)
x00 <- c(alpha_start,delta_start)
wt_matrix <- diag(J*2)

out4 <- BBoptim(f = ssbargain_calibrate, par = x00, 
                own = own_pre, price = p4, 
                shares = share4, cost  = c_j,
                weight = wt_matrix, lambda = 0.5)

# check if we recovered correct demand parameters
# finding good initial values is important.
alpha4 <- out4$par[1]
delta4 <- out4$par[2:4]
alpha_true
alpha4
delta_true
delta4

## -----------------------------------------------------------------------------
K1 <- 2

B1 <- 1 * matrix( c(1,0,
               1,0,
               0,1),
             ncol = K1, nrow = J, byrow = TRUE)
a1 <- B1    # rows of a should sum to 1 to facilitate interpretation.
mu1 <- rep(1.0,K1) # nesting parameters all 1 simplifies to logit
mu2 <- c(0.8,0.8)

## ----results = "hide"---------------------------------------------------------

x0 <- p1*1.1

out1 <- BBoptim(fn = bertrand_foc, par = x0, 
                own = own_pre, alpha = alpha_true, 
                delta = delta_true, cost = c_j, sumFOC = TRUE)

p_R1 <- out1$par

## ----eval = FALSE-------------------------------------------------------------
#  # Equilibrium prices are the same as we had from logit model before
#  p_R1
#  p1
#  

## ----results = "hide"---------------------------------------------------------

x0 <- p1*1.1

out2 <- BBoptim(fn = bertrand_foc, par = x0, 
                own = own_pre, alpha = alpha_true, 
                delta = delta_true, cost = c_j,
                nest_allocation=a1, mu=mu2,
                sumFOC = TRUE)

p_R2 <- out2$par

## -----------------------------------------------------------------------------
# Equilibrium prices are different than logit result
p_R2
  

## -----------------------------------------------------------------------------
diversions <- diversion_calc(price=p_R2,alpha=alpha_true,delta=delta_true,
                             nest_allocation=a1,mu=mu2)

diversions
rowSums(diversions) # each row should sum to <1 because of outside option

## -----------------------------------------------------------------------------
diversions2 <- diversion_calc(price=p_R2,alpha=alpha_true, delta=delta_true,
                              nest_allocation=a1, mu=mu2, marginal = TRUE)

diversions2
rowSums(diversions2) # each row should sum to <1 because of outside option

