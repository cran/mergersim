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
delta_true <- c(0.9, 0.9, 1.5, rep(3.0, 2)) + 3.5
c_j <- c(rep(2.3,2), 2.7, rep(3.7,2))
J <- 5

own_pre = diag(J)
own_post <- own_pre
own_post[1,2] <- 1
own_post[2,1] <- 1
own_post[1,3] <- 1
own_post[3,1] <- 1
own_post[2,3] <- 1
own_post[3,2] <- 1

K1 <- 2

B1 <- 1 * matrix( c(1,0,
                    1,0,
                    1,0,
                    0,1,
                    0,1),
                  ncol = K1, nrow = J, byrow = TRUE)
a1 <- B1    # rows of a should sum to 1 to facilitate interpretation.
mu1 <- rep(1.0,K1) # nesting parameters all 1 simplifies to logit
mu2 <- c(0.7,0.7)


x0 <- c_j*1.1

out1 <- BBoptim(fn = bertrand_foc, par = x0,
                own = own_pre, alpha = alpha_true,
                delta = delta_true, cost = c_j,
                nest_allocation=a1,mu=mu2,
                sumFOC = TRUE)

p_R1 <- out1$par
p_R1

shares1 <- share_calc(price=p_R1,delta=delta_true,alpha=alpha_true,
                      nest_allocation=a1,mu=mu2)
shares1
sum(shares1)

true_div <- diversion_calc(p=p_R1,alpha=alpha_true,delta=delta_true,
                           nest_allocation=a1,mu=mu2,
                         marginal = TRUE)


## ----results = "hide", warning=FALSE------------------------------------------


wt_vector <- c(1,1,1,1,1)

bertrand_calibrate_gnl(param = c(alpha_true,mu2),
                            own = own_pre, price = p_R1+.1,
                            shares = shares1, cost  = c_j,
                            weight = wt_vector, nest_allocation=a1,
                            div_matrix = true_div)


x00 <- c(alpha_true,mu2)-.1

out2 <- optim(f = bertrand_calibrate_gnl, par = x00,
               own = own_pre, price = p_R1,
               shares = shares1, cost  = c_j,
               weight = wt_vector, nest_allocation=a1, 
               div_matrix = true_div)

## -----------------------------------------------------------------------------
out2$par[1]
alpha_true
out2$par[2:3]
mu2



## -----------------------------------------------------------------------------

delta_start <- log(shares1) - log(1-sum(shares1)) - out2$par[1]*p_R1

out2_d <- multiroot(f = match_share, start = delta_start,
                   price=p_R1,alpha=out2$par[1],nest_allocation=a1,
                   mu=out2$par[2:3],
                   shares_obs = shares1)

out2_d$root
delta_true

delta_cal <- out2_d$root
alpha_cal <- out2$par[1]
mu_cal <- out2$par[2:3]


## -----------------------------------------------------------------------------

# Evaluate matching to observables

out_cal <- BBoptim(f = bertrand_foc, par = p_R1,
                own = own_pre, alpha= alpha_cal,
                delta = delta_cal, cost = c_j,
                nest_allocation = a1, mu =  mu_cal,
                sumFOC = TRUE)

price_cal <- out_cal$par
price_cal
p_R1

share_cal <- share_calc(price=price_cal,delta=delta_cal,alpha=alpha_cal,
                        nest_allocation=a1,mu=mu_cal)
share_cal
shares1

diversions_cal <- diversion_calc(price=price_cal,alpha=alpha_cal,
                                 delta=delta_cal,
                                 nest_allocation=a1,mu=mu_cal)

diversions_cal
true_div


## -----------------------------------------------------------------------------
alpha_true  <- -0.9
delta_true <- rep(c(0.9, 0.9, 1.5), times = 2) + 3.5
c_j <- rep(c(2.3, 2.3, 2.7), times = 2)
J <- length(delta_true)

own_pre <- diag(J)

K3 <- 5

B3 <- 1 * matrix( c(1,0,1,0,0,
                    1,0,0,1,0,
                    1,0,0,0,1,
                    0,1,1,0,0,
                    0,1,0,1,0,
                    0,1,0,0,1),
                  ncol = K3, nrow = J, byrow = TRUE)
a3 <- 0.5 * B3    # rows of a should sum to 1 to facilitate interpretation.
mu3 <- rep(0.7,K3) # nesting parameters all 1 simplifies to logit


## Generate observed objects: prices, shares, diversions.

x0 <- c_j*1.1

out3 <- BBoptim(fn = bertrand_foc, par = x0,
                own = own_pre, alpha = alpha_true,
                delta = delta_true, cost = c_j,
                nest_allocation=a3,mu=mu3,
                sumFOC = TRUE)

p_R3 <- out3$par
p_R3

shares3 <- share_calc(price=p_R3,delta=delta_true,alpha=alpha_true,
                      nest_allocation=a3,mu=mu3)
shares3
sum(shares3)

true_div <- diversion_calc(price=p_R3,alpha=alpha_true,delta=delta_true,
                           nest_allocation=a3,mu=mu3,
                           marginal = TRUE)


## ----results = "hide", warning=FALSE------------------------------------------

wt_vector <- c(1,1,1,1,1,1)

x00 <- c(alpha_true,mu3[1])-.2


out3 <- optim(f = bertrand_calibrate_gnl, par = x00,
              own = own_pre, price = p_R3,
              shares = shares3, cost  = c_j,
              weight = wt_vector, nest_allocation=a3, div_matrix = true_div,
              div_calc_marginal = TRUE,
              control = list(maxit = 1000))


## -----------------------------------------------------------------------------

out3$par[1]
alpha_true
out3$par[2:length(out3$par)]
mu3



## -----------------------------------------------------------------------------
mu_constraint_matrix <- matrix(1, nrow = K3, ncol = 1)
mu_cal <- mu_constraint_matrix %*% out3$par[2]

delta_start <- log(shares3) - log(1-sum(shares3)) - out3$par[1]*p_R3

out3b <- multiroot(f = match_share, start = delta_start,
                    price=p_R3,alpha=out3$par[1],nest_allocation=a3,
                   mu=mu_cal,
                    shares_obs = shares3)

out3b$root
delta_true

delta_cal <- out3b$root
alpha_cal <- out3$par[1]



## evaluate matching to observables


out_cal <- BBoptim(f = bertrand_foc, par = p_R3 - .4,
                   own = own_pre, alpha= alpha_cal,
                   delta = delta_cal, cost = c_j,
                   nest_allocation = a3, mu =  mu_cal,
                   sumFOC = TRUE)

price_cal <- out_cal$par
price_cal
p_R3

share_cal <- share_calc(price=p_R3,delta=delta_cal,alpha=alpha_cal,
                        nest_allocation=a3,mu=mu_cal)
share_cal
shares3

diversions_cal <- diversion_calc(price=p_R3,alpha=alpha_cal,delta=delta_cal,
                                 nest_allocation=a3,mu=mu_cal)

diversions_cal
true_div


## -----------------------------------------------------------------------------

# Create limited data objects
c_j_NA <- rep(NA, J)
c_j_NA[1] <- c_j[1]

true_div_NA <- NA * true_div
true_div_NA[1,] <- true_div[1,]


## ----results = "hide", warning=FALSE------------------------------------------

# calibration
x00 <- c(alpha_true,0.7,0.7)-.2

mu_constraints <- matrix(c(1,1,0,0,0,0,0,1,1,1), nrow = 5, ncol = 2)
wt_vector <- c(1)

out4 <- optim(f = bertrand_calibrate_gnl, par = x00,
              own = own_pre, price = p_R3,
              shares = shares3, cost  = c_j_NA,
              weight = wt_vector, nest_allocation=a3, 
              div_matrix = true_div_NA,
              mu_constraint_matrix = mu_constraints,
              control = list(maxit = 1000) )


## -----------------------------------------------------------------------------
out4$par[1]
alpha_true
out4$par[2:3]
mu3


