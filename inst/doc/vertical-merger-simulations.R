## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## -----------------------------------------------------------------------------
library(mergersim)

library(BB)
library(rootSolve)
library(antitrust)
library(numDeriv)  # for jacobian() function

## -----------------------------------------------------------------------------
alpha  <- -0.9
R <- 3
W <- 2
delta <- c(0.2, 0.3, 0.9, 1.0, 0.8, 0.9)
c_R_vec <- matrix(.1, nrow = (R*W), ncol = 1)
c_W_vec <- matrix(.2, nrow = (R*W), ncol = 1)

## -----------------------------------------------------------------------------
# Define market structure/ownership
own_down_pre <- paste0("R",rep(c(1,2,3),each=2))
own_up_pre <- paste0("W",rep(c(1,2),3))
own_down_pre
own_up_pre

## -----------------------------------------------------------------------------
p_W_vec <- matrix(.25, nrow = (R*W), ncol = 1)  # dummy wholesale prices

p_r_start <- c_R_vec*1.1
out1 <- multiroot(f = mergersim:::bertrand_foc_novert ,start = p_r_start, 
                  own_down = own_down_pre, alpha= alpha, 
                  delta = delta, cost = c_R_vec,
                  price_w = p_W_vec)

p1 <- out1$root
p1


## -----------------------------------------------------------------------------

p_r_start <- as.numeric(c_R_vec*1.1)
out1b <- BBoptim(f = mergersim:::bertrand_foc_novert, par = p_r_start, 
                  own_down = own_down_pre, alpha= alpha, 
                  delta = delta, cost = c_R_vec,
                  price_w = p_W_vec, sumFOC = TRUE)

p1b <- out1b$par
p1b


## -----------------------------------------------------------------------------
lambda <-  0.5

x0 <- as.numeric(c_W_vec*1.0)

mergersim:::bargain_foc_novert_sim(price_w = x0, own_down = own_down_pre, own_up = own_up_pre, 
         alpha= alpha, delta = delta, 
         cost_w = c_W_vec, cost_r = c_R_vec, 
         lambda = lambda, price_r = (p1*2))

out2 <- BBoptim(par = x0, fn = mergersim:::bargain_foc_novert_sim, 
                own_down = own_down_pre, own_up = own_up_pre,
                alpha= alpha, delta = delta, 
                cost_w = c_W_vec, cost_r = c_R_vec, 
                lambda = lambda, price_r = (p1*1))

out2$par

p_W2 <- out2$par
shares2 <- (exp(delta + alpha*p1))/(1+sum(exp(delta + alpha*p1)))


## ----results = "hide", warning=FALSE, include=FALSE---------------------------
tol <- .0001
error <- 1

p_W0 <- matrix(.25, nrow = (R*W), ncol = 1)  # dummy wholesale prices
p_R0 <- p_r_start
  
while (error > tol) {

  # use EITHER multiroot OR BBoptim to solve for downstream prices
  # multiroot:
  # out1 <- multiroot(f = bertrand_foc_novert ,start = p_R0, 
  #                 own_down = own_down_pre, alpha= alpha, 
  #                 delta = delta, cost = c_R_vec,
  #                 p_W = p_W0)
  # 
  # p_R1 <- out1$root
  # BBoptim:
  out1 <- BBoptim(f = mergersim:::bertrand_foc_novert, par = p_R0, 
                  own_down = own_down_pre, alpha= alpha, 
                  delta = delta, cost = c_R_vec,
                  price_w = p_W0, sumFOC = TRUE)

  p_R1 <- out1$par
  
  out2 <- BBoptim(par = as.numeric(p_W0), fn = mergersim:::bargain_foc_novert_sim, 
                  own_down = own_down_pre, own_up = own_up_pre, 
                  alpha= alpha, delta = delta, 
                  cost_w = c_W_vec, cost_r = c_R_vec, lambda = lambda, 
                  price_r = p_R1)

  p_W1 <- out2$par
  
  error <- max(abs(c(p_W1-p_W0,p_R1-p_R0)))
  print(error)
  
  p_W0 <- p_W1
  p_R0 <- p_R1
}

## -----------------------------------------------------------------------------
# Check if values make sense
print(p_R1)
print(p_W1)

shares1 <- (exp(delta + alpha*p_R1))/(1+sum(exp(delta + alpha*p_R1)))
as.numeric(shares1)
sum(shares1)

# GFT
mergersim:::bargain_foc_novert_sim(price_w = p_W1, own_down = own_down_pre, own_up = own_up_pre, 
         alpha= alpha, 
         delta = delta, cost_w = c_W_vec, cost_r = c_R_vec, lambda = lambda, 
         price_r = p_R1,
         returnGFT = TRUE)

# wholesaler profits
#as.numeric(own_W_pre %*% ((p_W1 - c_W_vec)*shares1) )

# retailer profits
#as.numeric(own_R_pre %*% ((p_R1 - p_W1 - c_R_vec)*shares1) )

## -----------------------------------------------------------------------------

# Create post ownership
own_up_post <- own_up_pre
own_down_post <- own_down_pre
own_down_post[own_down_post == "R1"] <- "W1"

bertrand_foc_vert(price_r = p_R1,own_down=own_down_post,own_up=own_up_post,
          alpha=alpha,delta=delta,
          cost_r = c_R_vec, price_w = p_W1, cost_w = c_W_vec)


bargain_foc_vert_sim(price_w = p_W1,own_down=own_down_post,own_up=own_up_post,
          alpha=alpha,delta=delta, cost_w =c_W_vec,
          cost_r =c_R_vec, lambda=0.5, price_r = p_R1)


# Test that foc_vert gives same result as before when VI = 0.
mergersim:::bargain_foc_novert_sim(price_w = p_W1, own_down = own_down_pre, own_up = own_up_pre, 
         alpha= alpha, delta = delta, cost_w = c_W_vec, 
         cost_r = c_R_vec, lambda = lambda, price_r = p_R1)

bargain_foc_vert_sim(price_w = p_W1,own_down=own_down_pre,own_up=own_up_pre,
          alpha=alpha,delta=delta, cost_w =c_W_vec,
          cost_r =c_R_vec, lambda=0.5, price_r = p_R1)

# And for upstream
bertrand_foc_vert(price_r = p_R1, own_down = own_down_pre, own_up = own_up_pre,
          alpha = alpha, delta = delta,
          cost_r = c_R_vec, price_w = p_W1, 
          cost_w = c_W_vec)

mergersim:::bertrand_foc_novert(price_r = p_R1, own_down = own_down_pre, alpha = alpha, 
         delta = delta, cost = c_R_vec, price_w = p_W1)
  

## ----results = "hide", warning=FALSE------------------------------------------
J <- length(p_R1)
alpha_start <- -1
delta_start <- rep(0.5,J)
x00b <- -1.1
wt_matrix <- diag(c(rep(1,J),rep(10,J)))

c_j <- c_R_vec  # eventually remove c_j from code

out3 <- BBoptim(f = bertrand_vert_calibrate, par = x00b, 
                   own_down = own_down_pre, price = p_R1, 
                   shares = shares1, cost  = c_j, price_w = p_W1)



## -----------------------------------------------------------------------------

alpha3 <- out3$par
delta3 <- log(shares1) - log(1-sum(shares1)) - alpha3*p_R1

# recover true parameters
alpha
delta
alpha3
delta3


## -----------------------------------------------------------------------------
# note that this optimization recovers the true demand parameters and shares

shares3 <- (exp(delta3 + alpha3*p_R1))/(1+sum(exp(delta3 + alpha3*p_R1)))
shares3
as.numeric(shares1)

## ----results = "hide", warning=FALSE------------------------------------------
lambda_start <- 0.4   # starting value for lambda 

out4 <- BBoptim(f = bargain_vert_sim_calibrate, par = lambda_start, 
                price_w = p_W1,own_down = own_down_pre,
                own_up =own_up_pre,alpha=alpha3,delta=delta3,
                cost_w = c_W_vec, cost_r = c_R_vec, price_r = p_R1,
                lower = 0, upper = 1)

## -----------------------------------------------------------------------------
lambda_end <- out4$par
lambda_end

## -----------------------------------------------------------------------------
# Set some costs to NA -- (this line not fully general yet)
c_R_vec_NA <- c(c_R_vec[own_down_pre == "R1"], NA, NA, NA, NA)


x00b <- -1.1
bertrand_vert_calibrate(param = x00b, 
                   own_down = own_down_pre, price = p_R1, 
                   shares = shares1, cost  = c_R_vec_NA, price_w = p_W1 )


out3d <- BBoptim(f = bertrand_vert_calibrate, par = x00b, 
                   own_down = own_down_pre, price = p_R1, 
                   shares = shares1, cost  = c_R_vec_NA, price_w = p_W1)

alpha3d <- out3d$par
delta3d <- log(shares1) - log(1-sum(shares1)) - alpha3d*p_R1

# still recover true parameters
alpha
delta
alpha3d
delta3d


## ----results = "hide", warning=FALSE------------------------------------------


x00b <- rep(1,J)
mergersim:::bertrand_vert_calibrate_costs(param = x00b, 
                   own_down = own_down_pre, price = p_R1, 
                   shares = shares1, alpha = alpha, delta = delta,
                   price_w = p_W1 )


out3e <- BBoptim(f = mergersim:::bertrand_vert_calibrate_costs, par = x00b, 
                   own_down = own_down_pre, price = p_R1, 
                   shares = shares1, alpha = alpha, delta = delta,
                   price_w = p_W1)

## -----------------------------------------------------------------------------
out3e$par
c(c_R_vec)


## -----------------------------------------------------------------------------
shares1
m_downstream <- p_R1 - c_R_vec - p_W1

cor(shares1,m_downstream)
# plot(shares1,m_downstream)


## ----results = "hide", warning=FALSE------------------------------------------
# ( Ideally, first recover pre-merger prices, using VI functions.)

bertrand_foc_vert(price_r=p_R1,own_down=own_down_post,own_up=own_up_post,
          alpha=alpha,delta=delta,
          cost_r=c_R_vec,price_w=p_W1, cost_w=c_W_vec)

bargain_foc_vert_sim(price_w=p_W1,own_down = own_down_post,own_up=own_up_post,
          alpha=alpha,delta=delta, cost_w=c_W_vec,
          cost_r=c_R_vec, lambda=0.5, price_r = p_R1)

tol <- .0001
error <- 1

p_W0_post <- matrix(.25, nrow = (R*W), ncol = 1)  # dummy wholesale prices
p_R0_post <- p_r_start
  
while (error > tol) {

  # Use EITHER multiroot or BBoptim
  # multiroot:
  # out1 <- multiroot(f = bertrand_foc_vert ,start = p_R0_post, 
  #                 own_down = own_down_post, own_up = own_up_post,
  #                 alpha = alpha3, delta = delta3, c_R = c_R_vec, 
  #                 p_W = p_W0_post, c_W = c_W_vec)
  # 
  # p_R1_post <- out1$root
  # BBoptim:
  out1 <- BBoptim(f = bertrand_foc_vert, par = p_R0_post, 
                  own_down = own_down_post, own_up = own_up_post,
                  alpha = alpha3, delta = delta3, cost_r = c_R_vec, 
                  price_w = p_W0_post, cost_w = c_W_vec, sumFOC = TRUE)

  p_R1_post <- out1$par
  
  out2 <- BBoptim(par = as.numeric(p_W0_post), fn = bargain_foc_vert_sim, 
                  own_down = own_down_post, own_up = own_up_post, 
                  alpha= alpha3, delta = delta3, 
                  cost_w = c_W_vec, cost_r = c_R_vec, lambda = lambda_end, 
                  price_r = p_R1_post)
  
  p_W1_post <- out2$par
  
  error <- max(abs(c(p_W1_post-p_W0_post,p_R1_post-p_R0_post)))
  print(error)
  
  p_W0_post <- p_W1_post
  p_R0_post <- p_R1_post
}

## -----------------------------------------------------------------------------
p_R1
p_W1
p_R1_post
p_W1_post

# retail price change with vertical integration
(p_R1_post-p_R1)/p_R1


## -----------------------------------------------------------------------------
shareDown <- as.numeric(shares1)
priceDown <- as.numeric(p_R1)
ownerPreDown <- own_down_pre
priceUp <- as.numeric(p_W1)
ownerPreUp <- own_up_pre
priceOutSide <- 0

# Downstream margin as a percent
marginDown <- (priceDown - c_R_vec - priceUp)/ priceDown
# Upstream margin as a percent
marginUp <- (priceUp - c_W_vec) / priceUp

## Simulate a vertical merger
ownerPostUp <- own_up_post
ownerPostDown <- own_down_post

simres_vert <- vertical.barg(sharesDown =shareDown,
                             pricesDown = priceDown,
                             marginsDown = as.numeric(marginDown),
                             ownerPreDown = ownerPreDown,
                             ownerPostDown = ownerPostDown,
                             pricesUp = priceUp,
                             marginsUp = as.numeric(marginUp),
                             ownerPreUp = ownerPreUp,
                             ownerPostUp = ownerPostUp,
                             priceOutside = priceOutSide)

summary(simres_vert)

simres_vert@down@slopes

(simres_vert@down@pricePost - simres_vert@down@pricePre) / simres_vert@down@pricePre

