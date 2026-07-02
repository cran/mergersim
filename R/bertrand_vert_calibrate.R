#' Calibrate Bertrand model demand parameters
#'
#' @param param Price coefficient, alpha
#' @param own_down Ownership matrix for downstream firms
#' @param price_r Downstream prices
#' @param shares Market shares
#' @param cost Marginal costs of downstream firm
#' @param price_w Upstream or wholesale prices, treated as costs by downstream
#' firms
#'
#' @returns Objective function value
#'
#' @details This function can be used to calibrate the demand parameters of a
#' Bertrand model of competition that is the downstream market of a vertical
#' supply chain. Wholesale prices are treated as costs by the downstream firms.
#' Assumes logit demand.
#'
#' @export



##################################################################
# Bertrand model downstream calibration
##################################################################
# Calibrates the downstream Bertrand model of a vertical supply chain. Accounts
# for price_w, which is treated as a cost by downstream firms.

bertrand_vert_calibrate <- function(param,own_down,price_r,shares,cost,price_w){

  J <- length(price_r)
  alpha <- param[1]

  s_0 <- 1 - sum(shares)
  delta_j <- log(shares) - log(s_0) - alpha*price_r

  x0 <- price_r
  out1 <- BB::BBoptim(f = bertrand_foc_novert, par = x0,
                  own_down = own_down, alpha= alpha,
                  delta = delta_j, cost = cost,
                  price_w = price_w, sumFOC = TRUE,
                  quiet = TRUE)

  p1 <- out1$par
  share1 <- (exp(delta_j + alpha*p1))/(1+sum(exp(delta_j + alpha*p1)))

  pdiff <- price_r - p1
  # sdiff <- shares - share1

  objfxn <- c(pdiff) %*% diag(J) %*% c(pdiff)
  return(objfxn)
}



#' Calibrate costs in downstream Bertrand model
#'
#' @param param Price coefficient, alpha
#' @param own_down Ownership matrix for downstream firms
#' @param price Downstream prices
#' @param shares Market shares
#' @param cost Marginal costs of downstream firm
#' @param price_w Upstream or wholesale prices, treated as costs by downstream
#' firms
#'
#' @returns Objective function value
#'
#' @details This function can be used to calibrate downstream costs in a
#' Bertrand model of competition that is the downstream market of a vertical
#' supply chain. Wholesale prices are treated as costs by the downstream firms.
#'
#' @noRd

bertrand_vert_calibrate_gnl <- function(cost, price_r,own_down,alpha,delta,price_w,
                                        nest_allocation=NA,mu=NA,sumFOC = FALSE){

  J <- length(price_r)

  # If GNL, define GNL objects
  a_jk <- nest_allocation
  B <- 1*(a_jk > 0)

  # If no GNL parameters, treat as standard logit. One nest with mu=1.
  if (any(is.na(nest_allocation))) {
    K <- 1
    B <- matrix(1, ncol = 1, nrow = J)
    a_jk <- B
    mu <- rep(1,K)
  }

  # first create ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )

  # then calculate foc
  shares <- share_calc(price=price_r,delta = delta, alpha = alpha,
                       nest_allocation=a_jk,mu=mu)
  m <- price_r - cost - price_w
  dd <- numDeriv::jacobian(share_calc, x = price_r, delta = delta, alpha = alpha,
                 nest_allocation=a_jk,mu=mu)
  omega <- (own_R * t(dd))

  # if all costs aviailable:
  # foc <- omega %*% m + shares

  # To deal with possible missing/NA costs
  exl_na <- which(!is.na(m))
  margin_tilde <- m[exl_na]
  omega_tilde <- omega[exl_na,exl_na]
  shares_tilde <- shares[exl_na]

  foc <- omega_tilde %*% margin_tilde + shares_tilde

  if (sumFOC == FALSE) {
    return(foc)
  } else {
    out <- sum(foc^2)
    return(out)
  }

}


#' Calibrate downstream Bertrand model
#'
#' @noRd

## A bertrand calibration function that accounts for vertical structure,
## calibrates all deltas, and requires all costs be provided. Uses Bert_foc().
## compared to bertrand_vert_calibrate, v2 takes weight as input, calibrates
## delta rather than use berry inversion to get it, uses multiroot to solve for
## prices rather than BBoptim, and the objective function includes prices and
## shares rather than only prices.

bertrand_vert_calibrate_v2 <- function(param,own_down,price,shares,cost,p_W,weight){

  J <- length(price)
  alpha <- param[1]
  delta <- param[2:(1+J)]

  x0 <- price
  out1 <- rootSolve::multiroot(f = bertrand_foc_novert ,start = x0,
                    own_down = own_down, alpha= alpha,
                    delta = delta, cost = cost, price_w = p_W)

  p1 <- out1$root
  share1 <- (exp(delta + alpha*p1))/(1+sum(exp(delta + alpha*p1)))

  pdiff <- price - p1
  sdiff <- shares - share1

  objfxn <- c(pdiff,sdiff) %*% weight %*% c(pdiff,sdiff)
  return(objfxn)
}


