#' Calibrate Bertrand model demand parameters
#'
#' @param param Costs
#' @param own_down Ownership matrix for downstream firms
#' @param price Downstream prices
#' @param shares Market shares
#' @param alpha Price coefficient
#' @param delta Mean values for each product
#' @param price_w Upstream or wholesale prices, treated as costs by downstream
#' firms
#'
#' @returns Objective function value
#'
#' @details This function can be used to calibrate the costs of a
#' Bertrand model of competition that is the downstream market of a vertical
#' supply chain. Wholesale prices are treated as costs by the downstream firms.
#' Assumes logit demand.
#'
#' @noRd



##################################################################
# Bertrand model downstream calibration
##################################################################
# Calibrates the downstream Bertrand model of a vertical supply chain. Accounts
# for price_w, which is treated as a cost by downstream firms.

bertrand_vert_calibrate_costs <- function(param,own_down,price,shares,alpha,
                                          delta,price_w){

  J <- length(price)
  cost <- param

  s_0 <- 1 - sum(shares)
  delta_j <- log(shares) - log(s_0) - alpha*price

  x0 <- price
  out1 <- BB::BBoptim(f = bertrand_foc_novert, par = x0,
                  own_down = own_down, alpha= alpha,
                  delta = delta_j, cost = cost,
                  price_w = price_w, sumFOC = TRUE)

  p1 <- out1$par
  share1 <- (exp(delta_j + alpha*p1))/(1+sum(exp(delta_j + alpha*p1)))

  pdiff <- price - p1
  # sdiff <- shares - share1

  objfxn <- c(pdiff) %*% diag(J) %*% c(pdiff)
  return(objfxn)
}
