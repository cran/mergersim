#' Downstream Bertrand first-order conditions
#'
#' @param price_r Downstream or retail prices
#' @param own_down Ownership matrix for downstream firms
#' @param own_up Ownership matrix for upstream firms
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost_r Marginal costs for downstream firm for each product
#' @param price_w Upstream or wholesale prices, treated as costs by downstream
#' firms
#' @param cost_w Marginal costs for upstream firm for each product
#' @param sumFOC logical; whether to return the sum of squares of
#' the first-order conditions. Defaults to FALSE, in which case it returns each
#' product first-order condition as a vector.
#'
#' @returns The first-order conditions
#'
#' @details This function calculates the first-order conditions from a Bertrand
#' price-setting model of competition as the downstream market of a vertical
#' supply chain. Assumes logit demand.
#'
#' @export



##################################################################
# Bertrand model first-order conditions
##################################################################
# This version is for a vertical supply chain, and accepts wholesale prices
# price_w as an input.
# This function can handle multiple integrated goods, and missing downstream costs.
# still need to add missing upstream costs (through omega_up)

bertrand_foc_vert <- function(price_r,own_down,own_up,alpha,delta,
                              cost_r,price_w,cost_w, sumFOC = FALSE){

  # first create ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )

  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )

  VI_D <- as.numeric(own_down == own_up)
  #k <- which(VI_D == 1)
  VI_D_idx <- which(own_down == own_up)

  # VI_id identifies integrated firm. Should be singleton or empty:
  #VI_id <- unique(own_down[which(own_down %in% own_up)])
  #if (length(VI_id) > 1) { stop('Function cannot handle multiple integrated firms')}
  #if (length(VI_id) == 0) { VI_id <- ""}

  #VI_U_w <- as.numeric(VI_id == own_up) - VI_D
  #VI_U_r <- as.numeric(VI_id == own_down) - VI_D
  VI_U_w <- as.numeric(own_up %in% own_down) - VI_D
  VI_U_r <- as.numeric(own_down %in% own_up) - VI_D
  # matrix form
  own_R_up <- t(sapply(own_down,
                       function(i) as.numeric(own_up %in% i) ))
  # remove integrated goods from partner profits
  own_R_up[,VI_D_idx] <- 0

  # calculate Bertrand FOCs
  shares <- (exp(delta + alpha*price_r))/(1+sum(exp(delta + alpha*price_r)))
  m <- (1-VI_D) * (price_r - price_w - cost_r) + (VI_D) * (price_r - cost_w - cost_r) #EDM Effect
  exl_na <- which(!is.na(m))

  ownd <- alpha*shares*(1-shares)
  crossd <- -alpha*shares%*%t(shares)
  dd <- crossd
  diag(dd) <- ownd
  omega <- (own_R * t(dd))

  #out <- omega %*% m + shares
  #UPP_effect <- (own_R_up * t(dd)) %*% (price_w - cost_w)

  #foc <- out + UPP_effect

  # This to deal with possible missing/NA costs
  margin_tilde <- m[exl_na]
  omega_tilde <- omega[exl_na,exl_na]

  out <- omega_tilde %*% margin_tilde + shares[exl_na]
  UPP_effect <- (own_R_up * t(dd)) %*% (price_w - cost_w)

  foc <- out + UPP_effect[exl_na]


  if (sumFOC == FALSE) {
    return(foc)
  } else {
    out <- sum(foc^2)
    return(out)
  }
}





#' Downstream Bertrand first-order conditions
#'
#' @param price_r Downstream or retail prices
#' @param own_down Ownership matrix for downstream firms
#' @param own_up Ownership matrix for upstream firms
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost_r Marginal costs for downstream firm for each product
#' @param price_w Upstream or wholesale prices, treated as costs by downstream
#' firms
#' @param cost_w Marginal costs for upstream firm for each product
#' @param nest_allocation For generalized nested logit demand, a J-by-K matrix
#' where each element (j,k) designates the membership of good j in nest k. Rows
#' should sum to 1.
#' @param mu Nesting parameters for each nest
#' @param sumFOC logical; whether to return the sum of squares of
#' the first-order conditions. Defaults to FALSE, in which case it returns each
#' product first-order condition as a vector.
#'
#' @returns The first-order conditions
#'
#' @details This function calculates the first-order conditions from a Bertrand
#' price-setting model of competition as the downstream market of a vertical
#' supply chain. This version allows for generalized nested logit demand.
#'
#' @export


bertrand_foc_vert_gnl <- function(price_r,own_down,own_up,alpha,delta,
                                  cost_r,price_w,cost_w,
                                  nest_allocation=NA, mu=NA, sumFOC = FALSE){

  J <- length(price_r)

  # If GNL, define GNL objects
  a_jk <- nest_allocation
  B <- 1*(a_jk > 0)

  # If no GNL parameters, treat as standard logit. One nest. mu=1.
  if (any(is.na(nest_allocation))) {
    K <- 1
    B <- matrix(1, ncol = 1, nrow = J)
    a_jk <- B
    mu <- rep(1,K)
  }

  # first create ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )

  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )

  VI_D <- as.numeric(own_down == own_up)
  k <- which(VI_D == 1)

  # VI_id identifies integrated firm. Should be singleton or empty:
  VI_id <- unique(own_down[which(own_down %in% own_up)])
  if (length(VI_id) > 1) { stop('Function cannot handle multiple integrated firms')}
  if (length(VI_id) == 0) { VI_id <- ""}

  VI_U_w <- as.numeric(VI_id == own_up) - VI_D
  VI_U_r <- as.numeric(VI_id == own_down) - VI_D

  # calculate Bertrand FOCs
  shares <- share_calc(price = price_r, alpha = alpha, delta = delta,
                       nest_allocation = a_jk, mu = mu)
  m <- (1-VI_D) * (price_r - price_w - cost_r) + (VI_D) * (price_r - cost_w - cost_r) #EDM Effect
  dd <- numDeriv::jacobian(share_calc, x = price_r, delta = delta, alpha = alpha,
                 nest_allocation=a_jk,mu=mu)
  omega <- (own_R * t(dd))

  # if all costs available:
  foc1 <- omega %*% m + shares

  UPP_effect <- (VI_D+VI_U_r) * (dd %*% (VI_U_w*(price_w - cost_w)) )

  foc <- foc1 + UPP_effect
  if (sumFOC == FALSE) {
    return(foc)
  } else {
    out <- sum(foc^2)
    return(out)
  }
}

