#' Downstream Bertrand first-order conditions
#'
#' @param price_r Downstream or retail prices
#' @param own_down Ownership matrix for downstream firms
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost Marginal costs for each product
#' @param price_w Upstream or wholesale prices, treated as costs by downstream
#' firms
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
#' @noRd



##################################################################
# Bertrand model first-order conditions
##################################################################
# This version is for a vertical supply chain, and thus accepts wholesale prices
# price_w as an input.
# This function is not exported, bertrand_foc_vert is a generalization.

bertrand_foc_novert <- function(price_r,own_down,alpha,delta,cost,price_w,
                                sumFOC = FALSE){
  # first create ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  ownership <- t(sapply(own_down, own_fun_down) )

  # then calculate foc
  shares <- (exp(delta + alpha*price_r))/(1+sum(exp(delta + alpha*price_r)))
  m <- price_r - cost - price_w
  exl_na <- which(!is.na(m))

  ownd <- alpha*shares*(1-shares)
  crossd <- -alpha*shares%*%t(shares)
  dd <- crossd
  diag(dd) <- ownd
  omega <- (ownership * t(dd))

  # This to deal with possible missing/NA costs
  margin_tilde <- m[exl_na]
  omega_tilde <- omega[exl_na,exl_na]

  foc <- omega_tilde %*% margin_tilde + shares[exl_na]

  if (sumFOC == FALSE) {
    return(foc)
  } else {
    out <- sum(foc^2)
    return(out)
  }
}


