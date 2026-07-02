#' Bertrand first-order conditions
#'
#' @param price Price
#' @param own Ownership matrix
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost Marginal costs for each product
#' @param sumFOC logical; whether to return the sum of squares of
#' the first-order conditions. Defaults to FALSE, in which case it returns each
#' product first-order condition as a vector.
#'
#' @returns The first-order conditions
#'
#' @details This function calculate the first-order conditions from a Bertrand
#' price-setting model of competition
#'
#' @examples
#' TO BE ADDED.
#' @noRd



##################################################################
# Bertrand model first-order conditions
##################################################################
# This is the simplest Bertrand model. This function is undocumented
# and not exported to the package (no export tag). Just for testing
# purposes. This function can handle only simple logit, and cannot
# handle missing costs.

bertrand_foc_logit <- function(price, own, alpha, delta, cost,
                         sumFOC = FALSE){
  own_R <- own

  shares <- (exp(delta + alpha*price))/(1+sum(exp(delta + alpha*price)))

  m <- price - cost

  ownd <- alpha*shares*(1-shares)
  crossd <- -alpha*shares%*%t(shares)
  dd <- crossd
  diag(dd) <- ownd

  omega <- (own_R * t(dd))

  foc <- omega %*% m + shares

  if (sumFOC == FALSE) {
    return(foc)
  } else {
    out <- sum(foc^2)
    return(out)
  }

}

