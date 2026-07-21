#' Second score bargaining first-order conditions
#'
#' @param price Price
#' @param own Ownership matrix
#' @param param Price coefficient alpha and mean values delta, parameters to
#'  calibrate
#' @param shares Observed market shares
#' @param cost Marginal costs for each product
#' @param lambda Bargaining power of the buyer
#' @param includeMUI logical; whether to include marginal utility of income
#' in buyer's payoff, thereby translating dollars to utility. Default is True,
#' interpreted as buyer maximizing utility. Setting equal to False would have
#' interpretation that buyer maximizes profits.
#' @param weight Weighting vector of length 2*J
#'
#' @returns The first-order conditions
#'
#' @details This function calculate the first-order conditions from a bargaining
#' model that nests the second score auction.
#'
#' @examples
#' alpha  <- -0.9
#' delta <- c(.81,.93,.82)
#' own_pre = diag(3)
#' p0 <- c(.05, .34, .33)
#' c_j <- c(.05,.31,.30)
#' wt_vector <- c(1,1,1,1000,1000,1000)
#' share1 <- c( 0.31, 0.27, 0.25)
#'
#' ssbargain_calibrate(param = c(alpha,delta),own = own_pre,price = p0,
#'                     shares = share1, cost = c_j, weight = wt_vector,
#'                     lambda = 0.5)
#'
#' @export


##################################################################
# Second score bargaining first-order conditions
##################################################################

ssbargain_calibrate <- function(param,own,price,shares,cost,weight = NA,
                                   lambda,includeMUI=TRUE){

  J <- length(price)

  if (anyNA(weight)) {
    weight <- c(rep(1, times = J), rep(1000, times = J))
  }

  alpha <- param[1]
  delta <- param[2:(1+J)]

  x0 <- price
  out <- rootSolve::multiroot(f = ssbargain_foc,start = x0, own = own,
                   alpha= alpha, delta = delta, cost = cost,
                   lambda = lambda, includeMUI = includeMUI)

  price_m <- out$root
  share_m <- (exp(delta + alpha*cost))/(1+sum(exp(delta + alpha*cost)))

  pdiff <- price - price_m
  sdiff <- shares - share_m

  objfxn <- c(pdiff,sdiff) %*% diag(weight) %*% c(pdiff,sdiff)
  return(objfxn)
}

