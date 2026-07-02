#' Nash bargaining calibration
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
#' interpretation that buyer maximizes profts.
#' @param weight Weighting matrix
#'
#' @returns The first-order conditions
#'
#' @details This function calculate the first-order conditions from a Bertrand
#' price-setting model of competition
#'
#' @examples
#' alpha  <- -0.9
#' delta <- c(.81,.93,.82)
#' c_j <- c(.05,.31,.30)
#' own_pre = diag(3)
#' p0 <- c_j*1.1
#' share1 <- (exp(delta + alpha*p0))/(1+sum(exp(delta + alpha*p0)))
#' wt_matrix <- diag(c(1,1,1,1000,1000,1000))
#'
#' bargain_calibrate(param = c(alpha,delta),own = own_pre,price = p0,
#' shares = share1,cost = c_j, weight = wt_matrix,
#' lambda = 0.5)
#'
#' @export



##################################################################
# Nash bargaining calibration
##################################################################

bargain_calibrate <- function(param,own,price,shares,cost,weight,
                               lambda,includeMUI=TRUE){

  J <- length(price)
  alpha <- param[1]
  delta <- param[2:(1+J)]

  x0 <- price
  out <- rootSolve::multiroot(f = bargain_foc, start = x0, own = own,
                   alpha= alpha, delta = delta, cost = cost,
                   lambda = lambda, includeMUI = includeMUI)
  price_m <- out$root
  share_m <- (exp(delta + alpha*price_m))/(1+sum(exp(delta + alpha*price_m)))

  pdiff <- price - price_m
  sdiff <- shares - share_m

  objfxn <- c(pdiff,sdiff) %*% weight %*% c(pdiff,sdiff)
  return(objfxn)
}

