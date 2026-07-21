#' Bertrand model calibration
#'
#' @param param Price coefficient alpha parameter to calibrate
#' @param price Price
#' @param shares Observed market shares
#' @param own Ownership matrix
#' @param cost Marginal costs for each product
#' @param weight Weighting vector of length equal to number of margins provided
#' @param returnOutcomes logical; should equilibrium objects be returned (mean
#' value parameter, prices, shares, costs) as a list.
#'
#' @returns Distance between observed values and model predicted values for
#' prices and shares
#'
#' @details This function calculate the first-order conditions from a Bertrand
#' price-setting model of competition. This function is only for standard logit
#' demand. For nested logit or generalized nested logit, see
#' bertrand_calibrate_gnl().
#'
#' @examples
#' alpha  <- -0.9
#' delta <- c(.81,.93,.82)
#' c_j <- c(.05,.31,.30)
#'
#' own_pre = diag(3)
#'
#' p0 <- c_j*1.1
#'
#' share1 <- (exp(delta + alpha*p0))/(1+sum(exp(delta + alpha*p0)))
#' x00 <- c(-1)
#' wt_vector <- c(1,1,1)
#'
#' bertrand_calibrate( param = x00,
#'                     own = own_pre, price = p0,
#'                     shares = share1, cost  = c_j,
#'                     weight = wt_vector)
#'
#' @export



##################################################################
# Bertrand model calibration
##################################################################

bertrand_calibrate <- function(param, own, price, shares, cost, weight = NA,
                                  returnOutcomes = FALSE) {

  num_c <- sum(!is.na(cost))

  if (anyNA(weight)) {
    weight <- rep(1, times = num_c)
  } else if (length(weight) != num_c) {
    stop("weight must have length equal to the number of margins (", num_c, ")")
  }

  alpha <- param[1]

  delta <- log(shares) - log(1 - sum(shares)) - alpha * price

  dd <- numDeriv::jacobian(share_calc, x = price,
                           delta = delta, alpha = alpha,
                           nest_allocation = matrix(1, nrow = length(price), ncol = 1),
                           mu = 1)
  omega <- own * t(dd)

  # Back out implied cost from FOC at observed prices: Omega*(p-c) + s = 0
  # => c_implied = p + solve(Omega) %*% s
  c_implied <- as.numeric(price + solve(omega) %*% shares)

  # Only compute residuals where cost is observed
  cdiff <- cost - c_implied
  keep <- !is.na(cdiff)

  if (returnOutcomes == FALSE) {
    objfxn <- sum(weight * cdiff[keep]^2)
    return(objfxn)
  }

  if (returnOutcomes == TRUE) {
    return(list("cost_cal" = c_implied,
                "delta_cal" = delta) )
  }

}
