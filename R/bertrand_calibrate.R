#' Bertrand model calibration
#'
#' @param param Price coefficient alpha parameter to calibrate
#' @param price Price
#' @param shares Observed market shares
#' @param own Ownership matrix
#' @param cost Marginal costs for each product
#' @param weight Weighting matrix
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
#' wt_matrix <- diag(c(1,1,1,1000,1000,1000))
#'
#' bertrand_calibrate( param = x00,
#'                     own = own_pre, price = p0,
#'                     shares = share1, cost  = c_j,
#'                     weight = wt_matrix)
#'
#' @export



##################################################################
# Bertrand model calibration
##################################################################

bertrand_calibrate <- function(param,own,price,shares,cost,weight,
                               returnOutcomes = FALSE){

  J <- length(price)
  alpha <- param[1]

  delta <- log(shares) - log(1-sum(shares)) - alpha*price

  x0 <- price
  out1 <- BB::BBoptim(f = bertrand_foc, par = x0,
                  own = own, alpha = alpha,
                  delta = delta, cost = cost,
                  sumFOC = TRUE)

  p_model <- out1$par
  share_m <- (exp(delta + alpha*p_model))/(1+sum(exp(delta + alpha*p_model)))

  pdiff <- price - p_model
  sdiff <- shares - share_m

  if (returnOutcomes == FALSE) {
    objfxn <- c(pdiff,sdiff) %*% weight %*% c(pdiff,sdiff)
    return(objfxn)
  }

  if (returnOutcomes == TRUE) {

    meancost <- mean(cost, na.rm = TRUE)
    x00 <- ifelse(!is.na(cost), cost, meancost)

    out_cost <- stats::optim(f = bertrand_foc_c, par = x00,
                      price = price, own = own, alpha = alpha,
                      delta = delta, sumFOC = TRUE,
                      control = list(maxit = 1500) )

    cost_cal <- out_cost$par

    if (out_cost$convergence != 0) {
      warning(paste0("Cost calibration did not converge with code "),
              out_cost$convergence)
    }


    return(list("FOCs" = c(pdiff,sdiff),
                "delta_cal" = delta,
                "p_model" = p_model,
                "share_m" = share_m,
                "cost_cal" = cost_cal) )
  }

}
