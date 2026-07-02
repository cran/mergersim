#' Calibrate lambda for Nash bargain in vertical model with simultaneous timing
#'
#' @param lambda Bargaining power of the buyer/retailer
#' @param price_w Upstream or wholesale prices
#' @param own_down Ownership matrix for downstream firms
#' @param own_up Ownership matrix for upstream firms
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost_r Marginal costs for downstream firm for each product
#' @param cost_w Marginal costs for upstream firm for each product
#' @param price_r Retail prices
#'
#' @returns The value of objective function
#'
#' @details This function can be used to calibrate the bargaining parameter
#' in a Nash bargain which is the upstream market of a vertical supply chain.
#' Assumes logit demand.
#'
#' @examples
#' bargain_vert_sim_calibrate(lambda = 0.4,
#' price_w = c(1.6, 1.6, 1.6, 1.6, 1.6, 1.6),
#' own_down = paste0("R",rep(c(1,2,3),each=2)),
#' own_up = paste0("W",rep(c(1,2),3)),
#' alpha = -0.9,
#' delta = c(0.2, 0.3, 0.9, 1.0, 0.8, 0.9),
#' cost_w = rep(.2, times = 6),
#' cost_r = rep(.1, times = 6),
#' price_r = c(2.9, 2.9, 3.0, 3.0, 3.0, 3.0))
#'
#' @export



##################################################################
# Nash Product
##################################################################

## A bargaining calibration function for simple logit demand and calibrates
## lambda. Uses bargain_foc_novert_sim().

bargain_vert_sim_calibrate <- function(lambda,price_w,own_down,own_up,
                                       alpha,delta,cost_w,cost_r,price_r){

  out <- BB::BBoptim(par = as.numeric(price_w), fn = bargain_foc_novert_sim,
                 own_down = own_down,
                 own_up = own_up, alpha= alpha, delta = delta,
                 cost_w = cost_w, cost_r = cost_r, lambda = lambda,
                 price_r = price_r)

  p_W1 <- out$par
  pdiff <- price_w - p_W1

  objfxn <- c(pdiff) %*% diag(length(pdiff)) %*% c(pdiff)
  return(objfxn)
}


