#' Second score auction calibration
#'
#' @param param The price coefficient alpha to be calibrated
#' @param own Ownership matrix
#' @param price Price
#' @param share Share
#' @param cost Marginal costs for each product
#' @param weight Weighting vector of length equal to number of observed costs
#'
#' @returns The first-order conditions
#'
#' @details This function calculate the first-order conditions from a second
#' score auction model of competition
#'
#' @examples
#' own_pre = diag(3)
#' p0 <- c(.05, .34, .33)
#' share1 <- c( 0.31, 0.27, 0.25)
#' c_j <- c(.05,.31,.30)
#' wt_vector <- c(1,1,1)
#'
#' ssa_calibrate(param = -1,own = own_pre,price=p0,share=share1,cost=c_j,
#'               weight = wt_vector)
#'
#' @export



##################################################################
# Second score auction calibration
##################################################################

ssa_calibrate <- function(param,own,price,share,cost,weight = NA){

  num_c <- sum(!is.na(cost))

  if (anyNA(weight)) {
    weight <- rep(1, times = num_c)
  } else if (length(weight) != num_c) {
    stop("weight must have length equal to the number of margins (", num_c, ")")
  }

  alpha <- param[1]

  p_hat <- cost + log(1 - own%*%share) / (alpha*own%*%share)

  FOC <- p_hat - price
  FOC <- stats::na.omit(FOC)

  objfxn <- c(FOC) %*% diag(weight) %*% c(FOC)

  return(objfxn)
}
