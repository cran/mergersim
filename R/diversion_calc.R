#' Market share calculations
#'
#' Calculates choice probabilities based on logit or generalized nested logit
#' demand parameters
#'
#' @param price Price
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param nest_allocation For generalized nested logit demand, a J-by-K matrix
#' where each element (j,k) designates the membership of good j in nest k. Rows
#' should sum to 1.
#' @param mu Nesting parameters for each nest
#' @param marginal logical; if True, diversions are calculated as diversions from
#' a small price increase. If False, diversions are calculated as second choice
#' diversions.
#' @param outsideOption logical; whether to include an outside option in choice
#' set. Default is TRUE.
#'
#' @returns Returns vector of choice probabilities for each good
#'
#' @details This function calculates choice probabilities based on demand
#' parameters for a logit or generalized nested logit demand system
#'
#' @examples
#'
#' diversion_calc(price=c(2.1,2.4,2.1),alpha=-0.9,delta=c(.81,.93,.82))
#'
#' @export


##################################################################
# Calculate choice probabilities
##################################################################

diversion_calc <- function(price,alpha,delta,nest_allocation=NA,mu=NA,
                                           marginal = FALSE,
                           outsideOption = TRUE){

  # Define GNL objects
  a_jk <- nest_allocation
  B <- 1*(a_jk > 0)

  obs_share <- share_calc(price = price, alpha = alpha, delta = delta,
                          nest_allocation = a_jk, mu = mu,
                          outsideOption = outsideOption)

  J <- length(delta)
  div_out <- matrix(0, nrow = J, ncol = J)

  if (marginal == FALSE) {
    for (i in 1:J) {
      delta_cf <- delta
      delta_cf[i] <- -Inf
      cf_sharei <- share_calc(price = price, alpha = alpha, delta = delta_cf,
                              nest_allocation = a_jk, mu = mu,
                              outsideOption = outsideOption)
      div_out[i,] <- (cf_sharei - obs_share)/obs_share[i]
    }

    diag(div_out) <- 0
    return(div_out)
  }

  if (marginal == TRUE) {
    for (i in 1:J) {
      p_cf <- price
      p_cf[i] <- p_cf[i] * 1.01   # increase price by 1%
      cf_sharei <- share_calc(price = p_cf, alpha = alpha, delta = delta,
                              nest_allocation = a_jk, mu = mu,
                              outsideOption = outsideOption)
      div_out[i,] <- (cf_sharei - obs_share)/(obs_share[i] - cf_sharei[i])
    }

    diag(div_out) <- 0
    return(div_out)
  }

}

