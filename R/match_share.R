#' Market share calculations
#'
#' Calculates choice probabilities based on logit or generalized nested logit
#' demand parameters
#'
#' @param price Price
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param shares_obs Observed market shares
#' @param nest_allocation For generalized nested logit demand, a J-by-K matrix
#' where each element (j,k) designates the membership of good j in nest k. Rows
#' should sum to 1.
#' @param mu Nesting parameters for each nest
#'
#' @returns Returns vector of difference between predicted shares and observed
#' shares
#'
#' @details This function calculates the difference between model predicted
#' choice probabilities and observed market shares.
#'
#' @examples
#' match_share(price=c(2.1,2.4,2.1), delta=c(.81,.93,.82), alpha=-0.9,
#' shares_obs = c(.2, .2, .2))
#'
#' @export


##################################################################
# Match shares
##################################################################


match_share <- function(price, delta, alpha, nest_allocation=NA, mu=NA,
                        shares_obs){

  out <- share_calc(price=price, delta=delta, alpha=alpha,
                    nest_allocation=nest_allocation, mu=mu) -
    shares_obs

  return(out)
}



#' Market share calculations
#'
#' Calculates choice probabilities based on logit or generalized nested logit
#' demand parameters
#'
#' @param price Price
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param shares Observed shares to match
#'
#' @returns Returns vector of difference between predicted shares and observed
#' shares
#'
#' @details This function calculates the difference between model predicted
#' choice probabilities and observed market shares.
#'
#' @examples
#' TO BE ADDED.
#'

#' @noRd


##################################################################
# Match shares
##################################################################


match_share_logit <- function(delta,shares,alpha,price){
  output <- shares - (exp(delta + alpha*price))/(1+sum(exp(delta + alpha*price)))
  return(output)
}

