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
#' @param returnLogsum logical; whether to return the denominator of the choice
#' probabilities (also known as the log-sum term). Defaults to FALSE, in which
#' case the function returns a vector with each product's choice probability.
#' @param outsideOption logical; whether to include an outside option in choice
#' set. Default is TRUE.
#'
#' @returns Returns vector of choice probabilities for each good
#'
#' @details This function calculates choice probabilities based on demand
#' parameters for a logit or generalized nested logit demand system
#'
#' @examples
#' share_calc(price=c(2.1,2.4,2.1),delta=c(.81,.93,.82),alpha=-0.9)
#'
#' @export



##################################################################
# Calculate choice probabilities
##################################################################

share_calc <- function(price,delta,alpha,nest_allocation=NA,mu=NA,
                           returnLogsum=FALSE, outsideOption = TRUE){
  # a is a J-by-K matrix of allocation parameters
  # B is a J-by-K matrix of indicators designating nests
  # mu is a vector length K of nesting parameters

  J <- length(price)

  # If GNL, define GNL objects
  a_jk <- nest_allocation
  B <- 1*(a_jk > 0)

  # If no GNL parameters, treat as standard logit. One nest with mu=1.
  if (any(is.na(nest_allocation))) {
    K <- 1
    B <- matrix(1, ncol = 1, nrow = J)
    a_jk <- B
    mu <- rep(1,K)
  }


  # check that length(p) == dim(B)[1] == J
  J <- dim(B)[1]
  K <- dim(B)[2]

  V_j <- delta + alpha * price

  temp1 <- a_jk * matrix(rep(exp(V_j),K), ncol = K, nrow = J)
  temp2 <- temp1 ^ matrix(rep((1/mu),J), ncol = K, nrow = J, byrow = TRUE)
  P_k_num <- colSums(temp2) ^ mu  # P_k_num should be length K
  if (outsideOption == TRUE) {
    P_k_denom <- 1+sum(P_k_num)  # WITH outside option
  } else {
    P_k_denom <- sum(P_k_num)  # NO outside option
    }

  P_k <- P_k_num / P_k_denom
  # P_k should be a length K vector that sums to 1

  P_j_Bk <- temp2 / matrix(rep(colSums(temp2),J), ncol = K, nrow = J, byrow = TRUE)
  # this should be a JxK matrix with columns that sum to 1, unless empty nest
  # fix for empty nests. Set to zero, just cannot be NaN:
  P_j_Bk[is.na(P_j_Bk)] <- 0

  P_j <- rowSums(P_j_Bk * matrix(rep((P_k),J), ncol = K, nrow = J, byrow = TRUE))
  # P_j should be a vector of length J that sums to one.
  # Add in check for both those conditions

  if (returnLogsum == FALSE) {
    out <- P_j
    return(out)
  }
  if (returnLogsum == TRUE) {
    out <- log(P_k_denom)
    return(out)
  }

}




#' Market share calculations - logit
#'
#' Calculates choice probabilities based on logit
#' demand parameters
#'
#' @param price Price
#' @param alpha Price coefficient
#' @param delta Mean values
#'
#' @returns Returns vector of choice probabilities for each good
#'
#' @details This function calculates choice probabilities based on demand
#' parameters for a logit demand system
#'
#' @noRd


## simplest version of share calculation for logit demand with outside option
## not exported nor documented, this version is for convenience only.

share_calc_logit <- function(price,delta,alpha){
  out <- (exp(delta + alpha*price))/(1+sum(exp(delta + alpha*price)))
  return(out)
}
