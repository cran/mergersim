#' Bertrand first-order conditions
#'
#' @param price Price
#' @param own Ownership matrix
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost Marginal costs for each product
#' @param nest_allocation For generalized nested logit demand, a J-by-K matrix
#' where each element (j,k) designates the membership of good j in nest k. Rows
#' should sum to 1.
#' @param mu Nesting parameters for each nest
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
#' alpha  <- -0.9
#' delta <- c(.81,.93,.82)
#' c_j <- c(.05,.31,.30)
#'
#' own_pre = diag(3)
#'
#' p0 <- c_j*1.1
#'
#' bertrand_foc(price = p0,
#'              own = own_pre, alpha= alpha,
#'              delta = delta, cost = c_j)
#'
#' @export bertrand_foc



##################################################################
# Bertrand model first-order conditions
##################################################################

bertrand_foc <- function(price, own, alpha, delta, cost,
                         nest_allocation=NA, mu=NA, sumFOC = FALSE){

  J <- length(price)

  # If GNL, define GNL objects
  a_jk <- nest_allocation
  B <- 1*(a_jk > 0)

  # If no GNL parameters, treat as standard logit. One nest. mu=1.
  if (any(is.na(nest_allocation))) {
    K <- 1
    B <- matrix(1, ncol = 1, nrow = J)
    a_jk <- B
    mu <- rep(1,K)
  }

  # first create ownership matrices
  #own_fun_down <- function(x) {as.numeric(x == own_down)}
  #own_R <- t(sapply(own_down, own_fun_down) )
  # Since ownership is given by matrix directly:
  own_R <- own

  # then calculate foc
  shares <- share_calc(price = price, alpha = alpha, delta = delta,
                       nest_allocation = a_jk, mu = mu)
  m <- price - cost
  dd <- numDeriv::jacobian(share_calc, x = price, delta = delta, alpha = alpha,
                 nest_allocation=a_jk, mu=mu)
  omega <- (own_R * t(dd))

  # if all costs available:
  # foc <- omega %*% m + shares

  # To deal with possible missing/NA costs
  exl_na <- which(!is.na(m))
  margin_tilde <- m[exl_na]
  omega_tilde <- omega[exl_na,exl_na]
  shares_tilde <- shares[exl_na]

  foc <- omega_tilde %*% margin_tilde + shares_tilde

  if (sumFOC == FALSE) {
    return(foc)
  } else {
    out <- sum(foc^2)
    return(out)
  }

}



#' @rdname bertrand_foc
#' @export bertrand_foc_c

##################################################################
# Second version of function that puts cost as first parameter
# Used for backing out costs in bertrand_calibrate_gnl().
##################################################################

bertrand_foc_c <- function(cost, price, own, alpha, delta,
                           nest_allocation=NA, mu=NA, sumFOC = FALSE){

  J <- length(price)

  # If GNL, define GNL objects
  a_jk <- nest_allocation
  B <- 1*(a_jk > 0)

  # If no GNL parameters, treat as standard logit. One nest. mu=1.
  if (any(is.na(nest_allocation))) {
    K <- 1
    B <- matrix(1, ncol = 1, nrow = J)
    a_jk <- B
    mu <- rep(1,K)
  }

  # first create ownership matrices
  #own_fun_down <- function(x) {as.numeric(x == own_down)}
  #own_R <- t(sapply(own_down, own_fun_down) )
  # Since ownership is given by matrix directly:
  own_R <- own

  # then calculate foc
  shares <- share_calc(price = price, alpha = alpha, delta = delta,
                       nest_allocation = a_jk, mu = mu)
  m <- price - cost
  dd <- numDeriv::jacobian(share_calc, x = price, delta = delta, alpha = alpha,
                 nest_allocation=a_jk, mu=mu)
  omega <- (own_R * t(dd))

  # since all costs available:
  foc <- omega %*% m + shares

  if (sumFOC == FALSE) {
    return(foc)
  } else {
    out <- sum(foc^2)
    return(out)
  }

}


