#' Calibrate second score auction model with GNL demand
#'
#' @param param Vector of demand parameters (alpha,mu)
#' @param price Observed prices
#' @param own Ownership matrix
#' @param shares Observed market shares
#' @param cost Marginal costs for each product
#' @param weight Weighting matrix of dimensions J-by-J
#' @param nest_allocation For generalized nested logit demand, a J-by-K matrix
#' where each element (j,k) designates the membership of good j in nest k. Rows
#' should sum to 1.
#' @param mu_constraint_matrix is a (K-by-K') matrix indicating which nesting
#' parameters are constrained to be equal to each other, where K is the
#' number of nests and K' is the number of freely varying nesting
#' parameters. mu_full = mu_constraint_matrix %*% mu_prime. Where mu_full
#' is a vector of length K of the nesting parameter value for each nest,
#' and mu_prime is a vector of length K' of parameters to be calculated.
#' It must be the case that K is greater than K'.
#'
#' @returns Difference between model predicted and observed values of
#' prices, shares, and diversions.
#'
#' @details This function calibrates a second score auction model with generalized
#' nested logit (GNL) demand
#'
#' @examples
#' nest1 <- matrix( c(1, 0, 0, 0, 1, 1), ncol = 2, nrow = 3)
#'
#' ssa_calibrate_gnl(param = c(-0.9, 1, 1),
#' own = diag(3),
#' price = c(.05, .34, .33),
#' shares = c( 0.31, 0.27, 0.25),
#' cost = c(.05,.31,.30),
#' weight = diag(c(1,1,1)),
#' nest_allocation = nest1,
#' mu_constraint_matrix = NA)
#'
#' @export


##################################################################
# Second score auction pricing equations for calibration with GNL demand
##################################################################


ssa_calibrate_gnl <- function(param, own, price,
                              shares, cost, weight,
                              nest_allocation,
                              mu_constraint_matrix = NA){

  # If GNL, define GNL objects
  a_jk <- nest_allocation
  B <- 1*(a_jk > 0)
  K_val <- dim(a_jk)[2]

  alpha <- param[1]

  #### checks on mu_constraint matrix ####

  mu_prime <- param[2:length(param)]
  K_prime <- length(mu_prime)

  ## Intended nesting constraints partly
  ## implied by length of param. Check to make sure consistent, and if not
  ## then throw error that more information needs to be supplied.

  ## give error if K_prime > K_val. This maybe should be a stop()
  if (K_prime > K_val) {warning("K' should not be greater than K")}

  ## if no matrix provided, but K_prime implied by length of param is 1,
  ## then we can assume just one nesting parameter
  if (anyNA(mu_constraint_matrix) & K_prime == 1 ) {
    mu_constraint_matrix <- matrix(1, nrow = K_val, ncol = 1)
  }

  ## if no matrix provided, but K_prime implied by length of param is K,
  ## then we can assume full flexibility intended
  if (anyNA(mu_constraint_matrix) & K_prime == K_val ) {
    mu_constraint_matrix <- diag(K_val)
  }

  ## if still no mcm determined, give error that more information is need
  if (anyNA(mu_constraint_matrix) ) {
    warning("Please provide more information on nesting parameter calibration
            in mu_constraint_matrix")
  }


  mu <- mu_constraint_matrix %*% mu_prime

  #### calculate delta
  delta0 <- log(shares) - log(1-sum(shares)) - alpha*price

  find_d <- rootSolve::multiroot(f = match_share, start = delta0,
                      price=cost, alpha=alpha, nest_allocation=a_jk,
                      mu=mu,
                      shares_obs = shares)

  delta <- find_d$root

  #### calculated expected values of maximum
  J <- length(price)

  E_z <- share_calc(price = cost, delta = delta, alpha = alpha,
                    nest_allocation = nest_allocation, mu = mu,
                    returnLogsum = TRUE)

  E_z_prime <- rep(0,J)

  for (i in 1:J) {
    delta_cf <- delta

    delta_cf[own[i,] == 1] <- -Inf

    logsum_i <- share_calc(price = cost, delta = delta_cf,
                           alpha = alpha,
                           nest_allocation = nest_allocation, mu = mu,
                           returnLogsum = TRUE)
    E_z_prime[i] <- logsum_i
  }


  p_hat <- cost - 1/(alpha*own%*%shares) * (E_z - E_z_prime)

  FOC <- p_hat - price
  FOC <- stats::na.omit(FOC)

  objfxn <- c(FOC) %*% diag(weight) %*% c(FOC)

  return(objfxn)
}

