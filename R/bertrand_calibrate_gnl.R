#' Calibrate Bertrand model with GNL demand
#'
#' @param param Vector of demand parameters (alpha,mu)
#' @param price Observed prices
#' @param own Ownership matrix
#' @param shares Observed market shares
#' @param cost Marginal costs for each product
#' @param weight Weighting vector of length equal to number of margins provided;
#' if diversions are provided, these weights are relative to weight on matching
#' diversions.
#' @param nest_allocation For generalized nested logit demand, a J-by-K matrix
#' where each element (j,k) designates the membership of good j in nest k. Rows
#' should sum to 1.
#' @param div_matrix A matrix of observed diversions from product in row j to
#' product in column k.
#' @param mu_constraint_matrix is a (K-by-K') matrix indicating which nesting
#' parameters are constrained to be equal to each other, where K is the
#' number of nests and K' is the number of freely varying nesting
#' parameters. mu_full = mu_constraint_matrix %*% mu_prime. Where mu_full
#' is a vector of length K of the nesting parameter value for each nest,
#' and mu_prime is a vector of length K' of parameters to be calculated.
#' It must be the case that K is greater than K'.
#' @param div_calc_marginal is a logical if function should match to marginal
#' diversions (if TRUE) or second choice diversions (if FALSE). Default
#' to TRUE.
#' @param returnOutcomes logical; should equilibrium objects be returned (mean
#' value parameter, prices, shares, costs) as a list.
#'
#' @returns Difference between model predicted and observed values of
#' prices, shares, and diversions.
#'
#' @details This function calibrates a Bertrand model with generalized nested
#' logit (GNL) demand
#'
#' @examples
#' nest1 <- matrix( c(1, 0, 0, 0, 1, 1), ncol = 2, nrow = 3)
#' divmat <- matrix( c(0, .4, .4, .4, 0, .4, .4, .4, 0), ncol = 3, nrow = 3)
#'
#' bertrand_calibrate_gnl(param = c(-0.9, 1, 1),
#'                        own = diag(3),
#'                        price = c(.05, .34, .33),
#'                        shares = c( 0.31, 0.27, 0.25),
#'                        cost = c(.05,.31,.30),
#'                        weight = c(1,1,1),
#'                        nest_allocation = nest1, div_matrix = divmat)
#'
#' @export


##################################################################
# Bertrand model first-order conditions for calibration with GNL demand
##################################################################


bertrand_calibrate_gnl <- function(param,own,price,shares,cost,
                                      weight = NA,
                                      nest_allocation, div_matrix = NA,
                                      mu_constraint_matrix = NA,
                                      div_calc_marginal = TRUE,
                                      returnOutcomes = FALSE){

  # checks for length of weighting vector
  num_c <- sum(!is.na(cost))

  if (anyNA(weight)) {
    weight <- rep(1, times = num_c)
  } else if (length(weight) != num_c) {
    stop("weight must have length equal to the number of margins (", num_c, ")")
  }

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
  ## then assume just one nesting parameter
  if (anyNA(mu_constraint_matrix) & K_prime == 1 ) {
    mu_constraint_matrix <- matrix(1, nrow = K_val, ncol = 1)
  }

  ## if no matrix provided, but K_prime implied by length of param is K,
  ## then assume full flexibility intended
  if (anyNA(mu_constraint_matrix) & K_prime == K_val ) {
    mu_constraint_matrix <- diag(K_val)
  }

  ## if still no mcm determined, give error that more information is need
  if (anyNA(mu_constraint_matrix) ) {
    warning("Please provide more information on nesting parameter calibration
            in mu_constraint_matrix")
  }


  mu <- mu_constraint_matrix %*% mu_prime


  # Find delta that matches shares
  delta0 <- log(shares) - log(1-sum(shares)) - alpha*price

  find_d <- rootSolve::multiroot(f = match_share, start = delta0,
                                 price=price,alpha=alpha,nest_allocation=a_jk,mu=mu,
                                 shares_obs = shares)

  delta <- find_d$root

  # calculate matrix of derivatives

  dd <- numDeriv::jacobian(share_calc, x = price,
                           delta = delta, alpha = alpha,
                           nest_allocation = a_jk,
                           mu = mu)
  omega <- own * t(dd)

  # Back out implied cost from FOC at observed prices: Omega*(p-c) + s = 0
  # => c_implied = p + solve(Omega) %*% s
  c_implied <- as.numeric(price + solve(omega) %*% shares)

  # Only compute residuals where cost is observed
  cdiff <- cost - c_implied
  keep <- !is.na(cdiff)

  # Compute diversions
  diversions_m <- diversion_calc(price=price,alpha=alpha,delta=delta,
                                 nest_allocation=a_jk,mu=mu,
                                 marginal = div_calc_marginal)

  div_diff <- (as.numeric(div_matrix - diversions_m)^2)*100

  if (returnOutcomes == FALSE) {
    objfxn <- sum(weight * cdiff[keep]^2) + sum(div_diff, na.rm = TRUE)
    return(objfxn)
  }

  if (returnOutcomes == TRUE) {
    return(list("cost_cal" = c_implied,
                "delta_cal" = delta) )
  }

}




#' Calibrate Bertrand model with GNL demand
#'
#' @param param Vector of nesting parameters (mu)
#' @param alpha Value of price coefficient, to be held fixed during calibration
#' @param price Observed prices
#' @param own Ownership matrix
#' @param share Observed market shares
#' @param cost Marginal costs for each product
#' @param weight Vector of weights given to prices, shares, diversions, respectively
#' @param nest_allocation For generalized nested logit demand, a J-by-K matrix
#' where each element (j,k) designates the membership of good j in nest k. Rows
#' should sum to 1.
#' @param div_matrix A matrix of observed diversions from product in row j to
#' product in column k.
#' @param mu_constraint_matrix is a (K-by-K') matrix indicating which nesting
#' parameters are constrained to be equal to each other, where K is the
#' number of nests and K' is the number of freely varying nesting
#' parameters. mu_full = mu_constraint_matrix %*% mu_prime. Where mu_full
#' is a vector of length K of the nesting parameter value for each nest,
#' and mu_prime is a vector of length K' of parameters to be calculated.
#' It must be the case that K is greater than K'.
#' @param div_calc_marginal is a logical if function should match to marginal
#' diversions (if TRUE) or second choice diversions (if FALSE). Default
#' to TRUE.
#' @param optimizer Which optimization routine should be used to find
#' equilibrium prices, either BBoptim or multiroot
#' @param returnOutcomes logical; should equilibrium objects be returned (mean
#' value parameter, prices, shares, costs) as a list.
#'
#' @returns Difference between model predicted and observed values of
#' prices, shares, and diversions.
#'
#' @details This function calibrates a Bertrand model with generalized nested
#' logit (GNL) demand, using only first-order conditions that are available,
#' i.e. first-order conditions for products that have non-missing costs.
#'
#' @examples
#' TO BE ADDED.
#' @noRd


## useOldWeight is a legacy option in case want to use old weighting in the
## objective function

##################################################################
# Bertrand model first-order conditions for calibration with GNL demand
##################################################################

## NOTE: These functions are not currently exported in mergersim


bertrand_calibrate_mu <- function(param,own,alpha,price,shares,cost,
                                  weight,nest_allocation,div_matrix,
                                  mu_constraint_matrix = NA,
                                  div_calc_marginal = TRUE,
                                  optimizer="BBoptim",
                                  useOldWeight = FALSE,
                                  returnOutcomes = FALSE){

  ## this function takes alpha as given. Finds best nesting parameters, mu.

  # If GNL, define GNL objects
  a_jk <- nest_allocation
  B <- 1*(a_jk > 0)
  K_val <- dim(a_jk)[2]


  #### checks on mu_constraint matrix ####

  mu_prime <- param
  K_prime <- length(mu_prime)

  ## intended nesting constraints partly
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

  J <- length(price)

  delta0 <- log(shares) - log(1-sum(shares)) - alpha*price

  find_d <- rootSolve::multiroot(f = match_share, start = delta0,
                      price=price,alpha=alpha,nest_allocation=a_jk,mu=mu,
                      shares_obs = shares)

  delta <- find_d$root

  # sometimes, multiroot inexplicably fails inside of optimization. Set to delta0.
  # in those cases.

  #useOld <- TRUE  # if want to use old version of function, which had no correction
  useOld <- FALSE # if want to use new version, with delta NA correction

  if (useOld == FALSE) {
    if (anyNA(delta)) {
      delta <- delta0
      warning("Multiroot failed to find mean values that matched shares.")
    }
  }

  x0 <- price


  ## Given observed prices, and current guess of demand parameter values,
  ## back out costs consistent with FOCs.
  x06 <- price * 0.5

  out_cost <- stats::optim(f = bertrand_foc_c, par = x06,
                    price = price, own = own, alpha = alpha,
                    delta = delta,
                    nest_allocation=a_jk, mu=mu, sumFOC = TRUE,
                    control = list(maxit = 1000) )

  cost_cal <- out_cost$par


  ## BBoptim or multiroot. If there are missing costs, use BBoptim
  if (optimizer == "BBoptim") {
    out1 <- BB::BBoptim(f = bertrand_foc, par = x0,
                    own = own, alpha= alpha,
                    delta = delta, cost = cost_cal,
                    nest_allocation = a_jk, mu = mu,
                    sumFOC = TRUE, control = list(trace=FALSE))

    p_model <- out1$par
  }
  if (optimizer == "multiroot") {
    out1 <- rootSolve::multiroot(f = bertrand_foc, start = x0,
                      own = own, alpha= alpha,
                      delta = delta, cost = cost_cal,
                      nest_allocation = a_jk, mu = mu)

    p_model <- out1$root
  }

  share_m <- share_calc(price=price,delta=delta,alpha=alpha,nest_allocation=a_jk,
                        mu=mu)
  diversions_m <- diversion_calc(price=price,alpha=alpha,delta=delta,
                                 nest_allocation=a_jk,mu=mu,
                                 marginal = div_calc_marginal)

  if (useOldWeight == TRUE) {
    pdiff <- price - p_model
    sdiff <- shares - share_m
    div_diff <- sum((as.numeric(div_matrix - diversions_m)^2), na.rm = TRUE)
    cost_diff <- sum((cost - cost_cal)^2, na.rm = TRUE)

    objfxn <- c(pdiff,sdiff,div_diff,cost_diff) %*% weight %*%
      c(pdiff,sdiff,div_diff,cost_diff)
  }


  if (useOldWeight == FALSE) {
    pdiff <- ((price - p_model)^2) * weight[1]
    sdiff <- ((shares - share_m)^2) * 1000 * weight[2]
    div_diff <- ((as.numeric(div_matrix - diversions_m)^2)*100 * weight[3])
    cost_diff <- ((cost - cost_cal)^2) * weight[4]

    objfxn <- sum(pdiff) + sum(sdiff) + sum(div_diff, na.rm = TRUE) +
      sum(cost_diff, na.rm = TRUE)
  }


  if (returnOutcomes == FALSE) {
    return(objfxn)
  }
  if (returnOutcomes == TRUE) {
    return(list("FOCs" = c(pdiff,sdiff,div_diff,cost_diff),
                "delta_cal" = delta,
                "p_model" = p_model,
                "share_m" = share_m,
                "cost_cal" = cost_cal) )
  }
}


#' Calibrate Bertrand model with GNL demand
#'
#' @param param Vector of nesting parameters (alpha)
#' @param alpha Value of price coefficient, to be held fixed during calibration
#' @param price Observed prices
#' @param own Ownership matrix
#' @param share Observed market shares
#' @param cost Marginal costs for each product
#' @param weight Vector of weights given to prices, shares, diversions, respectively
#' @param nest_allocation For generalized nested logit demand, a J-by-K matrix
#' where each element (j,k) designates the membership of good j in nest k. Rows
#' should sum to 1.
#' @param div_matrix A matrix of observed diversions from product in row j to
#' product in column k.
#' @param mu_constraint_matrix is a (K-by-K') matrix indicating which nesting
#' parameters are constrained to be equal to each other, where K is the
#' number of nests and K' is the number of freely varying nesting
#' parameters. mu_full = mu_constraint_matrix %*% mu_prime. Where mu_full
#' is a vector of length K of the nesting parameter value for each nest,
#' and mu_prime is a vector of length K' of parameters to be calculated.
#' It must be the case that K is greater than K'.
#' @param div_calc_marginal is a logical if function should match to marginal
#' diversions (if TRUE) or second choice diversions (if FALSE). Default
#' to TRUE.
#' @param optimizer Which optimization routine should be used to find
#' equilibrium prices, either BBoptim or multiroot
#' @param returnOutcomes logical; should equilibrium objects be returned (mean
#' value parameter, prices, shares, costs) as a list.
#'
#' @returns Difference between model predicted and observed values of
#' prices, shares, and diversions.
#'
#' @details This function calibrates a Bertrand model with generalized nested
#' logit (GNL) demand, using only first-order conditions that are available,
#' i.e. first-order conditions for products that have non-missing costs.
#'
#' @examples
#' TO BE ADDED.
#' @noRd

bertrand_calibrate_alpha <- function(param,own,price,shares,cost,
                                     weight,nest_allocation,div_matrix,
                                     mu_start,
                                     mu_constraint_matrix = NA,
                                     mu_lower = NA, mu_upper = NA,
                                     div_calc_marginal = TRUE,
                                     optimizer="BBoptim",
                                     useOldWeight = FALSE,
                                     returnOutcomes = FALSE){


  if (!(length(param)==1)) {warning("param should be length one")}

  J <- length(price)
  alpha <- param[1]

  mu_length <- length(mu_start)
  if (is.na(mu_lower)) {
    mu_lower <- rep(0, times = mu_length)
  }
  if (is.na(mu_upper)) {
    mu_upper <- rep(1, times = mu_length)
  }


  out_val <- stats::optim(f = bertrand_calibrate_mu, par = mu_start,
                   alpha = alpha,
                   own = own, price = price,
                   shares = shares, cost  = cost,
                   weight = weight, nest_allocation = nest_allocation,
                   div_matrix = div_matrix,
                   mu_constraint_matrix = mu_constraint_matrix,
                   div_calc_marginal = div_calc_marginal, optimizer=optimizer,
                   lower = mu_lower,
                   upper = mu_upper )

  mu <- out_val$par

  objfxn <- out_val$value

  if (returnOutcomes == FALSE) {
    return(objfxn)
  }
  if (returnOutcomes == TRUE) {
    return(list("mu_cal" = mu) )
  }
}



