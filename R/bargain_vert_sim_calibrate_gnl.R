#' Calibrate upstream Nash bargain in vertical model with simultaneous timing
#'
#' @param param vector of parameters to calibrate, will be either bargaining
#'  parameter or demand parameters, depending on what other information is
#'  supplied
#' @param lambda Bargaining power of the buyer/retailer
#' @param cost_w Marginal costs for upstream firm for each product
#' @param price_w Upstream or wholesale prices
#' @param own_down Ownership matrix for downstream firms
#' @param own_up Ownership matrix for upstream firms
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost_r Marginal costs for downstream firm for each product
#' @param price_r Retail prices
#' @param nest_allocation For generalized nested logit demand, a J-by-K matrix
#' where each element (j,k) designates the membership of good j in nest k. Rows
#' should sum to 1.
#' @param mu Nesting parameters for each nest
#' @param sumFOC logical; whether to return the sum of squares of
#' the first-order conditions. Defaults to FALSE, in which case it returns each
#' product first-order condition as a vector.
#'
#' @returns The value of objective function
#'
#' @details This function can be used to calibrate the bargaining parameter
#' in a Nash bargain which is the upstream market of a vertical supply chain.
#' Assumes logit demand.
#'
#' @examples
#' bargain_vert_sim_calibrate_gnl(param = 0.4, lambda=NA,
#' cost_w = rep(.2, times = 6),
#' own_down = paste0("R",rep(c(1,2,3),each=2)),
#' own_up = paste0("W",rep(c(1,2),3)),
#' alpha = -0.9,
#' delta = c(0.2, 0.3, 0.9, 1.0, 0.8, 0.9),
#' cost_r = rep(.1, times = 6),
#' price_r = c(2.9, 2.9, 3.0, 3.0, 3.0, 3.0),
#' price_w = c(1.6, 1.6, 1.6, 1.6, 1.6, 1.6),
#' nest_allocation=NA, mu=NA, sumFOC = TRUE)
#'
#' @export



##################################################################
# Nash Product
##################################################################

#### function for calibrating upstream in simultaneous model


bargain_vert_sim_calibrate_gnl <- function(param,lambda=NA,cost_w=NA,
                                           own_down,own_up,alpha,delta,
                                           cost_r,price_r,price_w,
                                           nest_allocation=NA, mu=NA,
                                           sumFOC = FALSE){

  # If cost_w is NA, assume we are calibrating it. If lambda is NA, assume we are
  # calibrating it. If both are NA, give an error. If neither are NA, give error.
  if (any(is.na(lambda)) & any(is.na(cost_w))) {
    message("Must provide either lambda or cost_w")
  }
  if (!any(is.na(lambda)) & !any(is.na(cost_w))) {
    message("Provide either lambda or cost_w, but not both")
  }
  if (any(is.na(lambda))) {
    lambda <- param
  }
  if (any(is.na(cost_w))) {
    cost_w <- param
  }

  # If no GNL parameters, treat as standard logit. One nest. mu=1.
  J <- length(price_w)

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

  # construct ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )

  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )

  # calculate shares
  shares <- share_calc(price = price_r, alpha = alpha, delta = delta,
                       nest_allocation = a_jk, mu = mu)

  # counterfactual shares
  shares_tilde <- vector("list",J)

  for (i in (1:J)) {
    delta_cf <- delta
    delta_cf[i] <- -Inf
    cf_sharei <- share_calc(price = price_r, alpha = alpha, delta = delta_cf,
                            nest_allocation = a_jk, mu = mu)
    shares_tilde[[i]] <- cf_sharei
  }

  shares_tilde <- matrix(unlist(shares_tilde), ncol = J, byrow = FALSE)

  pi_w <- own_W %*% ((price_w - cost_w)*shares)  # this won't work if some cost_w missing

  pi_w_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_W_tilde <- own_W
    own_W_tilde[j,j] <- 0
    temp <- own_W_tilde %*% ((price_w - cost_w)*shares_tilde[,j])
    pi_w_tilde[[j]] <- temp[j]
  }

  pi_r <- own_R %*% ((price_r - price_w - cost_r)*shares)

  pi_r_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_R_tilde <- own_R
    own_R_tilde[j,j] <- 0
    temp <- own_R_tilde %*% ((price_r - price_w - cost_r)*shares_tilde[,j])
    pi_r_tilde[[j]] <- temp[j]
  }

  foc <- lambda*(1)*(pi_w-pi_w_tilde) - (1-lambda)*(pi_r - pi_r_tilde)*(1)
  if (sumFOC == TRUE) {
    return(sum(foc^2))
  } else {
    return(foc)
  }

}
