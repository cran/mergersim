#' Nash Product FOCs for use in vertical model with simultaneous timing
#'
#' @param price_w Upstream or wholesale prices
#' @param own_down Ownership matrix for downstream firms
#' @param own_up Ownership matrix for upstream firms
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost_w Marginal costs for upstream firm for each product
#' @param cost_r Marginal costs for downstream firm for each product
#' @param lambda Bargaining power of the buyer/retailer
#' @param price_r Retail prices starting values
#' @param returnGFT logical; if true, returns gains from trade rather than FOC's
#'
#' @returns The first-order conditions
#'
#' @details This function calculate the first-order conditions from a Nash
#' bargaining model. For use in a vertical supply chain with simultaneous timing
#'
#' @examples
#' TO BE ADDED.
#'
#' @noRd



##################################################################
# Nash Product
##################################################################
# This function is unexported, generalized by bargain_foc_vert_sim
# This version has no vertical integration

bargain_foc_novert_sim <- function(price_w,own_down,own_up,alpha,delta,
                                   cost_w,cost_r,lambda,price_r,
                                   returnGFT = FALSE){
  J <- length(price_w)

  # construct ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )

  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )

  # calculate shares
  shares <- (exp(delta + alpha*price_r))/(1+sum(exp(delta + alpha*price_r)))

  # counterfactual shares
  shares_tilde <- vector("list",J)

  denom_tilde <- (1-diag(J)) %*% exp(delta + alpha*price_r)
  for (j in (1:J)) {
    shares_tilde[[j]] <- (exp(delta + alpha*price_r))/(1+denom_tilde[j])
    shares_tilde[[j]][j] <- 0
  }
  shares_tilde <- matrix(unlist(shares_tilde), ncol = J, byrow = FALSE)

  pi_w <- own_W %*% ((price_w - cost_w)*shares)

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
  if (returnGFT == FALSE) {
    return(sum(foc^2))
  } else {
    return(list("r_gft" = (pi_r - pi_r_tilde),"w_gft" = (pi_w-pi_w_tilde)))
  }

}




#' Nash Product FOCs for use in vertical model with simultaneous timing
#'
#' @param price_w Upstream or wholesale prices
#' @param own_down Ownership matrix for downstream firms
#' @param own_up Ownership matrix for upstream firms
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost_w Marginal costs for upstream firm for each product
#' @param cost_r Marginal costs for downstream firm for each product
#' @param lambda Bargaining power of the buyer/retailer
#' @param price_r Retail prices starting values
#' @param returnGFT logical; if true, returns gains from trade rather than FOC's
#'
#' @returns The first-order conditions
#'
#' @details This function calculate the first-order conditions from a Nash
#' bargaining model. For use in a vertical supply chain with simultaneous timing
#'
#' @examples
#' TO BE ADDED.
#'
#' @noRd


## Allowing for vertical integration to use with simultaneous
## timing model, and that can account for GNL.

## technically, this first function _novert_ should not be needed. Included only
## for completeness. Not used in the vignette.

bargain_foc_novert_sim_gnl <- function(p_W,own_down,own_up,alpha,delta,c_W,c_R,lambda,p_R,
                                       a_jk=NA, B=NA, mu=NA, returnGFT = FALSE){

  # If no GNL parameters, treat as standard logit. One nest. mu=1.
  J <- length(p_W)
  if (any(is.na(B))) {
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
  shares <- share_calc(price = p_R, alpha = alpha, delta = delta,
                       nest_allocation = a_jk, mu = mu)

  # counterfactual shares
  shares_tilde <- vector("list",J)

  for (i in (1:J)) {
    delta_cf <- delta
    delta_cf[i] <- -Inf
    cf_sharei <- share_calc(price = p_R, alpha = alpha, delta = delta_cf,
                            nest_allocation = a_jk, mu = mu)
    shares_tilde[[i]] <- cf_sharei
  }

  #denom_tilde <- (1-diag(J)) %*% exp(delta - alpha*p_R)
  #for (j in (1:J)) {
  #  shares_tilde[[j]] <- (exp(delta - alpha*p_R))/(1+denom_tilde[j])
  #  shares_tilde[[j]][j] <- 0
  #}
  shares_tilde <- matrix(unlist(shares_tilde), ncol = J, byrow = FALSE)

  pi_w <- own_W %*% ((p_W - c_W)*shares)

  pi_w_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_W_tilde <- own_W
    own_W_tilde[j,j] <- 0
    temp <- own_W_tilde %*% ((p_W - c_W)*shares_tilde[,j])
    pi_w_tilde[[j]] <- temp[j]
  }

  pi_r <- own_R %*% ((p_R - p_W - c_R)*shares)

  pi_r_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_R_tilde <- own_R
    own_R_tilde[j,j] <- 0
    temp <- own_R_tilde %*% ((p_R - p_W - c_R)*shares_tilde[,j])
    pi_r_tilde[[j]] <- temp[j]
  }

  foc <- lambda*(1)*(pi_w-pi_w_tilde) - (1-lambda)*(pi_r - pi_r_tilde)*(1)
  if (returnGFT == FALSE) {
    return(sum(foc^2))
  } else {
    return(list("r_gft" = (pi_r - pi_r_tilde),"w_gft" = (pi_w-pi_w_tilde)))
  }

}

