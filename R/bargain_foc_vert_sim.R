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
#' @param sumFOC logical; if true, returns the sum of first-order conditions, if
#' false, returns vector of FOC's for each product
#'
#' @returns The first-order conditions
#'
#' @details This function calculate the first-order conditions from a Nash
#' bargaining model. For use in a vertical supply chain with simultaneous timing
#'
#' @examples
#' bargain_foc_vert_sim(price_w = c(1.6, 1.6, 1.6, 1.6, 1.6, 1.6),
#' own_down = paste0("R",rep(c(1,2,3),each=2)),
#' own_up = paste0("W",rep(c(1,2),3)),
#' alpha = -0.9, delta = c(0.2, 0.3, 0.9, 1.0, 0.8, 0.9),
#' cost_w = rep(.2, times = 6),
#' cost_r = rep(.1, times = 6),
#' lambda = 0.5,
#' price_r = c(2.9, 2.9, 3.0, 3.0, 3.0, 3.0))
#'
#' @export



##################################################################
# Nash Product
##################################################################
# This function allows for zero or one integrated goods

bargain_foc_vert_sim <- function(price_w,own_down,own_up,alpha,delta,cost_w,cost_r,
                                 lambda,
                                 price_r,sumFOC = TRUE){

  J <- length(price_w)

  # first create ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )

  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )

  VI_D <- as.numeric(own_down == own_up)
  k <- which(VI_D == 1)

  # VI_id identifies integrated firm. Should be singleton or empty:
  VI_id <- unique(own_down[which(own_down %in% own_up)])
  if (length(VI_id) > 1) { stop('Function cannot handle multiple integrated firms')}
  if (length(VI_id) == 0) { VI_id <- ""}

  VI_U_w <- as.numeric(VI_id == own_up) - VI_D
  VI_U_r <- as.numeric(VI_id == own_down) - VI_D

  # calculate shares and disagreement shares
  shares <- (exp(delta + alpha*price_r))/(1+sum(exp(delta + alpha*price_r)))

  shares_tilde <- vector("list",J)

  denom_tilde <- (1-diag(J)) %*% exp(delta + alpha*price_r)
  for (j in (1:J)) {
    shares_tilde[[j]] <- (exp(delta + alpha*price_r))/(1+denom_tilde[j])
    shares_tilde[[j]][j] <- 0
  }
  shares_tilde <- matrix(unlist(shares_tilde), ncol = J, byrow = FALSE)

  # define margins, depends on whether vertically integrated
  margin_up <-   (1-VI_D)*(price_w - cost_w)       + VI_D*(price_r - cost_w - cost_r)
  margin_down <- (1-VI_D)*(price_r - price_w - cost_r) + VI_D*(price_r - cost_w - cost_r)

  # specify payoffs and disagreement payoffs
  pi_w <- own_W %*% (margin_up*shares)

  pi_w_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_W_tilde <- own_W
    own_W_tilde[j,j] <- 0
    temp <- own_W_tilde %*% (margin_up*shares_tilde[,j])
    pi_w_tilde[[j]] <- temp[j]
  }

  pi_r <- own_R %*% (margin_down*shares)

  pi_r_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_R_tilde <- own_R
    own_R_tilde[j,j] <- 0
    temp <- own_R_tilde %*% (margin_down*shares_tilde[,j])
    pi_r_tilde[[j]] <- temp[j]
  }

  # RRC Effect
  RRC_effect2 <- t(shares_tilde - replicate(J,shares)) %*% (VI_U_r * margin_down)
  RRC_effect <- VI_U_w * RRC_effect2

  # Recapture Leverage Effect
  Recap_effect <- t(shares_tilde - replicate(J,shares)) %*% (VI_U_w * margin_up)
  Recap_effect <- VI_U_r * Recap_effect

  foc <- lambda*(1)*(pi_w - pi_w_tilde - RRC_effect) -
    (1-lambda)*(pi_r - pi_r_tilde - Recap_effect)*(1)

  foc[VI_D] <- 0  # set foc value for integrated goods to zero

  if (sumFOC == FALSE) {
    return(foc)
  } else {
    out <- sum(foc^2)
    return(out)
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
#' @param sumFOC logical; if true, returns the sum of first-order conditions, if
#' false, returns vector of FOC's for each product
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


bargain_foc_vert_sim_gnl <- function(p_W,own_down,own_up,alpha,delta,
                                     c_W,c_R,lambda,p_R,
                                     a_jk=NA, B=NA, mu=NA,sumFOC = TRUE){

  # If no GNL parameters, treat as standard logit. One nest. mu=1.
  J <- length(p_W)
  if (any(is.na(B))) {
    K <- 1
    B <- matrix(1, ncol = 1, nrow = J)
    a_jk <- B
    mu <- rep(1,K)
  }

  # first create ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )

  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )

  VI_D <- as.numeric(own_down == own_up)
  k <- which(VI_D == 1)

  # VI_id identifies integrated firm. Should be singleton or empty:
  VI_id <- unique(own_down[which(own_down %in% own_up)])
  if (length(VI_id) > 1) { stop('Function cannot handle multiple integrated firms')}
  if (length(VI_id) == 0) { VI_id <- ""}

  VI_U_w <- as.numeric(VI_id == own_up) - VI_D
  VI_U_r <- as.numeric(VI_id == own_down) - VI_D

  # calculate shares and disagreement shares
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
  shares_tilde <- matrix(unlist(shares_tilde), ncol = J, byrow = FALSE)

  # define margins, depends on whether vertically integrated
  margin_up <-   (1-VI_D)*(p_W - c_W)       + VI_D*(p_R - c_W - c_R)
  margin_down <- (1-VI_D)*(p_R - p_W - c_R) + VI_D*(p_R - c_W - c_R)

  # specify payoffs and disagreement payoffs
  pi_w <- own_W %*% (margin_up*shares)

  pi_w_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_W_tilde <- own_W
    own_W_tilde[j,j] <- 0
    temp <- own_W_tilde %*% (margin_up*shares_tilde[,j])
    pi_w_tilde[[j]] <- temp[j]
  }

  pi_r <- own_R %*% (margin_down*shares)

  pi_r_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_R_tilde <- own_R
    own_R_tilde[j,j] <- 0
    temp <- own_R_tilde %*% (margin_down*shares_tilde[,j])
    pi_r_tilde[[j]] <- temp[j]
  }

  # RRC Effect
  RRC_effect2 <- t(shares_tilde - replicate(J,shares)) %*% (VI_U_r * margin_down)
  RRC_effect <- VI_U_w * RRC_effect2

  # Recapture Leverage Effect
  Recap_effect <- t(shares_tilde - replicate(J,shares)) %*% (VI_U_w * margin_up)
  Recap_effect <- VI_U_r * Recap_effect

  foc <- lambda*(1)*(pi_w - pi_w_tilde - RRC_effect) -
    (1-lambda)*(pi_r - pi_r_tilde - Recap_effect)*(1)

  foc[VI_D] <- 0  # set foc value for integrated goods to zero

  if (sumFOC == FALSE) {
    return(foc)
  } else {
    out <- sum(foc^2)
    return(out)
  }
}
