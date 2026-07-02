#' Nash Product for use in vertical model with sequential timing
#'
#' @param w_start Initial wholesale price for product for which to maximize Nash Product
#' @param product_max Index of product for which to maximize Nash Product
#' @param price_w Upstream or wholesale prices
#' @param own_down Ownership matrix for downstream firms
#' @param own_up Ownership matrix for upstream firms
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost_w Marginal costs for upstream firm for each product
#' @param cost_r Marginal costs for downstream firm for each product
#' @param lambda Bargaining power of the buyer/retailer
#' @param p_R0 Retail prices starting values
#' @param sigma Contract type; value between 0 and 1 where 0 is linear price,
#'  1 is two-part tariff
#' @param showAll logical; if true, returns gains for trade for every product
#' @param maxJointProfits Experimental option; if true, will maximize joint
#' profits rather than gains from trade
#'
#' @returns The first-order conditions
#'
#' @details This function calculate the Nash Product from a Nash bargaining
#' model
#'
#' @examples
#' bargain_NP_vert_seq(w_start = 1.5, product_max = 1,
#' price_w = c(1.6, 1.6, 1.6, 1.6, 1.6, 1.6),
#' own_down = paste0("R",rep(c(1,2,3),each=2)),
#' own_up = paste0("W",rep(c(1,2),3)),
#' alpha = -0.9, delta = c(0.2, 0.3, 0.9, 1.0, 0.8, 0.9),
#' cost_w = rep(.2, times = 6),
#' cost_r = rep(.1, times = 6),
#' lambda = 0.5, sigma = 0,
#' p_R0 = c(2.9, 2.9, 3.0, 3.0, 3.0, 3.0))
#'
#' @export



##################################################################
# Nash Product
##################################################################
#### downstream prices to adjust in disagreement scenarios
## maxJointProfits = TRUE is an experimental option

bargain_NP_vert_seq <- function(w_start,product_max,price_w,own_down,own_up,alpha,
                                delta,
                                cost_w,cost_r,lambda,p_R0,sigma,showAll = FALSE,
                                maxJointProfits = FALSE){


  price_w[product_max] <- w_start
  J <- length(price_w)

  # construct ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )

  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )

  # indicator of which goods are integrated
  VI_D <- as.numeric(own_down == own_up)
  VI_D_idx <- which(own_down == own_up)

  # VI_U_w = all non-integrated goods sold by upstream divisions of integrated firms. Calculate by finding all upstream goods that are owned by an integrated firm, and subtract off the integrated goods.
  # VI_U_r = all non-integrated goods sold by downstream divisions of integrated firms. Calculate by finding all downstream goods that are sold by an integrated firm, and subtract off the integrated goods.
  VI_U_w <- as.numeric(own_up %in% own_down) - VI_D
  VI_U_r <- as.numeric(own_down %in% own_up) - VI_D
  # matrix form
  own_R_up <- t(sapply(own_down,
                       function(i) as.numeric(own_up %in% i) ))
  own_W_down <- t(sapply(own_up,
                         function(i) as.numeric(own_down %in% i) ))
  # remove integrated goods from partner profits
  own_R_up[,VI_D_idx] <- 0
  own_W_down[,VI_D_idx] <- 0

  # Given the inputs, calculate optimal retail prices
  outtest <- BB::BBoptim(f = bertrand_foc_vert, par = p_R0,
                     own_down = own_down, own_up = own_up,
                     alpha= alpha,
                     delta = delta, cost_r = cost_r,
                     price_w = price_w, cost_w = cost_w, sumFOC = TRUE,
                     control = list(trace=FALSE),
                     quiet = TRUE)

  price_r <- outtest$par

  # calculate shares
  shares <- (exp(delta + alpha*price_r))/(1+sum(exp(delta + alpha*price_r)))

  # counterfactual shares
  shares_tilde <- matrix(data = 0, nrow = J, ncol = J)
  p_R_tilde <- matrix(data = 0, nrow = J, ncol = J)

  for (j in (1:J)) {
    delta_tilde <- delta
    delta_tilde[j] <- -Inf

    outtest_tilde <- BB::BBoptim(f = bertrand_foc_vert, par = p_R0,
                             own_down = own_down, own_up = own_up,
                             alpha= alpha,
                             delta = delta_tilde, cost_r = cost_r,
                             price_w = price_w, cost_w = cost_w, sumFOC = TRUE,
                             control = list(trace=FALSE),
                             quiet = TRUE)

    p_R_tilde[,j] <- outtest_tilde$par

    shares_tilde[,j] <- share_calc_logit(price=p_R_tilde[,j], delta = delta_tilde,
                                         alpha = alpha)
  }


  # define margins, depends on whether vertically integrated
  margin_up <-   (1-VI_D)*(price_w - cost_w)       + VI_D*(price_r - cost_w - cost_r)
  margin_down <- (1-VI_D)*(price_r - price_w - cost_r) + VI_D*(price_r - cost_w - cost_r)

  # specify payoffs and disagreement payoffs
  margin_up_tilde <- matrix(data = 0, nrow = J, ncol = J)
  margin_down_tilde <- matrix(data = 0, nrow = J, ncol = J)
  for (j in (1:J)) {
    margin_up_tilde[,j] <-   (1-VI_D)*(price_w - cost_w)       + VI_D*(p_R_tilde[,j] - cost_w - cost_r)
    margin_down_tilde[,j] <- (1-VI_D)*(p_R_tilde[,j] - price_w - cost_r) + VI_D*(p_R_tilde[,j] - cost_w - cost_r)
  }

  pi_w <- own_W %*% (margin_up*shares) + own_W_down %*% (margin_down*shares)

  pi_w_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_W_tilde <- own_W
    own_W_tilde[j,j] <- 0
    temp <- own_W_tilde %*% (margin_up_tilde[,j]*shares_tilde[,j]) +
      own_R %*% (VI_U_r*margin_down_tilde[,j]*shares_tilde[,j])
    pi_w_tilde[[j]] <- temp[j]
  }

  pi_r <- own_R %*% (margin_down*shares) + own_R_up %*% (margin_up*shares)

  pi_r_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_R_tilde <- own_R
    own_R_tilde[j,j] <- 0
    temp <- own_R_tilde %*% (margin_down_tilde[,j]*shares_tilde[,j]) +
      own_W %*% (VI_U_w*margin_up_tilde[,j]*shares_tilde[,j])
    pi_r_tilde[[j]] <- temp[j]
  }


  NP <- (pi_r - pi_r_tilde)^lambda * (pi_w-pi_w_tilde)^(1-lambda)
  L <- lambda^lambda * (1-lambda)^(1-lambda)
  NP_tpt <- L * (pi_r - pi_r_tilde + pi_w - pi_w_tilde)
  # could be small issue here where NP returns NaN and then tpt fails too

  # Experimental option that maximizes joint profits
  if (maxJointProfits == TRUE) {
    #NP_tpt <- L * (pi_r + pi_w)
    NP_tpt <- (pi_r + pi_w)
  }

  if (showAll == TRUE) {
    out <- (1-sigma) * NP + sigma * NP_tpt
    F_j <- (1-lambda) * (pi_r - pi_r_tilde) - lambda * (pi_w - pi_w_tilde)
    return(list("foc" = out,
                "r_gft" = (pi_r - pi_r_tilde),
                "w_gft" = (pi_w - pi_w_tilde),
                "F_j" = F_j) )
  } else {
    out <- -(1-sigma) * NP - sigma * NP_tpt  # BBoptim minimizes fn
    out <- out[product_max]
    return(out)
  }

}




#' Nash Product for use in vertical model with sequential timing
#'
#' @param w_start Initial wholesale price for product for which to maximize Nash Product
#' @param product_max Index of product for which to maximize Nash Product
#' @param price_w Upstream or wholesale prices
#' @param own_down Ownership matrix for downstream firms
#' @param own_up Ownership matrix for upstream firms
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost_w Marginal costs for upstream firm for each product
#' @param cost_r Marginal costs for downstream firm for each product
#' @param lambda Bargaining power of the buyer/retailer
#' @param p_R0 Retail prices starting values
#' @param sigma Contract type; value between 0 and 1 where 0 is linear price,
#'  1 is two-part tariff
#' @param nest_allocation For generalized nested logit demand, a J-by-K matrix
#' where each element (j,k) designates the membership of good j in nest k. Rows
#' should sum to 1.
#' @param mu Nesting parameters for each nest
#' @param showAll logical; if true, returns gains for trade for every product
#'
#' @returns The first-order conditions
#'
#' @details This function calculate the Nash Product from a Nash bargaining
#' model. This version allows for nested logit with overlapping nests.
#'
#' @examples
#'
#' own_down <- paste0("R",rep(c(1,2,3),each=2))
#' own_up <- paste0("W",rep(c(1,2),3))
#'
#' B1 <- 1 * matrix( c(own_down == "R1",
#'                     own_down == "R2",
#'                     own_down == "R3",
#'                     own_up == "W1",
#'                     own_up == "W2"),
#'                   ncol = 5, nrow = 6)
#' a1 <- B1 * 0.5
#' mu1 <- rep(1.0,5)
#'
#' bargain_NP_vert_seq_gnl(w_start = 1.5, product_max = 1,
#'                         price_w = c(1.6, 1.6, 1.6, 1.6, 1.6, 1.6),
#'                         own_down = own_down,
#'                         own_up = own_up,
#'                         alpha = -0.9, delta = c(0.2, 0.3, 0.9, 1.0, 0.8, 0.9),
#'                         cost_w = rep(.2, times = 6),
#'                         cost_r = rep(.1, times = 6),
#'                         lambda = 0.5, sigma = 0,
#'                         p_R0 = c(2.9, 2.9, 3.0, 3.0, 3.0, 3.0),
#'                         nest_allocation = a1, mu = mu1)
#'
#'
#' @export


bargain_NP_vert_seq_gnl <- function(w_start,product_max,price_w,own_down,own_up,
                                    alpha,delta,
                                    cost_w,cost_r,lambda,p_R0,sigma,
                                    nest_allocation,mu,showAll = FALSE){

  price_w[product_max] <- w_start
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

  # indicator of which goods are integrated
  VI_D <- as.numeric(own_down == own_up)
  VI_D_idx <- which(own_down == own_up)

  # VI_U_w = all non-integrated goods sold by upstream divisions of integrated firms. Calculate by finding all upstream goods that are owned by an integrated firm, and subtract off the integrated goods.
  # VI_U_r = all non-integrated goods sold by downstream divisions of integrated firms. Calculate by finding all downstream goods that are sold by an integrated firm, and subtract off the integrated goods.
  VI_U_w <- as.numeric(own_up %in% own_down) - VI_D
  VI_U_r <- as.numeric(own_down %in% own_up) - VI_D
  # matrix form
  own_R_up <- t(sapply(own_down,
                       function(i) as.numeric(own_up %in% i) ))
  own_W_down <- t(sapply(own_up,
                         function(i) as.numeric(own_down %in% i) ))
  # remove integrated goods from partner profits
  own_R_up[,VI_D_idx] <- 0
  own_W_down[,VI_D_idx] <- 0

  # Given the inputs, calculate optimal retail prices
  outtest <- BB::BBoptim(f = bertrand_foc_vert_gnl, par = p_R0,
                     own_down = own_down, own_up = own_up,
                     alpha= alpha,
                     delta = delta, cost_r = cost_r,
                     price_w = price_w, cost_w = cost_w,
                     nest_allocation=a_jk, mu=mu, sumFOC = TRUE,
                     control = list(trace=FALSE),
                     quiet = TRUE)

  price_r <- outtest$par

  # calculate shares
  shares <- share_calc(price=price_r,delta,alpha,nest_allocation=a_jk,mu)

  # counterfactual shares
  shares_tilde <- matrix(data = 0, nrow = J, ncol = J)
  p_R_tilde <- matrix(data = 0, nrow = J, ncol = J)

  for (j in (1:J)) {
    delta_tilde <- delta
    delta_tilde[j] <- -Inf

    outtest_tilde <- BB::BBoptim(f = bertrand_foc_vert, par = p_R0,
                             own_down = own_down, own_up = own_up,
                             alpha= alpha,
                             delta = delta_tilde, cost_r = cost_r,
                             price_w = price_w, cost_w = cost_w, sumFOC = TRUE,
                             control = list(trace=FALSE),
                             quiet = TRUE)

    p_R_tilde[,j] <- outtest_tilde$par

    shares_tilde[,j] <- share_calc(price=p_R_tilde[,j], delta = delta_tilde,
                                   alpha = alpha, nest_allocation=a_jk, mu=mu)
  }


  # define margins, depends on whether vertically integrated
  margin_up <-   (1-VI_D)*(price_w - cost_w)       + VI_D*(price_r - cost_w - cost_r)
  margin_down <- (1-VI_D)*(price_r - price_w - cost_r) + VI_D*(price_r - cost_w - cost_r)

  # specify payoffs and disagreement payoffs
  margin_up_tilde <- matrix(data = 0, nrow = J, ncol = J)
  margin_down_tilde <- matrix(data = 0, nrow = J, ncol = J)
  for (j in (1:J)) {
    margin_up_tilde[,j] <-   (1-VI_D)*(price_w - cost_w)       + VI_D*(p_R_tilde[,j] - cost_w - cost_r)
    margin_down_tilde[,j] <- (1-VI_D)*(p_R_tilde[,j] - price_w - cost_r) + VI_D*(p_R_tilde[,j] - cost_w - cost_r)
  }

  pi_w <- own_W %*% (margin_up*shares) + own_W_down %*% (margin_down*shares)

  pi_w_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_W_tilde <- own_W
    own_W_tilde[j,j] <- 0
    temp <- own_W_tilde %*% (margin_up_tilde[,j]*shares_tilde[,j]) +
      own_R %*% (VI_U_r*margin_down_tilde[,j]*shares_tilde[,j])
    pi_w_tilde[[j]] <- temp[j]
  }

  pi_r <- own_R %*% (margin_down*shares) + own_R_up %*% (margin_up*shares)

  pi_r_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_R_tilde <- own_R
    own_R_tilde[j,j] <- 0
    temp <- own_R_tilde %*% (margin_down_tilde[,j]*shares_tilde[,j]) +
      own_W %*% (VI_U_w*margin_up_tilde[,j]*shares_tilde[,j])
    pi_r_tilde[[j]] <- temp[j]
  }


  NP <- (pi_r - pi_r_tilde)^lambda * (pi_w-pi_w_tilde)^(1-lambda)
  L <- lambda^lambda * (1-lambda)^(1-lambda)
  NP_tpt <- L * (pi_r - pi_r_tilde + pi_w - pi_w_tilde)

  if (showAll == TRUE) {
    out <- (1-sigma) * NP + sigma * NP_tpt
    return(out)
  } else {
    out <- -(1-sigma) * NP - sigma * NP_tpt  # BBoptim minimizes fn
    out <- out[product_max]
    return(out)
  }

}




#' Nash Product for use in vertical model with sequential timing
#'
#' @noRd

## Vertical function where p^R cannot adjust in disagreement
## Unused as that is less appealing assumption for sequential timing

bargain_NP_vert_seq2 <- function(w_start,product_max,p_W,own_down,own_up,alpha,delta,
                                 cost_w,cost_r,lambda,p_R0,sigma,showAll = FALSE){

  p_W[product_max] <- w_start
  J <- length(p_W)

  # construct ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )

  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )

  # indicator of which goods are integrated
  VI_D <- as.numeric(own_down == own_up)
  VI_D_idx <- which(own_down == own_up)

  # VI_U_w = all non-integrated goods sold by upstream divisions of integrated firms. Calculate by finding all upstream goods that are owned by an integrated firm, and subtract off the integrated goods.
  # VI_U_r = all non-integrated goods sold by downstream divisions of integrated firms. Calculate by finding all downstream goods that are sold by an integrated firm, and subtract off the integrated goods.
  VI_U_w <- as.numeric(own_up %in% own_down) - VI_D
  VI_U_r <- as.numeric(own_down %in% own_up) - VI_D
  # matrix form
  own_R_up <- t(sapply(own_down,
                       function(i) as.numeric(own_up %in% i) ))
  own_W_down <- t(sapply(own_up,
                         function(i) as.numeric(own_down %in% i) ))
  # remove integrated goods from partner profits
  own_R_up[,VI_D_idx] <- 0
  own_W_down[,VI_D_idx] <- 0

  # Given the inputs, calculate optimal retail prices
  outtest <- BB::BBoptim(f = bertrand_foc_vert, par = p_R0,
                     own_down = own_down, own_up = own_up,
                     alpha= alpha,
                     delta = delta, cost_r = cost_r,
                     p_W = p_W, cost_w = cost_w, sumFOC = TRUE,
                     control = list(trace=FALSE),
                     quiet = TRUE)

  p_R <- outtest$par

  # calculate shares
  shares <- (exp(delta + alpha*p_R))/(1+sum(exp(delta + alpha*p_R)))

  # counterfactual shares
  shares_tilde <- vector("list",J)

  denom_tilde <- (1-diag(J)) %*% exp(delta + alpha*p_R)
  for (j in (1:J)) {
    shares_tilde[[j]] <- (exp(delta + alpha*p_R))/(1+denom_tilde[j])
    shares_tilde[[j]][j] <- 0
  }
  shares_tilde <- matrix(unlist(shares_tilde), ncol = J, byrow = FALSE)

  # define margins, depends on whether vertically integrated
  margin_up <-   (1-VI_D)*(p_W - cost_w)       + VI_D*(p_R - cost_w - cost_r)
  margin_down <- (1-VI_D)*(p_R - p_W - cost_r) + VI_D*(p_R - cost_w - cost_r)

  # specify payoffs and disagreement payoffs
  pi_w <- own_W %*% (margin_up*shares) + own_W_down %*% (margin_down*shares)

  pi_w_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_W_tilde <- own_W
    own_W_tilde[j,j] <- 0
    temp <- own_W_tilde %*% (margin_up*shares_tilde[,j]) +
      own_R %*% (VI_U_r*margin_down*shares_tilde[,j])
    pi_w_tilde[[j]] <- temp[j]
  }

  pi_r <- own_R %*% (margin_down*shares) + own_R_up %*% (margin_up*shares)

  pi_r_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_R_tilde <- own_R
    own_R_tilde[j,j] <- 0
    temp <- own_R_tilde %*% (margin_down*shares_tilde[,j]) +
      own_W %*% (VI_U_w*margin_up*shares_tilde[,j])
    pi_r_tilde[[j]] <- temp[j]
  }


  NP <- (pi_r - pi_r_tilde)^lambda * (pi_w-pi_w_tilde)^(1-lambda)


  if (showAll == TRUE) {
    out <- (1-sigma) * NP + sigma * (pi_w + pi_r)
    return(out)
  } else {
    out <- -(1-sigma) * NP - sigma * (pi_w + pi_r)  # BBoptim minimizes fn
    out <- out[product_max]
    return(out)
  }

}


#' Nash Product for use in vertical model with sequential timing
#'
#' @noRd

## Vertical function, if p^R cannot adjust in disagreement
## Unused as that is less appealing assumption for sequential timing
## This version allows for overlapping nests.

bargain_NP_vert_seq2_gnl <- function(w_start,product_max,p_W,own_down,own_up,alpha,delta,
                                     cost_w,cost_r,lambda,p_R0,sigma,
                                     a_jk,B,mu,showAll = FALSE){

  p_W[product_max] <- w_start
  J <- length(p_W)

  # construct ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )

  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )

  # indicator of which goods are integrated
  VI_D <- as.numeric(own_down == own_up)
  VI_D_idx <- which(own_down == own_up)

  # VI_U_w = all non-integrated goods sold by upstream divisions of integrated firms. Calculate by finding all upstream goods that are owned by an integrated firm, and subtract off the integrated goods.
  # VI_U_r = all non-integrated goods sold by downstream divisions of integrated firms. Calculate by finding all downstream goods that are sold by an integrated firm, and subtract off the integrated goods.
  VI_U_w <- as.numeric(own_up %in% own_down) - VI_D
  VI_U_r <- as.numeric(own_down %in% own_up) - VI_D
  # matrix form
  own_R_up <- t(sapply(own_down,
                       function(i) as.numeric(own_up %in% i) ))
  own_W_down <- t(sapply(own_up,
                         function(i) as.numeric(own_down %in% i) ))
  # remove integrated goods from partner profits
  own_R_up[,VI_D_idx] <- 0
  own_W_down[,VI_D_idx] <- 0

  # Given the inputs, calculate optimal retail prices
  outtest <- BB::BBoptim(f = bertrand_foc_vert_gnl, par = p_R0,
                     own_down = own_down, own_up = own_up,
                     alpha= alpha,
                     delta = delta, cost_r = cost_r,
                     p_W = p_W, cost_w = cost_w,
                     a_jk=a_jk, B=B, mu=mu, sumFOC = TRUE,
                     control = list(trace=FALSE),
                     quiet = TRUE)

  p_R <- outtest$par

  # calculate shares
  shares <- share_calc(price=p_R,delta,alpha,nest_allocation=a_jk,mu)

  # counterfactual shares
  shares_tilde <- matrix(data = 0, nrow = J, ncol = J)

  for (j in (1:J)) {
    delta_tilde <- delta
    delta_tilde[j] <- -Inf
    shares_tilde[,j] <- share_calc(price=p_R, delta = delta_tilde,
                                   alpha, nest_allocation=a_jk, mu)
  }


  # define margins, depends on whether vertically integrated
  margin_up <-   (1-VI_D)*(p_W - cost_w)       + VI_D*(p_R - cost_w - cost_r)
  margin_down <- (1-VI_D)*(p_R - p_W - cost_r) + VI_D*(p_R - cost_w - cost_r)

  # specify payoffs and disagreement payoffs
  pi_w <- own_W %*% (margin_up*shares) + own_W_down %*% (margin_down*shares)

  pi_w_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_W_tilde <- own_W
    own_W_tilde[j,j] <- 0
    temp <- own_W_tilde %*% (margin_up*shares_tilde[,j]) +
      own_R %*% (VI_U_r*margin_down*shares_tilde[,j])
    pi_w_tilde[[j]] <- temp[j]
  }

  pi_r <- own_R %*% (margin_down*shares) + own_R_up %*% (margin_up*shares)

  pi_r_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_R_tilde <- own_R
    own_R_tilde[j,j] <- 0
    temp <- own_R_tilde %*% (margin_down*shares_tilde[,j]) +
      own_W %*% (VI_U_w*margin_up*shares_tilde[,j])
    pi_r_tilde[[j]] <- temp[j]
  }


  NP <- (pi_r - pi_r_tilde)^lambda * (pi_w-pi_w_tilde)^(1-lambda)


  if (showAll == TRUE) {
    out <- (1-sigma) * NP + sigma * (pi_w + pi_r)
    return(out)
  } else {
    out <- -(1-sigma) * NP - sigma * (pi_w + pi_r)  # BBoptim minimizes fn
    out <- out[product_max]
    return(out)
  }

}
