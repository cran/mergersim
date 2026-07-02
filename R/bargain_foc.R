#' Nash bargaining first-order conditions
#'
#' @param price Price
#' @param own Ownership matrix
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost Marginal costs for each product
#' @param lambda Bargaining power of the buyer
#' @param includeMUI logical; whether to include marginal utility of income
#' in buyer's payoff, thereby translating dollars to utility. Default is True,
#' interpreted as buyer maximizing utility. Setting equal to False would have
#' interpretation that buyer maximizes profts.
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
#' own_pre = diag(3)
#' p0 <- c_j*1.1
#'
#' bertrand_foc(price = p0,
#' own = own_pre, alpha= alpha,
#' delta = delta, cost = c_j)
#'
#' @export



##################################################################
# Nash bargaining first-order conditions
##################################################################

bargain_foc <- function(price,own,alpha,delta,cost,
                     lambda,includeMUI=TRUE){

  wshares <- (exp(delta + alpha*price))/(1+sum(exp(delta + alpha*price)))
  m <- price - cost
  ownd <- alpha*wshares*(1-wshares)
  crossd <- -alpha*wshares%*%t(wshares)
  dd <- crossd
  diag(dd) <- ownd

  bert_foc <- (own * t(dd)) %*% m + wshares

  own_excl <- own
  diag(own_excl) <- 0

  J <- length(price)
  wshares_tilde <- vector("list",J) # counterfactual shares. Assume firm still offers other goods
  denom_tilde <- (1-diag(J)) %*% exp(delta + alpha*price)
  for (j in (1:J)) {
    wshares_tilde[[j]] <- (exp(delta + alpha*price))/(1+denom_tilde[j])
    wshares_tilde[[j]][j] <- 0
  }
  wshares_tilde <- matrix(unlist(wshares_tilde), ncol = J, byrow = FALSE)

  pi_w <- m*wshares

  delta_share <- wshares_tilde - wshares
  diag(delta_share) <- 0
  pi_w_tilde <- diag(own_excl %*% (matrix(m, ncol = J, nrow = J, byrow = FALSE) * delta_share))

  if (includeMUI == TRUE) {
    MUI <- (-1/alpha)
    dpi_r_dp <- -wshares
  } else {
    MUI <- 1
    dpi_r_dp <- alpha*wshares
  }

  pi_r <- MUI*log(1+sum(exp(delta + alpha*price)))
  pi_r_tilde <- MUI*log(1+sum(exp(delta + alpha*price))-exp(delta + alpha*price))

  foc <- lambda*dpi_r_dp*(pi_w-pi_w_tilde) +
    (1-lambda)*(pi_r - pi_r_tilde)*bert_foc

  return(foc)
}

