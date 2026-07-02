#' Calibrate costs for Nash bargain in vertical model with sequential timing
#'
#' @param c_w_val Wholesale/upstream costs
#' @param price_w Upstream or wholesale prices
#' @param own_down Ownership matrix for downstream firms
#' @param own_up Ownership matrix for upstream firms
#' @param alpha Price coefficient
#' @param delta Mean values
#' @param cost_r Marginal costs for downstream firm for each product
#' @param lambda Bargaining power of the buyer/retailer
#' @param price_r Retail prices
#' @param sigma Contract type; value between 0 and 1 where 0 is linear price,
#'  1 is two-part tariff
#' @param setTol tolerance for convergence
#' @param setMaxIter Maximum iterations for convergence loop
#' @param showAll logical; if true, returns gains for trade for every product
#' @param symmetricCosts Constrain costs to be symmetric
#'
#' @returns The value of objective function
#'
#' @details This function can be used to calibrate the costs in a Nash bargain
#' which is the upstream market of a vertical supply chain. Assumes logit demand.
#'
#' @examples
#' bargain_vert_seq_calibrate(c_w_val = rep(.2, times = 6),
#' price_w = c(1.6, 1.6, 1.6, 1.6, 1.6, 1.6),
#' own_down = paste0("R",rep(c(1,2,3),each=2)),
#' own_up = paste0("W",rep(c(1,2),3)),
#' alpha = -0.9, delta = c(0.2, 0.3, 0.9, 1.0, 0.8, 0.9),
#' lambda = 0.5, cost_r = rep(.1, times = 6),
#' price_r = c(2.9, 2.9, 3.0, 3.0, 3.0, 3.0),
#' sigma = 0, setTol = 1.5)
#'
#' @export



##################################################################
# Nash Product
##################################################################

#### Calibration functions based on seq, which allows p^R to adjust in disagreement
## Given demand parameters, what are wholesale costs that allow the model wholesale
## prices to closely match the observed wholesale prices.

bargain_vert_seq_calibrate <- function(c_w_val,price_w,own_down,
                                       own_up, alpha,delta,
                                       lambda,cost_r,price_r,sigma,
                                       setTol = 0.01,
                                       setMaxIter = 500,
                                       showAll = FALSE,
                                       symmetricCosts = FALSE){

  J <- length(price_w)

  if (symmetricCosts == FALSE) {

    error <- rep(1,J)
    tol <- setTol

    p_W0 <- price_w + .1
    p_R0 <- price_r
    iter <- 1

    while (max(error) > tol & iter < setMaxIter) {
      for (x in 1:J) {

        #w_start <- p_W0[x] - .001   # add eps to get cleaner convergence flags
        w_start <- p_W0[x]
        # for calibration function only, hard code p_R0 = price_r

        lbc <- c_w_val[x]
        ubc <- price_r[x]
        checktest <- stats::optimize(f = bargain_NP_vert_seq,
                              product_max = x, price_w = p_W0,
                              own_down = own_down, own_up = own_up,
                              alpha= alpha, delta = delta,
                              cost_w = c_w_val, cost_r = cost_r, lambda = lambda,
                              p_R0 = price_r, sigma = sigma, showAll = FALSE,
                              lower = lbc, upper = ubc)

        pw_test <- checktest$minimum

        error[x] <- abs(pw_test - p_W0[x])

        p_W0[x] <- pw_test
        #print(error)
        #print(p_W0)

        # recover price_r at these price_w and update r_R0
        outtest_r <- BB::BBoptim(f = bertrand_foc_vert, par = p_R0,
                             own_down = own_down, own_up = own_up,
                             alpha= alpha,
                             delta = delta, cost_r = cost_r,
                             price_w = p_W0, cost_w = c_w_val, sumFOC = TRUE,
                             control = list(trace=FALSE),
                             quiet = TRUE)

        p_R0 <- outtest_r$par

        iter <- iter + 1
        if (iter > setMaxIter) {warning("Max iterations reached")}
      }
    }

  }

  if (symmetricCosts == TRUE) {
    c_w_val2 <- rep(c_w_val, length(delta))
    error <- rep(1,J)
    tol <- setTol

    p_W0 <- price_w
    p_R0 <- price_r
    iter <- 1

    while (max(error) > tol & iter < setMaxIter) {
      for (x in 1:J) {

        #w_start <- p_W0[x] - .001   # add eps to get cleaner convergence flags
        w_start <- p_W0[x]
        # for calibration function only, hard code p_R0 = price_r
        checktest <- stats::optimize(f = bargain_NP_vert_seq,
                              product_max = x, price_w = p_W0,
                              own_down = own_down, own_up = own_up,
                              alpha= alpha, delta = delta,
                              cost_w = c_w_val2, cost_r = cost_r, lambda = lambda,
                              p_R0 = price_r, sigma = sigma, showAll = FALSE,
                              lower = 0, upper = 5)


        pw_test <- checktest$minimum

        error[x] <- abs(pw_test - p_W0[x])
        #print(error)

        p_W0[x] <- pw_test


        # recover price_r at these price_w and update r_R0
        outtest_r <- BB::BBoptim(f = bertrand_foc_vert, par = p_R0,
                             own_down = own_down, own_up = own_up,
                             alpha= alpha,
                             delta = delta, cost_r = cost_r,
                             price_w = p_W0, cost_w = c_w_val2, sumFOC = TRUE,
                             control = list(trace=FALSE),
                             quiet = TRUE)

        p_R0 <- outtest_r$par

        iter <- iter + 1
        if (iter > setMaxIter) {warning("Max iterations reached")}
      }
    }
  }

  # uncomment next two lines for debugging purposes.
  #print(iter)
  #print(error)
  #print(p_W0)
  #print(p_R0)

  price_w2 <- p_W0
  price_r2 <- p_R0

  shares2 <- (exp(delta + alpha*price_r2))/(1+sum(exp(delta + alpha*price_r2)))

  #print(sum((price_w2 - price_w)^2))

  if (showAll == TRUE) {
    return(list("price_w" = price_w2,"price_r" = price_r2,"shares" = as.numeric(shares2)) )
  } else {
    out <- sum((price_w2 - price_w)^2)
    return(out)
  }

}
