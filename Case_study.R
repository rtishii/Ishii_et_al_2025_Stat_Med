# -----------------------------------------------------------------------------
# R code for the case study presented by Ishii et al.

# N_g_j: number of patients in group g = 0,1,...,G and stage j = 1,2
# N_g = N_g_1 + N_g_2 for g = 0,1,...,G

# X_g_1: short-term endpoints from patients in group g = 0,1,...,G and stage 1
# X_g = X_g_1 + X_g_2 for g = 0,s

# Y_g_j: long-term endpoints from patients in group g = 0,s and stage j = 1,2
# Y_g = Y_g_1 + Y_g_2 for g = 0,s
# -----------------------------------------------------------------------------

G <- 3
s <- 3

# number of patients in stage 1
N_1_1 <- 28
N_2_1 <- 27
N_3_1 <- 27

# short- and long-term endpoints from patients in stage 1
X_1_1 <- 5
X_2_1 <- 8
X_3_1 <- 9

Y_3_1 <- 14

calc_Pr_Q_x <- function (xi, s) {
  pr <- data.frame(g = c(rep(1, N_1_1 + 1), 
                         rep(2, N_2_1 + 1), 
                         rep(3, N_3_1 + 1)),
                   x = c(0:N_1_1, 0:N_2_1, 0:N_3_1))
  
  
  pr_g1 <- rep(NA, N_1_1 + 1)
  pr_g2 <- rep(NA, N_2_1 + 1)
  pr_g3 <- rep(NA, N_3_1 + 1)
  
  if (s == 1) {
    pr_g1 <- pbinom(floor(0:N_1_1 * N_2_1 / N_1_1), N_2_1, xi[2]) * 
      pbinom(floor(0:N_1_1 * N_3_1 / N_1_1), N_3_1, xi[3])
    for (x in 0:N_2_1) {
      k <- ceiling(x * N_1_1 / N_2_1):N_1_1
      pr_g2[x + 1] <- sum(dbinom(k, N_1_1, xi[1]) * 
                            pbinom(floor(k * N_3_1 / N_1_1), N_3_1, xi[3]))
    }
    for (x in 0:N_3_1) {
      k <- ceiling(x * N_1_1 / N_3_1):N_1_1
      pr_g3[x + 1] <- sum(dbinom(k, N_1_1, xi[1]) * 
                            pbinom(floor(k * N_2_1 / N_1_1), N_2_1, xi[2]))
    }
  }
  
  if (s == 2) {
    pr_g2 <- pbinom(ceiling(0:N_2_1 * N_1_1 / N_2_1) - 1, N_1_1, xi[1]) * 
      pbinom(floor(0:N_2_1 * N_3_1 / N_2_1), N_3_1, xi[3])
    for (x in 0:N_1_1) {
      k <- (floor(x * N_2_1 / N_1_1) + 1):N_2_1
      pr_g1[x + 1] <- sum(dbinom(k, N_2_1, xi[2]) * 
                            pbinom(floor(k * N_3_1 / N_2_1), N_3_1, xi[3]))
    }
    for (x in 0:N_3_1) {
      k <- ceiling(x * N_2_1 / N_3_1):N_2_1
      pr_g3[x + 1] <- sum(dbinom(k, N_2_1, xi[2]) * 
                            pbinom(ceiling(k * N_1_1 / N_2_1) - 1, 
                                   N_1_1, xi[1]))
    }
  }
  
  if (s == 3) {
    pr_g3 <- pbinom(ceiling(0:N_3_1 * N_1_1 / N_3_1) - 1, N_1_1, xi[1]) * 
      pbinom(ceiling(0:N_3_1 * N_2_1 / N_3_1) - 1, N_2_1, xi[2])
    for (x in 0:N_1_1) {
      k <- (floor(x * N_3_1 / N_1_1) + 1):N_3_1
      pr_g1[x + 1] <- sum(dbinom(k, N_3_1, xi[3]) * 
                            pbinom(ceiling(k * N_2_1 / N_3_1) - 1, 
                                   N_2_1, xi[2]))
    }
    for (x in 0:N_2_1) {
      k <- floor(ceiling(x * N_3_1 / N_2_1) + 1):N_3_1
      pr_g2[x + 1] <- sum(dbinom(k, N_3_1, xi[3]) * 
                            pbinom(ceiling(k * N_1_1 / N_3_1) - 1,
                                   N_1_1, xi[1]))
    }
  }
  
  pr$prob <- c(pr_g1, pr_g2, pr_g3)
  
  return(pr)
}

alp_val <- function (xi, pi){
  LL <- max(-(1 - xi) / (1 + pi - xi), - xi / (1 + xi - pi))
  if (pi >  xi) UL <- xi / (pi - xi)
  if (pi <  xi) UL <- (1 - xi) / (xi - pi)
  if (pi == xi) UL <- 1e10
  
  c(LL, UL)
}

lik_alp <- function (xi, pi, alp) {
  beta <- xi + alp * (xi - pi)
  
  L <- pi ^ (p01 + p11) * (1 - pi) ^ (p00 + p10) *
    (1 + alp) ^ (-1) * (1 - beta) ^ p01 * (1 - beta + alp) ^ p00 *
    (beta + alp) ^ p11 * beta ^ p10
  
  if (!is.nan(L) & L > 0) loglik <- log(L)
  else loglik <- -10 ^ 10
  
  -loglik # negative log-likelihood
}

calc_pval <- function (x, delta) {
  pis <- x[G + 1]
  alp <- x[G + 2]
  y0 <- x[G + 3]
  ys <- x[G + 4]
  xis <- x[s]
  
  if (s == 1) {
    n1 <- N_1_1
    n <- N_1
  }
  if (s == 2) {
    n1 <- N_2_1
    n <- N_2
  }
  if (s == 3) {
    n1 <- N_3_1
    n <- N_3
  }
  
  pr_Q_y1 <- numeric(n1 + 1)
  beta <- xis + alp * (xis - pis)
  q1 <- (beta + alp) / (1 + alp)
  q2 <- beta / (1 + alp)
  
  c0 <- calc_Pr_Q_x(x[1:G], s)
  c <- c0$prob[c0$g == s]
  
  for (r in 0:n1) {
    for (x in 0:n1) {
      jdat <- max(0, x + r - n1):min(x, r)
      
      pr_Q_y1[r + 1] <- pr_Q_y1[r + 1] + 
        sum(dbinom(jdat, r, q1) * dbinom(x - jdat, n1 - r, q2)) * c[x + 1]
    }
  }
  
  t <- y0 + ys
  kr <- max(0, t - N_0):min(t, n)
  h <- numeric(length(kr))
  
  for(kn in 1:length(kr)) {
    rr <- max(0, kr[kn] - (n - n1)):min(kr[kn], n1)
    h[kn] <- sum(choose(n1, rr) * choose(n - n1, kr[kn] - rr) * 
                   pr_Q_y1[rr + 1])
  }
  
  h <- choose(N_0, t - kr) * h * exp(delta * kr)
  
  p0 <- h / sum(h)  
  pval_exact_lower <- sum(p0[kr >= ys])
  pval_exact_upper <- sum(p0[kr <= ys])
  pvals <- c(pval_exact_lower, pval_exact_lower - p0[kr == ys] * 0.5,
             pval_exact_upper, pval_exact_upper - p0[kr == ys] * 0.5)
  names(pvals) <- c("exact_lower", "midp_lower", "exact_upper", "midp-upper")
  
  return(pvals)
}

main <- function() {
  
  xi_MLE <- c(X_1_1 / N_1_1, X_2_1 / N_2_1, X_3 / N_3)
  pi_0_MLE <- Y_0 / N_0
  pi_s_MLE <- Y_3 / N_3
  alpha_MLE <- (p11 - xi_MLE[s] * pi_s_MLE) / 
    (pi_s_MLE * (1 - pi_s_MLE + xi_MLE[s]) - p11)
  
  cprob <- calc_Pr_Q_x(xi_MLE, s)
  Pr_Q <- sum(subset(cprob, g == s)$prob * dbinom(0:N_3_1, N_3_1, xi_MLE[3]))
  
  E <- rep(NA, G)
  E[1] <- sum(0:N_1_1 * dbinom(0:N_1_1, N_1_1, xi_MLE[1]) * 
                subset(cprob, g == 1)$prob) / (N_1_1 * Pr_Q)
  E[2] <- sum(0:N_2_1 * dbinom(0:N_2_1, N_2_1, xi_MLE[2]) * 
                subset(cprob, g == 2)$prob) / (N_2_1 * Pr_Q)
  E[3] <- (sum(0:N_3_1 * dbinom(0:N_3_1, N_3_1, xi_MLE[3]) * 
                 subset(cprob, g == 3)$prob) / Pr_Q + N_3_2 * xi_MLE[3]) / N_3
  
  xi_CMAE <- xi_MLE - (E - xi_MLE)
  
  
  beta <- xi_MLE[s] + alpha_MLE * (xi_MLE[s] - pi_s_MLE)
  q1 <- (beta + alpha_MLE) / (1 + alpha_MLE)
  q2 <- beta / (1 + alpha_MLE)
  
  pr_Q_y1 <- numeric(N_3_1 + 1)
  c <- subset(cprob, g == s)$prob
  for (r in 0:N_3_1) {
    for (x in 0:N_3_1) {
      jdat <- max(0, x + r - N_3_1):min(x, r)
      
      pr_Q_y1[r + 1] <- pr_Q_y1[r + 1] + 
        sum(dbinom(jdat, r, q1) * dbinom(x - jdat, N_1_1 - r, q2)) * c[x + 1]
    }
  }
  
  E <- (sum((0:N_3_1) * dbinom(0:N_3_1, N_3_1, pi_s_MLE) * pr_Q_y1) / Pr_Q + 
          N_3_2 * pi_s_MLE) / N_3
  pi_s_CMAE <- pi_s_MLE - (E - pi_s_MLE)
  
  xs <- X_3
  ys <- Y_3
  xr <- (0:N_3_1)[0:N_3_1 / N_3_1 > max(X_1_1 / N_1_1, X_2_1 / N_2_1)]
  
  k <- max(0, xs - N_3_2):min(N_3_1, xs)
  q <- choose(N_3_1, k) * choose(N_3_2, xs - k) * 
    (k / N_3_1 > max(X_1_1 / N_1_1, X_2_1 / N_2_1))
  xi_UMVCUE <- xi_MLE
  xi_UMVCUE[s] <- (xs - sum(k * q) / sum(q)) / N_3_2
  
  c1 <- yy <- c()
  for (x1 in xr) {
    x2 <- xs - x1
    yr <- max(0, ys - N_3_2, x2 + ys - zs - N_3_2):min(N_3_1, ys, N_3_1 + zs - x1)
    yr <- yr[pmin(yr, x1) + pmin(ys - yr, x2) >= zs]
    yy <- c(yy, yr)
    for (y1 in yr) {
      y2 <- ys - y1
      z1 <- (zs - min(x2, y2)):min(x1, y1, zs)
      c0 <- sum(choose(y1, z1) * 
                  choose(N_3_1 - y1, x1 - z1) * 
                  choose(y2, zs - z1) * 
                  choose(N_3_2 - y2, x2 - zs + z1))
      c1 <- c(c1, c0)
    }
  }
  
  c2 <- choose(N_3_2, ys - yy) * choose(N_3_1, yy)
  
  if (sum(is.infinite(c1 * c2)) > 0) {
    dmean <- mean(log10(c1) + log10(c2))
    c1 <- c1 * 10 ^ (-dmean)
    c2 <- c2 * 10 ^ (-dmean)
  }
  
  pi_s_UMVCUE <- (ys - sum(yy * c1 * c2) / sum(c1 * c2)) / N_3_2
  
  alpha_CMAE <- optimize(
    function (x) lik_alp(xi_CMAE[s], pi_s_CMAE, x),
    interval = alp_val(xi_CMAE[s], pi_s_CMAE))$minimum
  
  alpha_UMVCUE <- optimize(
    function (x) lik_alp(xi_UMVCUE[s], pi_s_UMVCUE, x),
    interval = alp_val(xi_UMVCUE[s], pi_s_UMVCUE))$minimum
  
  val <- c(xi_MLE, pi_s_MLE, alpha_MLE, Y_0, Y_3)
  
  pval_MLE <- calc_pval(val, delta = 0)
  
  LL_exact_MLE <- uniroot(function (d) calc_pval(val, d)[1] - 0.025, 
                          d_int)$root
  LL_midp_MLE <- uniroot(function (d) calc_pval(val, d)[2] - 0.025, 
                         d_int)$root
  UL_exact_MLE <- uniroot(function (d) calc_pval(val, d)[3] - 0.025, 
                          d_int)$root
  UL_midp_MLE <- uniroot(function (d) calc_pval(val, d)[4] - 0.025, 
                         d_int)$root
  
  val <- c(xi_CMAE, pi_s_CMAE, alpha_CMAE, Y_0, Y_3)
  
  pval_CMAE <- calc_pval(val, delta = 0)
  
  LL_exact_CMAE <- uniroot(function (d) calc_pval(val, d)[1] - 0.025, 
                           d_int)$root
  LL_midp_CMAE <- uniroot(function (d) calc_pval(val, d)[2] - 0.025, 
                          d_int)$root
  UL_exact_CMAE <- uniroot(function (d) calc_pval(val, d)[3] - 0.025, 
                           d_int)$root
  UL_midp_CMAE <- uniroot(function (d) calc_pval(val, d)[4] - 0.025, 
                          d_int)$root
  
  val <- c(xi_UMVCUE, pi_s_UMVCUE, alpha_UMVCUE, Y_0, Y_3)
  
  pval_UMVCUE <- calc_pval(val, delta = 0)
  
  LL_exact_UMVCUE <- uniroot(function (d) calc_pval(val, d)[1] - 0.025, 
                             d_int)$root
  LL_midp_UMVCUE <- uniroot(function (d) calc_pval(val, d)[2] - 0.025, 
                            d_int)$root
  UL_exact_UMVCUE <- uniroot(function (d) calc_pval(val, d)[3] - 0.025, 
                             d_int)$root
  UL_midp_UMVCUE <- uniroot(function (d) calc_pval(val, d)[4] - 0.025, 
                            d_int)$root
  
  
  res0 <-
    data.frame(t(xi_MLE), xi_MLE[s], 
               t(xi_CMAE), xi_CMAE[s], t(xi_UMVCUE), xi_UMVCUE[s],
               pi_0_MLE, pi_s_MLE, pi_s_CMAE, pi_s_UMVCUE,
               alpha_MLE, alpha_CMAE, alpha_UMVCUE,
               t(pval_MLE[1:2]), t(pval_CMAE[1:2]), t(pval_UMVCUE[1:2]), 
               LL_exact_MLE, LL_midp_MLE, UL_exact_MLE, UL_midp_MLE,
               LL_exact_CMAE, LL_midp_CMAE, UL_exact_CMAE, UL_midp_CMAE,
               LL_exact_UMVCUE, LL_midp_UMVCUE,
               UL_exact_UMVCUE, UL_midp_UMVCUE)
  
  
  fr <- function (x, d) format(round(x, d), nsmall = d)
  
  col1 <- c("MLE", "", "", "CMAE", "", "", "UMVCUE", "", "")
  col2 <- rep(c("Estimate", "CI (exact)", "CI (mid-p)"), 3)
  
  col3 <- rep(c(fr(pi_0_MLE * 100, 1), "", ""), 3)
  
  col4 <- c(fr(pi_s_MLE * 100, 1),
            paste0("[", fr(LL_exact_MLE, 2), ", ", fr(UL_exact_MLE, 2), "]"),
            paste0("[", fr(LL_midp_MLE, 2), ", ", fr(UL_midp_MLE, 2), "]"),
            fr(pi_s_CMAE * 100, 1),
            paste0("[", fr(LL_exact_CMAE, 2), ", ", fr(UL_exact_CMAE, 2), "]"),
            paste0("[", fr(LL_midp_CMAE, 2), ", ", fr(UL_midp_CMAE, 2), "]"),
            fr(pi_s_UMVCUE * 100, 1),
            paste0("[", fr(LL_exact_UMVCUE, 2), ", ",
                   fr(UL_exact_UMVCUE, 2), "]"),
            paste0("[", fr(LL_midp_UMVCUE, 2), ", ",
                   fr(UL_midp_UMVCUE, 2), "]")
  )
  
  pval <- c(NA, pval_MLE[1:2], NA, pval_CMAE[1:2], NA, pval_UMVCUE[1:2])
  
  col5 <- fr(pval, 2)
  col5[pval < 0.001] <- "<0.001"
  col5[is.na(pval)] <- ""
  
  out <- data.frame(Method = col1,
                    Parameter = col2,
                    Placebo = col3,
                    BARI = col4,
                    Pvalue = col5)

  out
}

calc_z <- function (N, X, Y, rho) {
  z <- max(0, X + Y - N):min(X, Y)
  rhos <- (N * z - X * Y) / sqrt(X * Y * (N - X) * (N - Y))
  z[abs(rho - rhos) == min(abs(rho - rhos))]
}

# -------------------------
# Original data
# -------------------------
N_0 <- 189
Y_0 <- 10

N_3_2 <- 254
X_3_2 <- 85
Y_3_2 <- 85

N_3 <- N_3_1 + N_3_2
X_3 <- X_3_1 + X_3_2
Y_3 <- Y_3_1 + Y_3_2

zs <- calc_z(N_3, X_3, Y_3, 0.6)

p00 <- (N_3 - X_3 - Y_3 + zs) / N_3
p01 <- (Y_3 - zs) / N_3
p10 <- (X_3 - zs) / N_3
p11 <- zs / N_3

d_int <- c(-4, 4)
main()


# -------------------------
# Small data (tau = 0.5)
# -------------------------
N_0 <- 56
Y_0 <- 3

N_3_2 <- 27
X_3_2 <- 9
Y_3_2 <- 9

N_3 <- N_3_1 + N_3_2
X_3 <- X_3_1 + X_3_2
Y_3 <- Y_3_1 + Y_3_2

zs <- calc_z(N_3, X_3, Y_3, 0.6)

p00 <- (N_3 - X_3 - Y_3 + zs) / N_3
p01 <- (Y_3 - zs) / N_3
p10 <- (X_3 - zs) / N_3
p11 <- zs / N_3

d_int <- c(-5, 5)
main()


# -------------------------
# Small data (tau = 0.75)
# -------------------------
N_0 <- 37
Y_0 <- 2

N_3_2 <- 9
X_3_2 <- 3
Y_3_2 <- 3

N_3 <- N_3_1 + N_3_2
X_3 <- X_3_1 + X_3_2
Y_3 <- Y_3_1 + Y_3_2

zs <- calc_z(N_3, X_3, Y_3, 0.6)

p00 <- (N_3 - X_3 - Y_3 + zs) / N_3
p01 <- (Y_3 - zs) / N_3
p10 <- (X_3 - zs) / N_3
p11 <- zs / N_3

d_int <- c(-5, 5)
main()
