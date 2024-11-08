# -----------------------------------------------------------------------------
# R code for the simulation study presented by Ishii et al.
# -----------------------------------------------------------------------------

library(doParallel)

# ----- scenario -----
n <- data.frame(n = c(50, 200))
rho1 <- data.frame(rho = c(0.2, 0.4, 0.6))
rho2 <- data.frame(rho = c(0.2, 0.4, 0.5))
tau <- data.frame(tau = c(0.25, 0.5, 0.75))
G <- 4:2

xi_H0_1 <- pi_H0_1 <- c(0.1, 0.5, 0.7)
xi_H0_2 <- c(0.3, 0.7)
pi_H0_2 <- c(0.5, 0.5)
xi_0_H1 <- pi_0_H1 <- 0.5
xi_G_H1 <- pi_G_H1 <- c(0.78, 0.645)

scenario <- c()
for (i in 1:length(G)) {
  dat1 <- matrix(NA, length(xi_H0_1), max(G) + 1)
  colnames(dat1) <- paste0("xi_", 0:max(G))
  dat1[1:(length(xi_H0_1) * (G[i] + 1))] <- xi_H0_1
  dat2 <- matrix(NA, length(pi_H0_1), max(G) + 1)
  colnames(dat2) <- paste0("pi_", 0:max(G))
  dat2[1:(length(pi_H0_1) * (G[i] + 1))] <- pi_H0_1
  
  dat3 <- merge(rho1, cbind(dat1, dat2))
  
  dat3 <- merge(dat3, tau)
  dat3 <- merge(dat3, n)
  dat3$G <- G[i]
  dat3$simn <- 100000
  dat3$H <- "H_0_1"
  
  scenario <- rbind(scenario, dat3)
}

for (i in 1:length(G)) {
  dat1 <- matrix(NA, length(xi_H0_2), max(G) + 1)
  colnames(dat1) <- paste0("xi_", 0:max(G))
  dat1[1:(length(xi_H0_2) * (G[i] + 1))] <- xi_H0_2
  dat2 <- matrix(NA, length(pi_H0_2), max(G) + 1)
  colnames(dat2) <- paste0("pi_", 0:max(G))
  dat2[1:(length(pi_H0_2) * (G[i] + 1))] <- pi_H0_2
  
  dat3 <- merge(rho2, cbind(dat1, dat2))
  
  dat3 <- merge(dat3, tau)
  dat3 <- merge(dat3, n)
  dat3$G <- G[i]
  dat3$simn <- 100000
  dat3$H <- "H_0_2"
  
  scenario <- rbind(scenario, dat3)
}

for (i in 1:length(G)) {
  for (j in 1:nrow(n)) {
    dat1 <- matrix(NA, 1, max(G) + 1)
    colnames(dat1) <- paste0("xi_", 0:max(G))
    dat1[1:G[i]] <- xi_0_H1
    dat1[G[i] + 1] <- xi_G_H1[j]
    dat2 <- matrix(NA, 1, max(G) + 1)
    colnames(dat2) <- paste0("pi_", 0:max(G))
    dat2[1:G[i]] <- pi_0_H1
    dat2[G[i] + 1] <- pi_G_H1[j]
    
    dat3 <- merge(rho1, cbind(dat1, dat2))
    
    dat3 <- merge(dat3, tau)
    dat3 <- merge(dat3, data.frame(n = n[j, 1]))
    dat3$G <- G[i]
    dat3$simn <- 10000
    dat3$H <- "H_1_1"
    
    scenario <- rbind(scenario, dat3)
  }
}

for (i in 1:length(G)) {
  for (j in 1:nrow(n)) {
    dat1 <- matrix(NA, 1, max(G) + 1)
    colnames(dat1) <- paste0("xi_", 0:max(G))
    dat1[1:(G[i] + 1)] <-  xi_0_H1 + (xi_G_H1[j] - xi_0_H1) * (0:G[i]) / G[i]
    dat2 <- matrix(NA, 1, max(G) + 1)
    colnames(dat2) <- paste0("pi_", 0:max(G))
    dat2[1:(G[i] + 1)] <-   pi_0_H1 + (pi_G_H1[j] - pi_0_H1) * (0:G[i]) / G[i]
    
    dat3 <- merge(rho1, cbind(dat1, dat2))
    
    dat3 <- merge(dat3, tau)
    dat3 <- merge(dat3, data.frame(n = n[j, 1]))
    dat3$G <- G[i]
    dat3$simn <- 10000
    dat3$H <- "H_1_2"
    
    scenario <- rbind(scenario, dat3)
  }
}

for (i in 1:length(G)) {
  for (j in 1:nrow(n)) {
    dat1 <- matrix(NA, 1, max(G) + 1)
    colnames(dat1) <- paste0("xi_", 0:max(G))
    dat1[1] <- xi_0_H1
    dat1[2:(G[i] + 1)] <- xi_G_H1[j]
    dat2 <- matrix(NA, 1, max(G) + 1)
    colnames(dat2) <- paste0("pi_", 0:max(G))
    dat2[1] <- pi_0_H1
    dat2[2:(G[i] + 1)] <- pi_G_H1[j]
    
    dat3 <- merge(rho1, cbind(dat1, dat2))
    
    dat3 <- merge(dat3, tau)
    dat3 <- merge(dat3, data.frame(n = n[j, 1]))
    dat3$G <- G[i]
    dat3$simn <- 10000
    dat3$H <- "H_1_3"
    
    scenario <- rbind(scenario, dat3)
  }
}
scenario$scenario <- 1:nrow(scenario)

scenario <- scenario[c("scenario", "simn", "H", "G", "n", "tau", "rho", 
                       paste0("xi_", 0:max(G)), paste0("pi_", 0:max(G)))]

rm(list = setdiff(ls(), "scenario"))


# ----- simulation -----

# cores: number of cores
cores <- detectCores()

# return Pr(Q = s | X_g^(1) = x) for x = 0,1,...,n1 and g = 1,...,G
calc_Pr_Q_x <- function (xi, G, s) {
  k <- 0:n1
  pr <- data.frame(g = rep(1:G, each = n1 + 1),
                   x = rep(k, G))
  
  if (G == 2) {
    if (s == 1) {
      pr_g1 <- pbinom(k, n1, xi[2])
      pr_g2 <- 1 - pbinom(k - 1, n1, xi[1])
    }
    
    if (s == 2) {
      pr_g1 <- 1 - pbinom(k, n1, xi[2])
      pr_g2 <- pbinom(k - 1, n1, xi[1])
    }
    
    pr$prob <- c(pr_g1, pr_g2)
  }
  
  if (G == 3) {
    pr_g1 <- pr_g2 <- pr_g3 <- rep(NA, n1 + 1)
    
    if (s == 1) {
      pr_g1 <- pbinom(k, n1, xi[2]) * pbinom(k, n1, xi[3])
      for (x in 0:n1) {
        pr_g2[x + 1] <- sum(dbinom(x:n1, n1, xi[1]) * pbinom(x:n1, n1, xi[3]))
        pr_g3[x + 1] <- sum(dbinom(x:n1, n1, xi[1]) * pbinom(x:n1, n1, xi[2]))
      }
    }
    
    if (s == 2) {
      pr_g2 <- pbinom(k - 1, n1, xi[1]) * pbinom(k, n1, xi[3])
      for (x in 0:n1) {
        pr_g1[x + 1] <- sum(dbinom((x + 1):n1, n1, xi[2]) * 
                              pbinom((x + 1):n1, n1, xi[3]))
        pr_g3[x + 1] <- sum(dbinom(x:n1, n1, xi[2]) * 
                              pbinom(x:n1 - 1, n1, xi[1]))
      }
    }
    
    if (s == 3) {
      pr_g3 <- pbinom(k - 1, n1, xi[1]) * pbinom(k - 1, n1, xi[2])
      for (x in 0:n1) {
        pr_g1[x + 1] <- sum(dbinom((x + 1):n1, n1, xi[3]) * 
                              pbinom((x + 1):n1 - 1, n1, xi[2]))
        pr_g2[x + 1] <- sum(dbinom((x + 1):n1, n1, xi[3]) * 
                              pbinom((x + 1):n1 - 1, n1, xi[1]))
      }
    }
    
    pr$prob <- c(pr_g1, pr_g2, pr_g3)
  }
  
  if (G == 4) {
    pr_g1 <- pr_g2 <- pr_g3 <- pr_g4 <- rep(NA, n1 + 1)
    
    if (s == 1) {
      pr_g1 <- pbinom(k, n1, xi[2]) * 
        pbinom(k, n1, xi[3]) * pbinom(k, n1, xi[4])
      for (x in 0:n1) {
        pr_g2[x + 1] <- 
          sum(dbinom(x:n1, n1, xi[1]) * 
                pbinom(x:n1, n1, xi[3]) * pbinom(x:n1, n1, xi[4]))
        pr_g3[x + 1] <- 
          sum(dbinom(x:n1, n1, xi[1]) * 
                pbinom(x:n1, n1, xi[2]) * pbinom(x:n1, n1, xi[4]))
        pr_g4[x + 1] <- 
          sum(dbinom(x:n1, n1, xi[1]) * 
                pbinom(x:n1, n1, xi[2]) * pbinom(x:n1, n1, xi[3]))
      }
    }
    
    if (s == 2) {
      pr_g2 <- pbinom(k - 1, n1, xi[1]) *
        pbinom(k, n1, xi[3]) * pbinom(k, n1, xi[4])
      for (x in 0:n1) {
        pr_g1[x + 1] <- 
          sum(dbinom((x + 1):n1, n1, xi[2]) * 
                pbinom((x + 1):n1, n1, xi[3]) * pbinom((x + 1):n1, n1, xi[4]))
        pr_g3[x + 1] <- 
          sum(dbinom(x:n1, n1, xi[2]) * 
                pbinom(x:n1 - 1, n1, xi[1]) * pbinom(x:n1, n1, xi[4]))
        pr_g4[x + 1] <- 
          sum(dbinom(x:n1, n1, xi[2]) * 
                pbinom(x:n1 - 1, n1, xi[1]) * pbinom(x:n1, n1, xi[3]))
      }
    }
    
    if (s == 3) {
      pr_g3 <- pbinom(k - 1, n1, xi[1]) * 
        pbinom(k - 1, n1, xi[2]) * pbinom(k, n1, xi[4])
      for (x in 0:n1) {
        pr_g1[x + 1] <- 
          sum(dbinom((x + 1):n1, n1, xi[3]) * 
                pbinom((x + 1):n1 - 1, n1, xi[2]) * 
                pbinom((x + 1):n1, n1, xi[4]))
        pr_g2[x + 1] <- 
          sum(dbinom((x + 1):n1, n1, xi[3]) * 
                pbinom((x + 1):n1 - 1, n1, xi[1]) * 
                pbinom((x + 1):n1, n1, xi[4]))
        pr_g4[x + 1] <- 
          sum(dbinom(x:n1, n1, xi[3]) * 
                pbinom(x:n1 - 1, n1, xi[1]) * pbinom(x:n1 - 1, n1, xi[2]))
      }
    }
    
    if (s == 4) {
      pr_g4 <- pbinom(k - 1, n1, xi[1]) * 
        pbinom(k - 1, n1, xi[2]) * pbinom(k - 1, n1, xi[3])
      for (x in 0:n1) {
        pr_g1[x + 1] <- sum(dbinom((x + 1):n1, n1, xi[4]) * 
                              pbinom((x + 1):n1 - 1, n1, xi[2]) * 
                              pbinom((x + 1):n1 - 1, n1, xi[3]))
        pr_g2[x + 1] <- sum(dbinom((x + 1):n1, n1, xi[4]) * 
                              pbinom((x + 1):n1 - 1, n1, xi[1]) * 
                              pbinom((x + 1):n1 - 1, n1, xi[3]))
        pr_g3[x + 1] <- sum(dbinom((x + 1):n1, n1, xi[4]) * 
                              pbinom((x + 1):n1 - 1, n1, xi[1]) * 
                              pbinom((x + 1):n1 - 1, n1, xi[2]))
      }
    }
    
    pr$prob <- c(pr_g1, pr_g2, pr_g3, pr_g4)
  }
  
  return(pr)
}

# return Pr(Q = s | X_s^(1) = x) for x = 0,1,...,n1
calc_Pr_Q_x_s <- function (xi, s) {
  if (G == 2){
    if (s == 1) return(pbinom(0:n1, n1, xi[2]))
    if (s == 2) return(pbinom(0:n1 - 1, n1, xi[1]))
  }
  if (G == 3){
    if (s == 1) 
      return(pbinom(0:n1, n1, xi[2]) * pbinom(0:n1, n1, xi[3]))
    if (s == 2) 
      return(pbinom(0:n1 - 1, n1, xi[1]) * pbinom(0:n1, n1, xi[3]))
    if (s == 3) 
      return(pbinom(0:n1 - 1, n1, xi[1]) * pbinom(0:n1 - 1, n1, xi[2]))
  }
  if (G == 4){
    if (s == 1) 
      return(pbinom(0:n1, n1, xi[2]) * pbinom(0:n1, n1, xi[3]) * 
               pbinom(0:n1, n1, xi[4]))
    if (s == 2) 
      return(pbinom(0:n1 - 1, n1, xi[1]) * pbinom(0:n1, n1, xi[3]) * 
               pbinom(0:n1, n1, xi[4]))
    if (s == 3) 
      return(pbinom(0:n1 - 1, n1, xi[1]) * pbinom(0:n1 - 1, n1, xi[2]) * 
               pbinom(0:n1, n1, xi[4]))
    if (s == 4) 
      return(pbinom(0:n1 - 1, n1, xi[1]) * pbinom(0:n1 - 1, n1, xi[2]) * 
               pbinom(0:n1 - 1, n1, xi[3]))
  }
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
  
  pr_Q_y1 <- numeric(n1 + 1)
  beta <- xis + alp * (xis - pis)
  if (abs(beta) < 1e-5) beta <- 0
  q1 <- (beta + alp) / (1 + alp)
  q2 <- beta / (1 + alp)
  if (q1 < 0) q1 <- 0
  if (q2 < 0) q2 <- 0
  if (q1 > 1) q1 <- 1
  if (q2 > 1) q2 <- 1
  
  c <- calc_Pr_Q_x_s(x[1:G], s)
  
  if (beta != 0 & beta != 1){
    for (r in 0:n1) {
      for (x in 0:n1) {
        jdat <- max(0, x + r - n1):min(x, r)
        
        pr_Q_y1[r + 1] <- pr_Q_y1[r + 1] + 
          sum(dbinom(jdat, r, q1) * dbinom(x - jdat, n1 - r, q2)) * c[x + 1]
      }
    }
  }
  
  if (beta == 0) {
    for (r in 0:n1) {
      x <- 0:r
      pr_Q_y1[r + 1] <- (1 + alp) ^ (-r) * 
        sum(choose(r, x) * alp ^ x * c[x + 1])
    }
  }
  
  if (beta == 1) {
    for (r in 0:n1) {
      x <- r:n1
      pr_Q_y1[r + 1] <- (1 + alp) ^ (-(n1 - r)) * 
        sum(choose(n1 - r, x - r) * alp ^ (n1 - x) * c[x + 1])
    }
  }
  
  t <- y0 + ys
  kr <- max(0, t - n):min(t, n)
  h <- numeric(length(kr))
  
  if (length(kr) > 1){
    for(kn in 1:length(kr)) {
      rr <- max(0, kr[kn] - (n - n1)):min(kr[kn], n1)
      h[kn] <- sum(choose(n1, rr) * choose(n - n1, kr[kn] - rr) * 
                     pr_Q_y1[rr + 1])
    }
    
    h <- choose(n, t - kr) * h * exp(delta * kr - max(log(h + 1) + delta * kr))
  } else h <- 1
  
  p0 <- h / sum(h)  
  pval_exact_lower <- sum(p0[kr >= ys])
  pval_exact_upper <- sum(p0[kr <= ys])
  pvals <- c(pval_exact_lower, pval_exact_lower - p0[kr == ys] * 0.5,
             pval_exact_upper, pval_exact_upper - p0[kr == ys] * 0.5)
  names(pvals) <- c("exact_lower", "midp_lower", "exact_upper", "midp-upper")
  
  return(pvals)
}


for (snum in nums) {
  time0 <- Sys.time()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  print(snum)
  print(date())
  
  ssub <- subset(scenario, scenario == snum)
  simn <- ssub$simn
  n <- ssub$n
  rho <- ssub$rho
  tau <- ssub$tau
  G <- ssub$G
  
  # xi = c(xi_1, ..., xi_G, xi_0)
  # pi = c(pi_1, ..., pi_G, pi_0)
  xi <- pi <- numeric(G + 1)
  for (g in 1:G) {
    eval(parse(text = paste0("xi[", g, "] <- ssub$xi_", g)))
    eval(parse(text = paste0("pi[", g, "] <- ssub$pi_", g)))
  }
  xi[G + 1] <- ssub$xi_0
  pi[G + 1] <- ssub$pi_0
  
  
  n1 <- ceiling(n * tau)
  n2 <- n - n1
  stage <- (1:n > n1) + 1
  alpha <- rho / (sqrt(pi * (1 - pi) / (xi * (1 - xi))) - rho)
  
  res <- foreach(sim = 1:simn, .combine = "rbind") %dopar% {
    
    set.seed(snum * 123 + sim * 456)
    Y_i <- rbinom(n * (G + 1), 1, rep(pi, each = n))
    X_i <- rbinom(n * (G + 1), 1, (rep(xi + alpha * (xi - pi), each = n) +
                                     rep(alpha, each = n) * Y_i) / (1 + alpha))
    Y_i <- matrix(Y_i, n, G + 1)
    X_i <- matrix(X_i, n, G + 1)
    
    # interim analysis
    X_IA <- colSums(X_i[stage == 1, 1:G])
    s <- (1:G)[X_IA == max(X_IA)][1]
    
    p00 <- mean(X_i[, s] == 0 & Y_i[, s] == 0)
    p01 <- mean(X_i[, s] == 0 & Y_i[, s] == 1)
    p10 <- mean(X_i[, s] == 1 & Y_i[, s] == 0)
    p11 <- mean(X_i[, s] == 1 & Y_i[, s] == 1)
    
    # ----- estimate binomial probabilities and alpha-----
    # MLE (Takahashi et al., 2022)
    xi_MLE <- X_IA / n1
    xi_MLE[s] <- mean(X_i[, s])
    pi_s_MLE <- mean(Y_i[, s])
    if (xi_MLE[s] == 0) xi_MLE[s] <- 0.001
    if (xi_MLE[s] == 1) xi_MLE[s] <- 0.999
    
    xi_0_MLE <- mean(X_i[, G + 1])
    pi_0_MLE <- mean(Y_i[, G + 1])
    
    alpha_MLE <- (p11 - xi_MLE[s] * pi_s_MLE) / 
      (pi_s_MLE * (1 - pi_s_MLE + xi_MLE[s]) - p11)
    
    if (xi_MLE[s] == pi_s_MLE & p01 == 0) alpha_MLE <- 1e+10
    if (xi_MLE[s] != pi_s_MLE & 
        (xi_MLE[s] %in% c(0, 1) | pi_s_MLE %in% c(0, 1))) alpha_MLE <- 1e+10
    
    
    # CMAE
    cprob <- calc_Pr_Q_x(xi_MLE, G, s)
    Pr_Q <- sum(subset(cprob, g == s)$prob * dbinom(0:n1, n1, xi_MLE[s]))
    
    E <- rep(NA, G)
    for (g0 in 1:G) {
      E[g0] <- sum(0:n1 * dbinom(0:n1, n1, xi_MLE[g0]) * 
                     subset(cprob, g == g0)$prob)
    }
    E[-s] <- E[-s] / (n1 * Pr_Q)
    E[s] <- (E[s] / Pr_Q + n2 * xi_MLE[s]) / n
    xi_CMAE <- xi_MLE - (E - xi_MLE)
    
    xi_CMAE[xi_CMAE < 0] <- 0.001
    xi_CMAE[xi_CMAE > 1] <- 0.999
    
    pr_Q_y1 <- numeric(n1 + 1)
    beta <- xi_MLE[s] + alpha_MLE * (xi_MLE[s] - pi_s_MLE)
    if (abs(beta) < 1e-5) beta <- 0
    q1 <- (beta + alpha_MLE) / (1 + alpha_MLE)
    q2 <- beta / (1 + alpha_MLE)
    if (q1 < 0) q1 <- 0
    if (q2 < 0) q2 <- 0
    if (q1 > 1) q1 <- 1
    if (q2 > 1) q2 <- 1
    
    c <- subset(cprob, g == s)$prob
    if (beta != 0 & beta != 1){
      for (r in 0:n1) {
        for (x in 0:n1) {
          jdat <- max(0, x + r - n1):min(x, r)
          
          pr_Q_y1[r + 1] <- pr_Q_y1[r + 1] + 
            sum(dbinom(jdat, r, q1) * dbinom(x - jdat, n1 - r, q2)) * c[x + 1]
        }
      }
    }
    
    if (beta == 0) {
      for (r in 0:n1) {
        x <- 0:r
        pr_Q_y1[r + 1] <- (1 + alpha_MLE) ^ (-r) * 
          sum(choose(r, x) * alpha_MLE ^ x * c[x + 1])
      }
    }
    
    if (beta == 1) {
      for (r in 0:n1) {
        x <- r:n1
        pr_Q_y1[r + 1] <- (1 + alpha_MLE) ^ (-(n1 - r)) * 
          sum(choose(n1 - r, x - r) * alpha_MLE ^ (n1 - x) * c[x + 1])
      }
    }
    
    E <- (sum((0:n1) * dbinom(0:n1, n1, pi_s_MLE) * pr_Q_y1) / Pr_Q + 
            n2 * pi_s_MLE) / n
    pi_s_CMAE <- pi_s_MLE - (E - pi_s_MLE)
    
    alpha_CMAE <- optimize(
      function (x) lik_alp(xi_CMAE[s], pi_s_CMAE, x),
      interval = alp_val(xi_CMAE[s], pi_s_CMAE))$minimum
    
    # UMVCUE
    xs <- sum(X_i[, s])
    ys <- sum(Y_i[, s])
    zs <- sum(X_i[, s] == Y_i[, s] & X_i[, s] == 1)
    xr <- max(X_IA[(1:G) < s] + 1, X_IA[(1:G) > s], xs - n2):min(n1, xs)
    
    q <- choose(n1, xr) * choose(n2, xs - xr)
    xi_UMVCUE <- xi_MLE
    xi_UMVCUE[s] <- (xs - sum(xr * q) / sum(q)) / n2
    if (xi_UMVCUE[s] <= 0) xi_UMVCUE[s] <- 0.001
    if (xi_UMVCUE[s] >= 1) xi_UMVCUE[s] <- 0.999
    
    c1 <- yy <- c()
    for (x1 in xr) {
      x2 <- xs - x1
      yr <- max(0, ys - n2, x2 + ys - zs - n2):min(n1, ys, n1 + zs - x1)
      yr <- yr[pmin(yr, x1) + pmin(ys - yr, x2) >= zs]
      yy <- c(yy, yr)
      for (y1 in yr) {
        y2 <- ys - y1
        z1 <- (zs - min(x2, y2)):min(x1, y1, zs)
        c0 <- sum(choose(y1, z1) * 
                    choose(n1 - y1, x1 - z1) * 
                    choose(y2, zs - z1) * 
                    choose(n2 - y2, x2 - zs + z1))
        c1 <- c(c1, c0)
      }
    }
    
    c2 <- choose(n2, ys - yy) * choose(n1, yy)
    
    if (sum(is.infinite(c1 * c2)) > 0) {
      dmean <- mean(log10(c1) + log10(c2))
      c1 <- c1 * 10 ^ (-dmean)
      c2 <- c2 * 10 ^ (-dmean)
    }
    
    pi_s_UMVCUE <- (ys - sum(yy * c1 * c2) / sum(c1 * c2)) / n2
    
    if (pi_s_UMVCUE <= 0) pi_s_UMVCUE <- 0.001
    if (pi_s_UMVCUE >= 1) pi_s_UMVCUE <- 0.999
    
    alpha_UMVCUE <- optimize(
      function (x) lik_alp(xi_UMVCUE[s], pi_s_UMVCUE, x),
      interval = alp_val(xi_UMVCUE[s], pi_s_UMVCUE))$minimum
    
    # ----- calculate confidence limits -----
    # MLE
    if (!is.na(alpha_MLE)) {
      val <- c(xi_MLE[1:G], pi_s_MLE, alpha_MLE, 
               sum(Y_i[, G + 1]), sum(Y_i[, s]))
      
      pval_MLE <- calc_pval(val, delta = 0)
      
      pval1 <- calc_pval(val, delta = -10)
      pval2 <- calc_pval(val, delta = 10)
      
      logOR <- log(pi_s_MLE / (1 - pi_s_MLE)) - 
        log(pi_0_MLE / (1 - pi_0_MLE))
      ase <- sqrt((1 / pi_s_MLE + 1 / (1 - pi_s_MLE) + 
                     1 / pi_0_MLE + 1 / (1 - pi_0_MLE)) / n)
      
      if (pval1[1] < 0.025 & pval2[1] > 0.025) {
        LL_exact_MLE <- 
          try(uniroot(function (d) calc_pval(val, d)[1] - 0.025, 
                      c(logOR - 3 * ase, logOR - ase))$root, silent = TRUE)
        if (class(LL_exact_MLE) == "try-error") {
          LL_exact_MLE <- uniroot(function (d) calc_pval(val, d)[1] - 0.025, 
                                  c(-10, 10))$root
        }
      } else {
        if (pval1[1] > 0.025) LL_exact_MLE <- -10
        if (pval2[1] < 0.025) LL_exact_MLE <- 10
      }
      
      if (pval1[2] < 0.025 & pval2[2] > 0.025) {
        LL_midp_MLE <- 
          try(uniroot(function (d) calc_pval(val, d)[2] - 0.025, 
                      c(LL_exact_MLE, logOR - ase))$root, silent = TRUE)
        if (class(LL_midp_MLE) == "try-error") {
          LL_midp_MLE <- uniroot(function (d) calc_pval(val, d)[2] - 0.025, 
                                 c(-10, 10))$root
        }
      } else {
        if (pval1[2] > 0.025) LL_midp_MLE <- -10
        if (pval2[2] < 0.025) LL_midp_MLE <- 10
      }
      
      if (pval1[3] > 0.025 & pval2[3] < 0.025) {
        UL_exact_MLE <- 
          try(uniroot(function (d) calc_pval(val, d)[3] - 0.025, 
                      c(logOR + ase, logOR + 3 * ase))$root, silent = TRUE)
        if (class(UL_exact_MLE) == "try-error") {
          UL_exact_MLE <- uniroot(function (d) calc_pval(val, d)[3] - 0.025, 
                                  c(-10, 10))$root
        }
      } else {
        if (pval1[3] < 0.025) UL_exact_MLE <- -10
        if (pval2[3] > 0.025) UL_exact_MLE <- 10
      }
      
      if (pval1[4] > 0.025 & pval2[4] < 0.025) {
        UL_midp_MLE <- 
          try(uniroot(function (d) calc_pval(val, d)[4] - 0.025, 
                      c(logOR + ase, UL_exact_MLE))$root, silent = TRUE)
        if (class(UL_midp_MLE) == "try-error") {
          UL_midp_MLE <- uniroot(function (d) calc_pval(val, d)[4] - 0.025, 
                                 c(-10, 10))$root
        }
      } else {
        if (pval1[4] < 0.025) UL_midp_MLE <- -10
        if (pval2[4] > 0.025) UL_midp_MLE <- 10
      }
      
    } else {
      pval_MLE <- c(NA, NA)
      LL_exact_MLE <- LL_midp_MLE <- UL_exact_MLE <- UL_midp_MLE <- NA
    }
    
    # CMAE
    if (!is.na(pi_s_CMAE) & !is.na(alpha_CMAE)) {
      val <- c(xi_CMAE, pi_s_CMAE, alpha_CMAE, 
               sum(Y_i[, G + 1]), sum(Y_i[, s]))
      
      pval_CMAE <- calc_pval(val, delta = 0)
      
      pval1 <- calc_pval(val, delta = -10)
      pval2 <- calc_pval(val, delta = 10)
      
      logOR <- log(pi_s_CMAE / (1 - pi_s_CMAE)) - 
        log(pi_0_MLE / (1 - pi_0_MLE))
      ase <- sqrt((1 / pi_s_CMAE + 1 / (1 - pi_s_CMAE) + 
                     1 / pi_0_MLE + 1 / (1 - pi_0_MLE)) / n)
      
      if (pval1[1] < 0.025 & pval2[1] > 0.025) {
        LL_exact_CMAE <- 
          try(uniroot(function (d) calc_pval(val, d)[1] - 0.025, 
                      c(logOR - 3 * ase, logOR - ase))$root, silent = TRUE)
        if (class(LL_exact_CMAE) == "try-error") {
          LL_exact_CMAE <- uniroot(function (d) calc_pval(val, d)[1] - 0.025, 
                                   c(-10, 10))$root
        }
      } else {
        if (pval1[1] > 0.025) LL_exact_CMAE <- -10
        if (pval2[1] < 0.025) LL_exact_CMAE <- 10
      }
      
      if (pval1[2] < 0.025 & pval2[2] > 0.025) {
        LL_midp_CMAE <- 
          try(uniroot(function (d) calc_pval(val, d)[2] - 0.025, 
                      c(LL_exact_CMAE, logOR - ase))$root, silent = TRUE)
        if (class(LL_midp_CMAE) == "try-error") {
          LL_midp_CMAE <- uniroot(function (d) calc_pval(val, d)[2] - 0.025, 
                                  c(-10, 10))$root
        }
      } else {
        if (pval1[2] > 0.025) LL_midp_CMAE <- -10
        if (pval2[2] < 0.025) LL_midp_CMAE <- 10
      }
      
      if (pval1[3] > 0.025 & pval2[3] < 0.025) {
        UL_exact_CMAE <- 
          try(uniroot(function (d) calc_pval(val, d)[3] - 0.025, 
                      c(logOR + ase, logOR + 3 * ase))$root, silent = TRUE)
        if (class(UL_exact_CMAE) == "try-error") {
          UL_exact_CMAE <- uniroot(function (d) calc_pval(val, d)[3] - 0.025, 
                                   c(-10, 10))$root
        }
      } else {
        if (pval1[3] < 0.025) UL_exact_CMAE <- -10
        if (pval2[3] > 0.025) UL_exact_CMAE <- 10
      }
      
      if (pval1[4] > 0.025 & pval2[4] < 0.025) {
        UL_midp_CMAE <- 
          try(uniroot(function (d) calc_pval(val, d)[4] - 0.025, 
                      c(logOR + ase, UL_exact_CMAE))$root, silent = TRUE)
        if (class(UL_midp_CMAE) == "try-error") {
          UL_midp_CMAE <- uniroot(function (d) calc_pval(val, d)[4] - 0.025, 
                                  c(-10, 10))$root
        }
      } else {
        if (pval1[4] < 0.025) UL_midp_CMAE <- -10
        if (pval2[4] > 0.025) UL_midp_CMAE <- 10
      }
      
    } else {
      pval_CMAE <- c(NA, NA)
      LL_exact_CMAE <- LL_midp_CMAE <- UL_exact_CMAE <- UL_midp_CMAE <- NA
    }
    
    # UMVCUE
    if (!is.na(pi_s_UMVCUE) & !is.na(alpha_UMVCUE)) {
      val <- c(xi_UMVCUE, pi_s_UMVCUE, alpha_UMVCUE, 
               sum(Y_i[, G + 1]), sum(Y_i[, s]))
      
      pval_UMVCUE <- calc_pval(val, delta = 0)
      
      pval1 <- calc_pval(val, delta = -10)
      pval2 <- calc_pval(val, delta = 10)
      
      logOR <- log(pi_s_UMVCUE / (1 - pi_s_UMVCUE)) - 
        log(pi_0_MLE / (1 - pi_0_MLE))
      ase <- sqrt((1 / pi_s_UMVCUE + 1 / (1 - pi_s_UMVCUE) + 
                     1 / pi_0_MLE + 1 / (1 - pi_0_MLE)) / n)
      
      if (pval1[1] < 0.025 & pval2[1] > 0.025) {
        LL_exact_UMVCUE <- 
          try(uniroot(function (d) calc_pval(val, d)[1] - 0.025, 
                      c(logOR - 3 * ase, logOR - ase))$root, silent = TRUE)
        if (class(LL_exact_UMVCUE) == "try-error") {
          LL_exact_UMVCUE <- uniroot(function (d) calc_pval(val, d)[1] - 0.025, 
                                     c(-10, 10))$root
        }
      } else {
        if (pval1[1] > 0.025) LL_exact_UMVCUE <- -10
        if (pval2[1] < 0.025) LL_exact_UMVCUE <- 10
      }
      
      if (pval1[2] < 0.025 & pval2[2] > 0.025) {
        LL_midp_UMVCUE <- 
          try(uniroot(function (d) calc_pval(val, d)[2] - 0.025, 
                      c(LL_exact_UMVCUE, logOR - ase))$root, silent = TRUE)
        if (class(LL_midp_UMVCUE) == "try-error") {
          LL_midp_UMVCUE <- uniroot(function (d) calc_pval(val, d)[2] - 0.025, 
                                    c(-10, 10))$root
        }
      } else {
        if (pval1[2] > 0.025) LL_midp_UMVCUE <- -10
        if (pval2[2] < 0.025) LL_midp_UMVCUE <- 10
      }
      
      if (pval1[3] > 0.025 & pval2[3] < 0.025) {
        UL_exact_UMVCUE <- 
          try(uniroot(function (d) calc_pval(val, d)[3] - 0.025, 
                      c(logOR + ase, logOR + 3 * ase))$root, silent = TRUE)
        if (class(UL_exact_UMVCUE) == "try-error") {
          UL_exact_UMVCUE <- uniroot(function (d) calc_pval(val, d)[3] - 0.025, 
                                     c(-10, 10))$root
        }
      } else {
        if (pval1[3] < 0.025) UL_exact_UMVCUE <- -10
        if (pval2[3] > 0.025) UL_exact_UMVCUE <- 10
      }
      
      if (pval1[4] > 0.025 & pval2[4] < 0.025) {
        UL_midp_UMVCUE <- 
          try(uniroot(function (d) calc_pval(val, d)[4] - 0.025, 
                      c(logOR + ase, UL_exact_UMVCUE))$root, silent = TRUE)
        if (class(UL_midp_UMVCUE) == "try-error") {
          UL_midp_UMVCUE <- uniroot(function (d) calc_pval(val, d)[4] - 0.025, 
                                    c(-10, 10))$root
        }
      } else {
        if (pval1[4] < 0.025) UL_midp_UMVCUE <- -10
        if (pval2[4] > 0.025) UL_midp_UMVCUE <- 10
      }
      
    } else {
      pval_UMVCUE <- c(NA, NA)
      LL_exact_UMVCUE <- LL_midp_UMVCUE <- 
        UL_exact_UMVCUE <- UL_midp_UMVCUE <- NA
    }
    
    res0 <-
      data.frame(sim, s, xi_0_MLE, t(xi_MLE), xi_MLE[s], 
                 t(xi_CMAE), xi_CMAE[s], t(xi_UMVCUE), xi_UMVCUE[s],
                 pi_0_MLE, pi_s_MLE, pi_s_CMAE, pi_s_UMVCUE,
                 alpha_MLE, alpha_CMAE, alpha_UMVCUE,
                 t(pval_MLE[1:2]), t(pval_CMAE[1:2]), t(pval_UMVCUE[1:2]), 
                 xi[s], pi[s], 
                 log(pi[s] / (1 - pi[s])) - log(pi[G + 1] / (1 - pi[G + 1])),
                 LL_exact_MLE, LL_midp_MLE, UL_exact_MLE, UL_midp_MLE,
                 LL_exact_CMAE, LL_midp_CMAE, UL_exact_CMAE, UL_midp_CMAE,
                 LL_exact_UMVCUE, LL_midp_UMVCUE,
                 UL_exact_UMVCUE, UL_midp_UMVCUE)
    
    colnames(res0) <- 
      c("sim", "s","xi_0_hat", 
        paste0("xi_", 1:G, "_MLE"), "xi_s_MLE",
        paste0("xi_", 1:G, "_CMAE"), "xi_s_CMAE",
        paste0("xi_", 1:G, "_UMVCUE"), "xi_s_UMVCUE",
        "pi_0_MLE", "pi_s_MLE", "pi_s_CMAE", "pi_s_UMVCUE",
        "alpha_MLE", "alpha_CMAE", "alpha_UMVCUE",
        "pval_exact_MLE", "pval_midp_MLE",
        "pval_exact_CMAE", "pval_midp_CMAE", 
        "pval_exact_UMVCUE", "pval_midp_UMVCUE",
        "xi_s", "pi_s", "delta_s",
        "LL_exact_MLE", "LL_midp_MLE", "UL_exact_MLE", "UL_midp_MLE",
        "LL_exact_CMAE", "LL_midp_CMAE", "UL_exact_CMAE", "UL_midp_CMAE",
        "LL_exact_UMVCUE", "LL_midp_UMVCUE",
        "UL_exact_UMVCUE", "UL_midp_UMVCUE")
    
    res0
  }
  
  parms <- as.data.frame(t(c(n, tau, rho, G,
                             xi[G + 1], xi[1:G], pi[G + 1], pi[1:G], alpha)))
  colnames(parms) <- c("n", "tau", "rho", "G", "xi_0",
                       paste0("xi_", 1:G), "pi_0", paste0("pi_", 1:G), 
                       paste0("alpha_", 0:G))
  
  
  res1 <- 
    data.frame(E_xi_s = mean(res$xi_s, na.rm = TRUE),
               E_pi_s = mean(res$pi_s, na.rm = TRUE),
               E_xi_s_MLE = mean(res$xi_s_MLE, na.rm = TRUE),
               E_xi_s_CMAE = mean(res$xi_s_CMAE, na.rm = TRUE),
               E_xi_s_UMVCUE = mean(res$xi_s_UMVCUE, na.rm = TRUE),
               E_pi_s_MLE = mean(res$pi_s_MLE, na.rm = TRUE),
               E_pi_s_CMAE = mean(res$pi_s_CMAE, na.rm = TRUE),
               E_pi_s_UMVCUE = mean(res$pi_s_UMVCUE, na.rm = TRUE),
               RMSE_xi_s_MLE = 
                 sqrt(mean((res$xi_s_MLE - res$xi_s) ^ 2, na.rm = TRUE)),
               RMSE_xi_s_CMAE = 
                 sqrt(mean((res$xi_s_CMAE - res$xi_s) ^ 2, na.rm = TRUE)),
               RMSE_xi_s_UMVCUE = 
                 sqrt(mean((res$xi_s_UMVCUE - res$xi_s) ^ 2, na.rm = TRUE)),
               RMSE_pi_s_MLE = 
                 sqrt(mean((res$pi_s_MLE - res$pi_s) ^ 2, na.rm = TRUE)),
               RMSE_pi_s_CMAE = 
                 sqrt(mean((res$pi_s_CMAE - res$pi_s) ^ 2, na.rm = TRUE)),
               RMSE_pi_s_UMVCUE = 
                 sqrt(mean((res$pi_s_UMVCUE - res$pi_s) ^ 2, na.rm = TRUE)),
               CP_exact_MLE = mean(res$LL_exact_MLE <= res$delta_s &
                                     res$UL_exact_MLE >= res$delta_s,
                                   na.rm = TRUE) * 100,
               CP_exact_CMAE = mean(res$LL_exact_CMAE <= res$delta_s &
                                      res$UL_exact_CMAE >= res$delta_s,
                                    na.rm = TRUE) * 100,
               CP_exact_UMVCUE = mean(res$LL_exact_UMVCUE <= res$delta_s &
                                        res$UL_exact_UMVCUE >= res$delta_s,
                                      na.rm = TRUE) * 100,
               CP_midp_MLE = mean(res$LL_midp_MLE <= res$delta_s &
                                    res$UL_midp_MLE >= res$delta_s,
                                  na.rm = TRUE) * 100,
               CP_midp_CMAE = mean(res$LL_midp_CMAE <= res$delta_s &
                                     res$UL_midp_CMAE >= res$delta_s,
                                   na.rm = TRUE) * 100,
               CP_midp_UMVCUE = mean(res$LL_midp_UMVCUE <= res$delta_s &
                                       res$UL_midp_UMVCUE >= res$delta_s,
                                     na.rm = TRUE) * 100,
               size_exact_MLE = mean(res$pval_exact_MLE < 0.025,
                                     na.rm = TRUE) * 100,
               size_exact_CMAE = mean(res$pval_exact_CMAE < 0.025,
                                      na.rm = TRUE) * 100,
               size_exact_UMVCUE = mean(res$pval_exact_UMVCUE < 0.025,
                                        na.rm = TRUE) * 100,
               size_midp_MLE = mean(res$pval_midp_MLE < 0.025,
                                    na.rm = TRUE) * 100,
               size_midp_CMAE = mean(res$pval_midp_CMAE < 0.025,
                                     na.rm = TRUE) * 100,
               size_midp_UMVCUE = mean(res$pval_midp_UMVCUE < 0.025,
                                       na.rm = TRUE) * 100
    )
  
  write.csv(cbind(parms, res),  paste0("all_", snum,".csv"))
  write.csv(cbind(parms, res1), paste0("Summary_", snum,".csv"))
  
  stopCluster(cl)
  
  print(paste0("snum = ", snum, ": completed"))
  print(Sys.time() - time0)
}

