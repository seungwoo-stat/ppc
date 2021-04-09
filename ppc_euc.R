## source code for PPC on the Euclidean space R^p.
## one-step procedure is made into a function

library(magrittr)
library(pracma)

## this function works only for 2d Euclidean space
normal.e.plot <- function(mu, sigma, conf.level = 0.95,...){
  k <- sqrt(qchisq(conf.level, df = 2))
  theta <- seq(0,2*pi, by = pi/360)
  e <- matrix(nrow = length(theta),ncol=2)
  for(t in seq_along(theta)){
    e[t,] <- mu + k*sigma*c(cos(theta[t]),sin(theta[t]))
  }
  points(mu[1], mu[2], pch = "+", ...)
  lines(e, ...)
}

## for smoothing splines
band.Q <- function(knots){
  k <- sort(knots)
  n <- length(k)
  h <- k[2:n] - k[1:(n-1)]
  Q <- matrix(0, nrow = n, ncol = n-2)
  for(j in 1:(n-2)){
    Q[j,j] <- 1/h[j]
    Q[j+1,j] <- -1/h[j] -1/h[j+1]
    Q[j+2,j] <- 1/h[j+1]
  }
  return(Q)
}

## ((n-2) x (n-2)) matrix, n: number of knots
## n > 2
band.R <- function(knots){
  k <- sort(knots)
  n <- length(k)
  h <- k[2:n] - k[1:(n-1)]
  R <- matrix(0, nrow = n-2, ncol = n-2)
  for(j in 1:(n-2)){
    if(j>=1) R[j-1,j] <- h[j]/6
    R[j,j] <- (h[j] + h[j+1])/3
    if(j<=n-3) R[j+1,j] <- h[j+1]/6
  }
  return(R)
}

mat.K <- function(x) band.Q(x) %*% solve(band.R(x), t(band.Q(x)))

## for 2D
ppc.iter.e <- function(data, w, a, v, sigma, sigma.max = Inf, lambda){
  n <- nrow(data)
  y1 <- data[,1]; y2 <- data[,2]
  # D1 <- diag(colSums(w)/sigma^2)
  # D2 <- diag(colSums(w)/sigma^2)
  d1 <- colSums(w)/sigma^2
  d2 <- colSums(w)/sigma^2
  y1.bar <- t(w) %*% y1/sigma^2
  y2.bar <- t(w) %*% y2/sigma^2
  # K <- mat.K(a)
  fit <- matrix(ncol = 2, nrow = n)
  # fit[,1] <- solve(D1 + 2*lambda* K, y1.bar)
  # fit[,2] <- solve(D2 + 2*lambda* K, y2.bar)
  
  ## REINSCH ALGORITHM
  k <- sort(a)
  h <- k[2:n] - k[1:(n-1)]
  QtDinvY1 <- sapply(2:(n-1), function(i) (y1.bar[i+1]/d1[i+1] - y1.bar[i]/d1[i])/h[i] 
                     - (y1.bar[i]/d1[i]-y1.bar[i-1]/d1[i-1])/h[i-1])
  QtDinvY2 <- sapply(2:(n-1), function(i) (y2.bar[i+1]/d2[i+1] - y2.bar[i]/d2[i])/h[i] 
                     - (y2.bar[i]/d2[i]-y2.bar[i-1]/d2[i-1])/h[i-1]) 
  Q <- band.Q(k)
  R <- band.R(k)
  tryCatch(
    {
      U1 <- pracma::lu(R + 2*lambda*t(Q) %*% (Q/d1))
      U2 <- pracma::lu(R + 2*lambda*t(Q) %*% (Q/d2))
      gamma1 <- backsolve(U1$U, forwardsolve(U1$L, QtDinvY1))
      gamma2 <- backsolve(U2$U, forwardsolve(U2$L,QtDinvY2))
      fit[,1] <- (y1.bar - 2*lambda*Q%*%gamma1)/d1
      fit[,2] <- (y2.bar - 2*lambda*Q%*%gamma2)/d2
    },
    error = function(cond){
      K <- Q %*% solve(R, t(Q))
      fit[,1] <- solve(diag(d1,n) + 2*lambda* K, y1.bar)
      fit[,2] <- solve(diag(d2,n) + 2*lambda* K, y2.bar)
      message(cond)
    }
  )
  
  return(fit)
}

## for 3D
ppc.iter.e.3d <- function(data, w, a, v, sigma, sigma.max = Inf, lambda){
  n <- nrow(data)
  y1 <- data[,1]; y2 <- data[,2]; y3 <- data[,3]
  # D1 <- diag(colSums(w)/sigma^2)
  # D2 <- diag(colSums(w)/sigma^2)
  d1 <- colSums(w)/sigma^2
  d2 <- d1
  d3 <- d1
  y1.bar <- t(w) %*% y1/sigma^2
  y2.bar <- t(w) %*% y2/sigma^2
  y3.bar <- t(w) %*% y3/sigma^2
  # K <- mat.K(a)
  fit <- matrix(ncol = 3, nrow = n)
  # fit[,1] <- solve(D1 + 2*lambda* K, y1.bar)
  # fit[,2] <- solve(D2 + 2*lambda* K, y2.bar)
  
  ## REINSCH ALGORITHM
  k <- sort(a)
  h <- k[2:n] - k[1:(n-1)]
  QtDinvY1 <- sapply(2:(n-1), function(i) (y1.bar[i+1]/d1[i+1] - y1.bar[i]/d1[i])/h[i] 
                     - (y1.bar[i]/d1[i]-y1.bar[i-1]/d1[i-1])/h[i-1])
  QtDinvY2 <- sapply(2:(n-1), function(i) (y2.bar[i+1]/d2[i+1] - y2.bar[i]/d2[i])/h[i] 
                     - (y2.bar[i]/d2[i]-y2.bar[i-1]/d2[i-1])/h[i-1]) 
  QtDinvY3 <- sapply(2:(n-1), function(i) (y3.bar[i+1]/d2[i+1] - y3.bar[i]/d2[i])/h[i] 
                     - (y3.bar[i]/d2[i]-y3.bar[i-1]/d2[i-1])/h[i-1]) 
  Q <- band.Q(k)
  R <- band.R(k)
  tryCatch(
    {
      U1 <- pracma::lu(R + 2*lambda*t(Q) %*% (Q/d1))
      U2 <- pracma::lu(R + 2*lambda*t(Q) %*% (Q/d2))
      U3 <- pracma::lu(R + 2*lambda*t(Q) %*% (Q/d3))
      gamma1 <- backsolve(U1$U, forwardsolve(U1$L,QtDinvY1))
      gamma2 <- backsolve(U2$U, forwardsolve(U2$L,QtDinvY2))
      gamma3 <- backsolve(U3$U, forwardsolve(U3$L,QtDinvY3))
      fit[,1] <- (y1.bar - 2*lambda*Q%*%gamma1)/d1
      fit[,2] <- (y2.bar - 2*lambda*Q%*%gamma2)/d2
      fit[,3] <- (y3.bar - 2*lambda*Q%*%gamma3)/d3
    },
    error = function(cond){
      K <- Q %*% solve(R, t(Q))
      fit[,1] <- solve(diag(d1,n) + 2*lambda* K, y1.bar)
      fit[,2] <- solve(diag(d2,n) + 2*lambda* K, y2.bar)
      fit[,3] <- solve(diag(d3,n) + 2*lambda* K, y3.bar)
      message(cond)
    }
  )
  
  return(fit)
}
