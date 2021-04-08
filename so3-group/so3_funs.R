source("../generics.R")
source("../ppc_euc.R")

library(expm)
library(rgl)
library(sphereplot)
library(MASS)
library(progress)

IDEN <- diag(1,ncol=3,nrow=3)

#' define class `so3`
#' 
#' @param l A list or matrix
#' @return A list or matrix of `so3` class
new_so3 <- function(l){
  if(is.matrix(l)){
    stopifnot(all(abs(l %*% t(l) - IDEN) < 1e-6))
    stopifnot(abs(det(l) - 1) < 1e-6)
    return(structure(l, class = "so3"))
  }
  stopifnot(all(sapply(1:length(l), function(i) all(abs(l[[i]] %*% t(l[[i]]) - IDEN) < 1e-6))))
  stopifnot(all(sapply(1:length(l), function(i) all(abs(det(l[[i]]) - 1) < 1e-6))))
  structure(l, class = "so3")
}

#' Constructor of class \code{so3} 
#' 
#' @param x A list of so3 elements or a matrix of single so3 element.
#' @return `so3` class object.
#' @export
#' @examples
#' R <- matrix(c(0,0,1, 0,-1,0, 1,0,0),ncol=3)
#' l <- list(IDEN, R)
#' so3(l)
#' so3(R)
so3 <- function(x){
  if(is.matrix(x)) return(new_so3(x))
  new_so3(x)
}

# `[[.so3` <- function(x, i, ...) {
#   out <- unclass(x)
#   out <- out[[i,...]]
#   return(so3(out))
# }

## plot methods for so3 --------------------------------------------------------

plot.mat_so3 <- function(mat){
  mat[1,] <- mat[1,]/2
  # mat[2,] <- mat[2,]
  mat[3,] <- mat[3,]*1.5
  c1 <- mat[1,]
  c2 <- mat[2,]
  c3 <- mat[3,]
  l1 <- rbind(c(0,0,0),c1,c1+c2,c2,c(0,0,0),
              c3,c3+c1,c3+c2+c1,c3+c2,c3)
  l2 <- rbind(c1,c1+c3)
  l3 <- rbind(c2,c2+c3)
  l4 <- rbind(c1+c2,c1+c2+c3)
  
  trans <- (c1+c2+c3)/2
  l1 <- t(t(l1)-trans)
  l2 <- t(t(l2)-trans)
  l3 <- t(t(l3)-trans)
  l4 <- t(t(l4)-trans)
  
  plot3d(t(t(rbind(c(0,0,0),mat))-trans),
         xlim = c(-1.5,1.5),ylim = c(-1.5,1.5),zlim = c(-1.5,1.5), 
         xlab = "", ylab= "", zlab = "", axes = FALSE)
  lines3d(l1, lwd = 3)
  lines3d(l2, lwd = 3)
  lines3d(l3, lwd = 3)
  lines3d(l4, lwd = 3)
  q1 <- rbind(c(0,0,0),c1,c1+c2,c2)
  q1 <- t(t(q1)-trans)
  # q2 <- rbind(c(0,0,0),c1,c1+c3,c3)
  # q2 <- t(t(q2)-trans)
  # q3 <- rbind(c(0,0,0),c2,c2+c3,c3)
  # q3 <- t(t(q3)-trans)
  quads3d(q1, col="lightgray", alpha = 0.4)
  # quads3d(q2, col="lightgray", alpha = 0.6)
  # quads3d(q3, col="lightgray", alpha = 0.9)
  # lines3d(t(t(rbind(c(0,0,0),c1+c2+c3))-trans), col = "red", lwd = 3)
}

#' Method of \code{plot} generic for \code{so3} class
#' 
#' @param data object of class \code{so3}
#' @param red additional data for drawing
#' @param black additional data for drawing
#' @param gray additional data for drawing
#' @param threeD 0 when drawing a 2d MDS plot, if bigger, draws `threeD` number of 3d plots. 
#' @return No return, draws a plot.
#' @export
#' @examples
#' 
plot.so3 <- function(data, red = NULL, black = NULL, gray = NULL, threeD = 0, ...){
  if(threeD >= 1L){
    n <- length(data)
    n_plot <- min(n, threeD)
    mfrow3d(ifelse(is.null(red),1,2),n_plot, sharedMouse = TRUE)
    for(i in round(seq(1,n,length.out=n_plot))){
      plot.mat_so3(data[[i]])
    }
    if(!is.null(red)){
      for(i in round(seq(1,min(n,length(red)),length.out=n_plot))){
        plot.mat_so3(red[[i]])
      }
    }
  }else{
    n <- length(data)
    data_merge <- data
    n_red <- 0; n_gray <- 0; n_black <- 0
    n_total <- length(data)
    if(!is.null(red)){
      data_merge <- c(data_merge, red)
      n_red <- (n_total+1):(n_total+length(red))
      n_total <- n_total + length(n_red)
    }
    if(!is.null(gray)){
      data_merge <- c(data_merge, gray)
      n_gray <- (n_total+1):(n_total+length(gray))
      n_total <- n_total + length(n_gray)
    }
    if(!is.null(black)){
      data_merge <- c(data_merge, black)
      n_black <- (n_total+1):(n_total+length(black))
      n_total <- n_total + length(n_black)
    }
    d <- matrix(0, ncol=n_total, nrow= n_total)
    for(i in 1:(n_total-1)){
      for(j in (i+1):n_total){
        d[i,j] <- DIST.so3(data_merge[[i]], data_merge[[j]])
        # if(is.nan(d[i,j])) print(paste(i, j))
      }
    }
    d <- d + t(d)
    d <- d + 1e-10
    diag(d) <- 0
    fit <- isoMDS(d, k=2)$points
    plot(fit[1:n,1],fit[1:n,2], xlab = "MDS 1", ylab = "MDS 2", pch = 20, ...)
    if(!is.null(red)){
      lines(fit[n_red,], col = "red", lwd=4)
    }
    if(!is.null(gray)){
      lines(fit[n_gray,], col = "gray", lwd=4)
    }
    if(!is.null(black)){
      lines(fit[n_black,], col = "black", lwd=4)
    }
  }
}

## -----------------------------------------------------------------------------

#' Frobenius norm of the matrix
#' 
#' @param mat matrix
#' @return Frobenius matrix
#' @export
#' @examples
#' frob.norm(IDEN) # = sqrt(3)
frob.norm <- function(mat){
  sqrt(sum(mat^2))
}

#' x, y vector
cross.vec <- function(x, y){
  c(x[2]*y[3]-x[3]*y[2],
    -x[1]*y[3]+x[3]*y[1],
    x[1]*y[2]-x[2]*y[1])
}

## axis: axis
## theta: rotation angle
rot.mat <- function(axis,theta){
  axis <- axis/sqrt(sum(axis^2))
  x <- axis[1]
  y <- axis[2]
  z <- axis[3]
  R <- matrix(c(cos(theta) + x^2*(1-cos(theta)), x*y*(1-cos(theta)) - z*sin(theta), x*z*(1-cos(theta))+y*sin(theta),
                y*x*(1-cos(theta))+z*sin(theta), cos(theta) + y^2*(1-cos(theta)), y*z*(1-cos(theta))-x*sin(theta),
                z*x*(1-cos(theta))-y*sin(theta), z*y*(1-cos(theta))+x*sin(theta), cos(theta)+z^2*(1-cos(theta))), 
              byrow = T, ncol = 3)
  return(R)
}

#' Geodesic distance between two points. 
#' 
#' Vectorized operations are not supported.
#' @param A 3*3 so(3) matrix
#' @param B 3*3 so(3) matrix
#' @return The geodesic distance between \code{A} and \code{B}. 
#' @export
#' @examples
#' DIST.so3(IDEN, H)
DIST.so3 <- function(A, B){
  trace <- sum(diag(t(A) %*% B))
  if(trace <= -1) return(sqrt(2)*pi)
  if(trace >= 3) return(0)
  return( sqrt(2) * acos((trace - 1)/2) )
}

#' Matrix exponential of so(3) Lie group
#' 
#' Use Rodrigues's formula (1840)
#' @param A 3*3 skew-symmetric matrix
#' @return exp(A) in SO(3)
exp.so3 <- function(A){
  theta <- sqrt(A[1,2]^2 + A[1,3]^2 + A[2,3]^2)
  if(abs(theta) < 1e-7) return(IDEN)
  IDEN + sin(theta)/theta * A + (1-cos(theta))/theta^2 * A%*%A
}

#' Matrix logarithm of SO(3) group
#' 
#' Since matrix log is multi-valued function, we restrict tr(A) to [-1,3].
#' @param A 3*3 orthogonal matrix
#' @return log(A) (skew-symmetric)
log.so3 <- function(A){
  if((sum(diag(A))-1)/2 >= 1) return(matrix(0,ncol=3,nrow=3))
  theta <- acos((sum(diag(A))-1)/2)
  theta/2/sin(theta)*(A - t(A))
  # logm(A)
}

#' Returns a logarithmic map
#' 
#' @param A SO(3) matrix.
#' @param B SO(3) matrix. 
#' @return Log_A(B), t(A) %*% Log_A(B) is skew-symmetric
#' @export
#' @examples
#' LOG(ll[[1]],ll[[2]])
LOG.so3 <- function(A, B){
  A %*% log.so3(t(A) %*% B)
}


#' Returns an exponential map
#' 
#' @param A SO(3) matrix
#' @param H A tangent vector on T_A(SO(3))
#' @return Exp_A(H). 
#' @export
#' @examples
#' EXP(ll[[1]],LOG(ll[[1]],ll[[2]]))
EXP.so3 <- function(A,H){
  A %*% exp.so3(t(A) %*% H)
}

#' Parallel transport a vector on the tangent space via geodesic
#' 
#' @param start SO(3) matrix.
#' @param end SO(3) matrix. 
#' @param vec T_start(SO(3)) matrix
#' @return Parallel transport \code{vec} from \code{start} to \code{end} via a geodesic. 
#' @export
#' @examples
#' PARALLEL(ll[[1]],ll[[2]],H)
PARALLEL.so3 <- function(start, end, vec){
  end %*% t(start) %*% vec
}


#' Project a point onto a geodesic
#' 
#' @param C A list of SO(3)
#' @param geo A 2-list of point and the direction list(R,H). Thus, R'H is skew-symm.
#' @return Projected point of \code{C} to \code{geo}.
#' @export
project.so3 <- function(C, geo){
  n <- length(C)
  A <- geo[[1]]
  H <- geo[[2]]
  H <- H/frob.norm(H)
  geo_curve <- lapply((1:360)/180*sqrt(2)*pi, function(i) EXP.so3(A,H*i))
  mins <- sapply(1:n, function(i) 
    which.min(sapply(1:360,function(j) DIST.so3(C[[i]],geo_curve[[j]]))))
  geo_curve[mins]
}

#' Unroll a piecewise geodesic curve
#' 
#' @param data knots of the curve.
#' @return Unrolled curve of the data on T_{data[1,]}(SO3)
#' @export
UNROLL.so3 <- function(data){
  n <- length(data)
  unroll.mat <- vector("list", length = n)
  unroll.mat[[1]] <- matrix(0,ncol=3,nrow=3) #data[[1]]
  for(i in 2:n){ 
    pp <- LOG.so3(data[[i-1]],data[[i]])
    if(i >= 3){
      for(j in (i-1):2){
        pp <- PARALLEL.so3(data[[j]],data[[j-1]],pp)
      }
    }
    unroll.mat[[i]] <- pp + unroll.mat[[i-1]]
  }
  return(unroll.mat)
}


#' Unwrap points with timestamp
#' 
#' @param data2 points on the SO(3)
#' @param timestamp timestamp of each data, corresponding to each row of \code{data2}.
#' @param base base curve to unwrap w.r.t.
#' @return Unwrapped data on T_{base[1,]}(SO(3)).
#' @export
UNWRAP.so3 <- function(data2, timestamp, base){
  n <- length(data2)
  base.unroll <- UNROLL.so3(base)
  data2.unwrap <- vector("list", length = n)
  
  for(i in 1:n){
    if(timestamp[i] == 1){
      data2.unwrap[[i]] <- LOG.so3(base[[1]],data2[[i]])
    }else{
      t <- timestamp[i]
      data2.log <- LOG.so3(base[[t]],data2[[i]])
      for(j in t:2){
        data2.log <- PARALLEL.so3(base[[j]],base[[j-1]],data2.log)
      }
      data2.unwrap[[i]] <- data2.log + base.unroll[[t]]
    }
  }
  return(data2.unwrap)
}


#' Wrap points with timestamp
#' 
#' @param data2 points on the SO(3)
#' @param timestamp timestamp of each data, corresponding to each row of \code{data2}.
#' @param base base curve to wrap w.r.t.
#' @return Wrapped data on SO(3)
#' @export
WRAP.so3 <- function(data2, timestamp, base){ 
  n <- length(data2)
  base.unroll <- UNROLL.so3(base)
  data2.wrap <- vector("list",length = n)
  
  for(i in 1:n){
    if(timestamp[i] == 1){
      data2.wrap[[i]] <- EXP.so3(base[[1]],data2[[i]])
    }else{
      t <- timestamp[i]
      v <- data2[[i]] - base.unroll[[t]]
      for(j in 2:t){
        v <- PARALLEL.so3(base[[j-1]],base[[j]],v)
      }
      v <- t(base[[t]]) %*% v
      diag(v) <- 0
      v <- (v - t(v))/2
      v <- base[[t]] %*% v
      data2.wrap[[i]] <- EXP.so3(base[[t]],v)
    }
  }
  return(so3(data2.wrap))
}


#' Normalizing constant of SO(3)
#' 
#' @param sigma scale param. of the norm. distn.
#' @return Normalizing constant of the Riemannian normal distribution of the given scale parameter.
#' @export
normal_const.so3 <- function(sigma){
  formula <- function(x) exp(-x^2/2/sigma^2) * sin(x/sqrt(8))^2
  integrate(formula,0,pi*sqrt(2))$value * 32 * pi
}


#' Frechet mean (Intrinsic mean)
#' 
#' @param data each row is data point in SO(3)
#' @param alpha step size
#' @return frechet mean
#' @export
FRECHET.so3 <- function (data, alpha = 1) {
  n <- length(data)
  mu <- data[[1]]
  delta.magnitude <- 1
  
  while(delta.magnitude > 1e-3){
    delta <- matrix(0,ncol=3,nrow=3)
    for(j in 1:n){
      delta <- delta + LOG.so3(mu, data[[j]])
    }
    delta.magnitude <- sqrt(sum(delta^2))
    mu <- EXP.so3(mu, alpha * delta/n)
  }
  return(mu)
}


#' Principal Geodesic Analysis (PGA)
#' 
#' @param data list of SO(3) objects
#' @param fmean to reduce computation time, plug in the precomputed frechet mean
#' @return geodesic parametrized by 2-list
#' @export
PGA.so3 <- function(data, f.mean=NULL){
  n <- length(data)
  if(is.null(f.mean)) f.mean <- FRECHET.so3(data)
  theta <- seq(0,2*pi,length.out = 30)
  phi <- seq(0,pi,length.out=30)
  dir <- matrix(0,ncol=3,nrow=3)
  res.sum <- 0
  for(i in seq_along(theta)){
    for(j in seq_along(phi)){
      xyz <- sph2car(theta[i], phi[j], radius = 1, deg = FALSE)
      H <- matrix(c(0,xyz[1],xyz[2], 0,0,xyz[3], 0,0,0),ncol=3)
      H <- H - t(H)
      H <- f.mean %*% H
      pp <- project.so3(data,list(f.mean,H))
      pv.sum <- sapply(1:n, function(k) DIST.so3(f.mean,pp[[k]])^2)
      pv.sum <- sum(pv.sum)
      if(pv.sum > res.sum){
        res.sum <- pv.sum
        dir <- H
      }
    }
    cat(theta[i], "\t")
  }
  return(list(f.mean,H))
}

#' Sigma: rule-of-thumb
#' 
#' 0.1n-nearest neighborhood, median*qchisq(0.8,3)
#' @param data list of SO(3)
#' @return rule-of-thumb sigma value
#' @export
sigma_thumb.so3 <- function(data){
  n <- length(data)
  dist.mat <- matrix(ncol = n, nrow = n)
  for(i in 1:n){
    for(j in 1:n){
      dist.mat[i,j] <- DIST.so3(data[[i]], data[[j]])
    }
  }
  
  # mean 0.1-nearest neigh distance
  nn.cov <- sapply(1:n, function(r) sqrt(sum(sort(dist.mat[r,])[2:(n/10+1)]^2)/(n/10)))
  sigma.thumb <- median(nn.cov) * sqrt(qchisq(0.8,3))
  sigma.thumb
}


PPC.so3 <- function(data, lambda, max.iter = 20, vis = FALSE, geo.pga = NULL){
  n <- length(data)
  if(is.null(geo.pga)) geo.pga <- PGA.so3(data)
  f <- project.so3(data, geo.pga)
  d <- matrix(0,ncol=n,nrow=n)
  for(i in 1:n){
    for(j in 1:n){
      d[i,j] <- DIST.so3(f[[i]],f[[j]])
    }
  }
  d <- d+1e-10
  diag(d) <- 0
  
  pbit <- progress_bar$new(
    format = paste("ITER :its / Gradient Ascent at :what in :elapsed"),
    clear = FALSE, total = 1e7, width = 60)
  pbit$tick(tokens = list(its = 1,what="waiting..."))
  
  fit <- isoMDS(d, k=1)$points
  f <- f[order(fit)]
  sigma.max <- Inf
  a.const <- 1
  a <- vector(length = n-1)
  for(i in seq_along(a)){
    a[i] <- DIST.so3(f[[i]],f[[i+1]])
  }
  a[a<1e-5] <- 1e-5 ## if not set, singular when fitting
  a <- c(0,cumsum(a))/sum(a) *a.const
  sigma.thumb <- sigma_thumb.so3(data)
  sigma <- rep(sigma.thumb,n)
  v <- rep(1/n,n)
  w <- matrix(ncol = n, nrow = n)
  for(i in 1:n){
    for(k in 1:n){
      d <- DIST.so3(data[[i]],f[[k]])
      g <- exp(-d^2/2/sigma[k]^2)/normal_const.so3(sigma[k])
      w[i,k] <- g * v[k]
    }
  }
  w <- w/rowSums(w)
  timestamp <- apply(w, 1, which.max)
  
  ## iteration 1
  data.unwrap <- UNWRAP.so3(data, timestamp, f)
  data.unwrap.skew <- lapply(1:n, function(i) t(f[[1]]) %*% data.unwrap[[i]])
  y1 <- sapply(1:n, function(i) data.unwrap.skew[[i]][1,2])
  y2 <- sapply(1:n, function(i) data.unwrap.skew[[i]][1,3])
  y3 <- sapply(1:n, function(i) data.unwrap.skew[[i]][2,3])
  iter1 <- ppc.iter.e.3d(cbind(y1,y2,y3), w, a, v, sigma, sigma.max = sigma.max, lambda)
  f.mat <- iter1
  f.tan <- vector("list", length = n)
  for(i in 1:n){
    f.tan[[i]] <- matrix(c(0,f.mat[i,1],f.mat[i,2],
                           -f.mat[i,1],0,f.mat[i,3],
                           -f.mat[i,2],-f.mat[i,3],0),byrow=T, ncol=3)
    # if(frob.norm(f.tan[[i]]) >= sqrt(2)*pi) stop()
    f.tan[[i]] <- f[[1]] %*% f.tan[[i]]
  }
  f.new <- WRAP.so3(f.tan, 1:n, f)
  
  if(vis) plot(data, red = f.new, gray = f)
  
  ## SGD
  tau <- 1/sigma^2
  sgd.const2 <-  sapply(1:n, function(k) 
    sum(sapply(1:n, function(i) -1/2*w[i,k]*DIST.so3(data[[i]],f.new[[k]])^2)))
  sgd.const1 <- vector(length = n)
  alpha <- rep(1, n) #step-size
  grad <- rep(1,n)
  iters <- 0
  while(sum(grad[alpha!=0]^2)>1e-7){
    pbit$tick(tokens = list(its=1,what = sum(grad[alpha!=0]^2)))
    iters <- iters + 1
    for(k in 1:n){
      formula <- function(r) r^2/2*exp(-r^2/2/sigma[k]^2)*sin(r/sqrt(8))^2
      sgd.const1[k] <- 32*pi*integrate(formula, 0, pi)$value/normal_const.so3(sigma[k])
    }
    sgd.const1 <- colSums(w) * sgd.const1
    grad.old <- grad
    grad <- sgd.const1 + sgd.const2
    if(sum((grad.old-grad)^2)<1e-7) alpha <- alpha * 1.1
    alpha[which(abs(grad) < 1e-7)] <- 0
    # alpha <- alpha*0.99
    tau <- tau + alpha*grad
    tau[tau < 0] <- min(abs(tau[tau!=0]))
    sigma <- sqrt(1/tau)
  }
  
  ## iteration >= 2
  kk <- 2
  while(TRUE){
    message(paste("iter ",kk))
    if(max(sapply(1:n, function(i) DIST.so3(f[[i]],f.new[[i]]))) < 1e-2 || kk > max.iter) break()
    pbit$tick(tokens = list(its = kk,what="waiting..."))
    f <- f.new
    a <- vector(length = n-1)
    for(i in seq_along(a)){
      a[i] <- DIST.so3(f[[i]],f[[i+1]])
    }
    a[a==0] <- 1e-5
    a <- c(0,cumsum(a))/sum(a) *a.const
    v <- colSums(w)/n
    w <- matrix(ncol = n, nrow = n)
    for(i in 1:n){
      for(k in 1:n){
        d <- DIST.so3(data[[i]],f[[k]])
        g <- exp(-d^2/2/sigma[k]^2)/normal_const.so3(sigma[k])
        w[i,k] <- g * v[k]
      }
    }
    w <- w/rowSums(w)
    if(!all(!is.nan(w)) || !all(colSums(w) != 0)){
      message("stopped because some latent variables are useless")
      break()
    }
    timestamp <- apply(w, 1, which.max)
    data.unwrap <- UNWRAP.so3(data, timestamp, f)
    data.unwrap.skew <- lapply(1:n, function(i) t(f[[1]]) %*% data.unwrap[[i]])
    y1 <- sapply(1:n, function(i) data.unwrap.skew[[i]][1,2])
    y2 <- sapply(1:n, function(i) data.unwrap.skew[[i]][1,3])
    y3 <- sapply(1:n, function(i) data.unwrap.skew[[i]][2,3])
    iter2 <- ppc.iter.e.3d(cbind(y1,y2,y3), w, a, v, sigma, sigma.max = sigma.max, lambda)
    f.mat <- iter2
    f.tan <- vector("list", length = n)
    for(i in 1:n){
      f.tan[[i]] <- matrix(c(0,f.mat[i,1],f.mat[i,2],
                             -f.mat[i,1],0,f.mat[i,3],
                             -f.mat[i,2],-f.mat[i,3],0),byrow=T, ncol=3)
      # if(frob.norm(f.tan[[i]]) >= sqrt(2)*pi) stop()
      f.tan[[i]] <- f[[1]] %*% f.tan[[i]]
    }
    f.new <- WRAP.so3(f.tan, 1:n, f)
    
    if(vis){
      plot(data,red=f.new,gray=f, main = paste("ITER",kk,"/ red:new, gray=old"))
    }
    
    tau <- 1/sigma^2
    tau.old <- tau
    sigma.old <- sigma
    sgd.const2 <-  sapply(1:n, function(k) 
      sum(sapply(1:n, function(i) -1/2*w[i,k]*DIST.so3(data[[i]],f.new[[k]])^2)))
    sgd.const1 <- vector(length = n)
    alpha <- rep(1, n) #step-size
    grad <- rep(1,n)
    iters <- 0
    while(sum(grad[alpha!=0]^2)>1e-7){
      pbit$tick(tokens = list(its=kk,what = sum(grad[alpha!=0]^2)))
      iters <- iters + 1
      for(k in 1:n){
        formula <- function(r) r^2/2*exp(-r^2/2/sigma[k]^2)*sin(r/sqrt(8))^2
        sgd.const1[k] <- 32*pi*integrate(formula, 0, pi)$value/normal_const.so3(sigma[k])
      }
      sgd.const1 <- colSums(w) * sgd.const1
      grad.old <- grad
      grad <- sgd.const1 + sgd.const2
      if(sum((grad.old-grad)^2)<1e-7) alpha <- alpha * 1.1
      alpha[which(abs(grad) < 1e-7)] <- 0
      if(iters > 1000) return(list(f,sigma.old))
      tau <- tau + alpha*grad
      tau[tau < 0] <- min(abs(tau[tau!=0]))
      sigma <- sqrt(1/tau)
    }
    
    kk <- kk + 1
  }
  return(list(f.new, sigma))
}
