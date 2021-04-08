source("../generics.R")
source("../ppc_euc.R")

library(progress)

#' define class `half_plane`
new_half_plane <- function(x = double(), y = double()){
  stopifnot(y > 0 && is.double(x) && is.double(y))
  structure(cbind(x, y), class = "half_plane")
}

#' Constructor of class \code{half_plane} 
#' 
#' @param x A vector, or a matrix with 2 columns.
#' @param y A vector. May be missing when x is a matrix.
#' @return n x 2 matrix with coordinates \code{x} and \code{y}.
#' @export
#' @examples
#' x <- rnorm(100)
#' y <- runif(100, 1, 2)
#' half_plane(x, y)
half_plane <- function(x, y=NULL){
  if(is.null(y)){
    if(is.vector(x)) return(new_half_plane(x[1],x[2]))
    return(new_half_plane(as.double(x[,1]), as.double(x[,2])))
  }
  y <- as.double(y)
  new_half_plane(x, y)
}

`[.half_plane` <- function(x, i, j, ...) {
  out <- unclass(x)
  if(is.vector(out)){
    out <- `[`(out,i)
    if(length(out) == 2){
      return(half_plane(out))
    }
    return(out)
  }
  if(dim(out)[1] == 1L){
    out <- out[i]
    if(length(out) == 2){
      return(half_plane(out))
    }
    return(out)
  }
  if(missing(j)){
    if(missing(i)) i <- 1:nrow(out)
    out <- out[i,,...]
    return(half_plane(out))
  }
  if(length(j) != 2){
    out <- out[i,j,...]
    return(out)
  }
  return(half_plane(out[i,j,...]))
}

## -----------------------------------------------------------------------------

#' Draws a geodesic between two points.
#' 
#' There should be an existing plot to draw the geodesic onto.
#' Vectorized operations are not supported. Use \code{geodesics.half_plane}.
#' @param a A 2-vector.
#' @param b A 2-vector.
#' @return No return. Draws a geodesic. 
#' @export
#' @examples
#' plot(x=c(1,2),y=c(1,2), xlim = c(1,3), ylim=c(0,2))
#' geodesic.half_plane(c(1,1),c(2,2))
geodesic.half_plane <- function(a,b,...){
  x1 <- a[1]; x2 <- b[1]
  y1 <- a[2]; y2 <- b[2]
  if(abs(x1 - x2) < 1e-10){
    lines(rbind(a,b),...)
  }else{
    z <- (y1^2-y2^2+x1^2-x2^2)/2/(x1-x2) # x-intersection
    r <- sqrt((x1-z)^2 + y1^2)
    theta1 <- acos((x1-z)/r)
    theta2 <- acos((x2-z)/r)
    theta <- seq(min(theta1,theta2),max(theta1,theta2),by=abs(theta2-theta1)/100)
    lines(z + r*cos(theta), r*sin(theta),...)
  }
}

#' Draw geodesics between data
#' 
#' Draw geodesics connecting each row of the given data onto the existing plot.
#' @param data An n x 2 matrix.
#' @return No return. Draws a geodesic. 
#' @export
#' @examples
#' plot(x=c(1,2,3),y=c(1,2,.5), xlim = c(1,3), ylim=c(0,2))
#' geodesics.half_plane(rbind(c(1,1),c(2,2),c(3,.5)), col = "red")
geodesics.half_plane <- function(data, ...){
  n <- nrow(data)
  for(i in 1:(n-1)){
    geodesic.half_plane(data[i,],data[i+1,], ...)
  }
}

#' Draws a full geodesic
#' 
#' @param geo A 2-vector consisting of parameters of a geodesic. Center of circle and the radius. If the geodesic is a vertical line, the x-coord of the line and Inf. 
#' @return No return. Just draws a geodesic onto the existing plot.
#' @export
#' @examples
#' geo1 <- c(z=-2,r=2)
#' geo2 <- c(z=-1,r=Inf)
#' plot(data)
#' geodesic_full.half_plane(geo1, col = "red")
#' geodesic_full.half_plane(geo2, col = "blue")
geodesic_full.half_plane <- function(geo, ...){
  z <- geo[1]; r <- geo[2]
  if(r == Inf){
    abline(v=z, ...)
  }else{
    theta <- seq(0,pi,by=pi/100)
    lines(z + r * cos(theta), r*sin(theta), ...)
  }
}


## -----------------------------------------------------------------------------

#' Geodesic distance between two points. 
#' 
#' Vectorized operations are not supported.
#' @param a A 2-vector.
#' @param b A 2-vector.
#' @return The geodesic distance between \code{a} and \code{b}. 
#' @export
#' @examples
#' DIST.half_plane(c(1,3), c(2,2))
DIST.half_plane <- function(a,b){
  x1 <- a[1]; x2 <- b[1]
  y1 <- a[2]; y2 <- b[2]
  2*log((sqrt((x2-x1)^2 + (y1-y2)^2) + sqrt((x2-x1)^2 + (y1+y2)^2))/2/sqrt(y1*y2))
}


#' Returns a logarithmic map
#' 
#' @param p A 2-vector on H^2.
#' @param v A 2-vector on H^2. 
#' @return Log_p(v). 
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100, -4,4)
#' x <- sort(x)[seq(1,100,by=10)]
#' y <- runif(100, 0, 4)[seq(1,100,by=10)]
#' data <- half_plane(x, y)
#' 
#' plot(data, main="LOG map"); 
#' geodesics.half_plane(data, col = "red")
#' for(i in 1:9){
#'   pp <- LOG(data[i+1,], data[i,])
#'   lines(x=c(data[i+1,1],data[i+1,1]+pp[1]),y=c(data[i+1,2],data[i+1,2]+pp[2]))
#' }
#' for(i in 1:9){
#'   pp <- LOG(data[i,],data[i+1,])
#'   lines(x=c(data[i,1],data[i,1]+pp[1]),y=c(data[i,2],data[i,2]+pp[2]))
#' }
LOG.half_plane <- function(p,v){
  x1 <- p[1]; x2 <- v[1]
  y1 <- p[2]; y2 <- v[2]
  d <- DIST.half_plane(p,v)
  if(abs(x1-x2)<1e-5) return(c(0,sign(y2-y1)*d))
  z <- (y1^2-y2^2+x1^2-x2^2)/2/(x1-x2)
  r <- sqrt((x1-z)^2 + y1^2)
  theta1 <- acos((x1-z)/r)
  if(x1 < x2) return(c(d*cos(theta1-pi/2), d*sin(theta1-pi/2)))
  return(c(-d*cos(theta1-pi/2), -d*sin(theta1-pi/2)))
}

#' Returns an exponential map
#' 
#' @param p A 2-vector on H^2.
#' @param v A 2-vector on T_p(H^2). 
#' @return Exp_p(v). 
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100, -4,4)
#' x <- sort(x)[seq(1,100,by=10)]
#' y <- runif(100, 0, 4)[seq(1,100,by=10)]
#' data <- half_plane(x, y)
#' 
#' plot(data, main="EXP map",ylim=c(0,5),xlim = c(-4,5)); 
#' for(i in 1:9){
#'   pp <- EXP(data[i+1,],data[i,]-data[i+1,])
#'   points(pp[1],pp[2], pch = "+")
#'   geodesic.half_plane(data[i+1,],pp)
#' }
EXP.half_plane <- function(p,v){
  v <- v+p
  x1 <- p[1]; x2 <- v[1]
  y1 <- p[2]; y2 <- v[2]
  d <- sqrt(sum((p-v)^2))
  if(abs(x1-x2)<1e-5) {
    return(c(x1, y1*exp(sign(y2-y1)*d)))
  }
  z <- y1*(y2-y1)/(x2-x1)+x1
  r <- sqrt((x1-z)^2 + y1^2)
  theta1 <- acos((x1-z)/r)
  if(x2 < x1){
    theta <- atan(exp(d)*tan(theta1/2))*2
    out <- c(z+r*cos(theta), r*sin(theta))
    class(out) <- "half_plane"
    return(out)
  }else{
    theta <- atan(exp(-d)*tan(theta1/2))*2
    out <- c(z+r*cos(theta), r*sin(theta))
    class(out) <- "half_plane"
    return(out)
  }
}

#' Parallel transport a vector on the tangent space via geodesic
#' 
#' @param start A 2-vector on H^2.
#' @param end A 2-vector on H^2. 
#' @param vec A 2-vector on T_start(H^2)
#' @return Parallel transport \code{vec} from \code{start} to \code{end} via a geodesic. 
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100, -4,4)
#' x <- sort(x)[seq(1,100,by=10)]
#' y <- runif(100, 0, 4)[seq(1,100,by=10)]
#' data <- half_plane(x, y)
#' 
#' move.vec <- c(1,2)
#' plot(data, main="parallel transport(red)", ylim = c(-1,5),xlim=c(-4,4))
#' geodesics.half_plane(data, col = "blue")
#' lines(rbind(data[1,],data[1,]+move.vec), col = "red")
#' for(i in 1:9){
#'   move.vec <- PARALLEL(data[i,],data[i+1,], move.vec)
#'   points(rbind(data[i+1,]+move.vec), pch = "+")
#'   lines(rbind(data[i+1,],data[i+1,]+move.vec))
#' }
PARALLEL.half_plane <- function(start, end, vec){
  if(sum((start-end)^2) < 1e-5) return(vec)
  x1 <- start[1]; x2 <- end[1]
  y1 <- start[2]; y2 <- end[2]
  z <- (y1^2-y2^2+x1^2-x2^2)/2/(x1-x2) # x-intersection
  r <- sqrt((x1-z)^2 + y1^2)
  theta1 <- acos((x1-z)/r)
  theta2 <- acos((x2-z)/r)
  theta <- theta2 - theta1
  v <- vec #- start
  return(c(v[1]*cos(theta) - v[2]*sin(theta), v[1]*sin(theta) + v[2]*cos(theta)))
}


#' Project a point onto a geodesic
#' 
#' @param a A 2-vector on H^2.
#' @param geo A 2-vector parametrization of a geodesic
#' @param draw If TRUE, draws the projection path onto the existing plot.
#' @return Projected point of \code{a} to \code{geo}.
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100, -4,4)
#' x <- sort(x)[seq(1,100,by=10)]
#' y <- runif(100, 0, 4)[seq(1,100,by=10)]
#' data <- half_plane(x, y)
#' 
#' geo1 <- c(z=-2,r=2)
#' plot(data); geodesic_full.half_plane(geo1, col = "red"); 
#' for(i in 1:10){
#'   points(x=project.half_plane(data[i,],geo1)[1], 
#'          y=project.half_plane(data[i,],geo1,TRUE)[2], 
#'          pch = as.character(i))
#' }
project.half_plane <- function(a, geo, draw = FALSE){
  z <- geo[1]; r <- geo[2]
  x <- a[1]; y <- a[2]
  if(r == Inf){
    rp <- sqrt(sum((a-c(z,0))^2))
    if(draw == TRUE) geodesic.half_plane(c(z, rp), a)
    return(c(x=z,y=rp))
  }
  if(abs(x-z)<1e-5) return(c(x=z, y=r))
  zp <- (x^2 + y^2 + r^2 - z^2)/2/(x-z)
  rp <- sqrt(sum((a-c(zp,0))^2))
  if(z < zp){
    theta <- atan(rp/r)
  }else{
    theta <- pi - atan(rp/r)
  }
  ret <- c(x=z + r * cos(theta), y= r * sin(theta))
  if(draw == TRUE) geodesic.half_plane(ret, a)
  return(ret)
}

#' Unroll a piecewise geodesic curve
#' 
#' @param data knots of the curve.
#' @return Unrolled curve of the data on T_{data[1,]}(H^2)
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100, -4,4)
#' x <- sort(x)[seq(1,100,by=10)]
#' y <- runif(100, 0, 4)[seq(1,100,by=10)]
#' data <- half_plane(x, y)
#' 
#' dt <- UNROLL(data) + matrix(rep(data[1,],10),byrow=TRUE,ncol=2)
#' plot(data, ylim = c(0,5), main = "black path = unrolled"); 
#' geodesics.half_plane(data,col = "red")
#' for(i in 1:9){
#'   pp <- LOG(data[i,],data[i+1,])
#'   lines(rbind(data[i,],data[i,]+pp))
#' }
#' points(dt, pch = "+"); lines(dt)
UNROLL.half_plane <- function(data){
  n <- nrow(data)
  unroll.mat <- matrix(ncol = 2, nrow = n)
  unroll.mat[1,] <- c(0,0)
  for(i in 2:n){
    pp <- LOG.half_plane(data[i-1,],data[i,])
    if(i >= 3){
      for(j in (i-1):2){
        pp <- PARALLEL.half_plane(data[j,],data[j-1,],pp)
      }
    }
    unroll.mat[i,] <- pp + unroll.mat[i-1,]
  }
  return(unroll.mat)
}

#' Unwrap points with timestamp
#' 
#' @param data2 points on the half plane
#' @param timestamp timestamp of each data, corresponding to each row of \code{data2}.
#' @param base base curve to unwrap w.r.t.
#' @return Unwrapped data on T_{base[1,]}(H^2)
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100, -4,4)
#' x <- sort(x)[seq(1,100,by=10)]
#' y <- runif(100, 0, 4)[seq(1,100,by=10)]
#' data <- half_plane(x, y)
#' set.seed(0)
#' data2 <- data + rnorm(20)
#' timestamp <- seq(1:10)
#' timestamp[4] <- 3
#' 
#' plot(data, ylim = c(0,5), xlim = c(-5,6)); 
#' points(data2, pch = as.character(timestamp), col = "red")
#' geodesics.half_plane(data,col = "red")
#' data2.unwrap <- UNWRAP(data2, timestamp, data)
#' dt <- UNROLL(data)
#' plot(data2.unwrap, pch = as.character(timestamp), main="Tangent space");
#' points(dt, pch = "+"); lines(dt)
UNWRAP.half_plane <- function(data2, timestamp, base){
  n <- nrow(data2)
  base.unroll <- UNROLL.half_plane(base)
  data2.unwrap <- matrix(ncol=2,nrow=n)
  
  for(i in 1:n){
    if(timestamp[i] == 1){
      data2.unwrap[i,] <- LOG.half_plane(base[1,],data2[i,])
    }else{
      t <- timestamp[i]
      data2.log <- LOG.half_plane(base[t,],data2[i,])
      for(j in t:2){
        data2.log <- PARALLEL.half_plane(base[j,],base[j-1,],data2.log)
      }
      data2.unwrap[i,] <- data2.log + base.unroll[t,]
    }
  }
  return(data2.unwrap)
}



#' Wrap points with timestamp
#' 
#' @param data2 points on the tangent space of half plane
#' @param timestamp timestamp of each data, corresponding to each row of \code{data2}.
#' @param base base curve to wrap w.r.t.
#' @return Wrapped data on H^2
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100, -4,4)
#' x <- sort(x)[seq(1,100,by=10)]
#' y <- runif(100, 0, 4)[seq(1,100,by=10)]
#' data <- half_plane(x, y)
#' set.seed(0)
#' data2 <- data + rnorm(20)
#' timestamp <- seq(1:10)
#' timestamp[4] <- 3
#' 
#' plot(data, xlim = c(-4,5), ylim = c(-1,5), main = "wrap points")
#' points(data2, pch = "+")
#' data2.unwrap <- UNWRAP(data2, timestamp, data)
#' dt.wrap <- WRAP.half_plane(data2.unwrap, timestamp, data)
#' points(dt.wrap, pch = as.character(timestamp), col = "red")
WRAP.half_plane <- function(data2, timestamp, base){
  n <- nrow(data2)
  base.unroll <- UNROLL.half_plane(base)
  data2.wrap <- matrix(ncol=2,nrow=n)
  
  for(i in 1:n){
    if(timestamp[i] == 1){
      data2.wrap[i,] <- EXP.half_plane(base[1,],data2[i,])
    }else{
      t <- timestamp[i]
      v <- data2[i,] - base.unroll[t,]
      for(j in 2:t){
        v <- PARALLEL.half_plane(base[j-1,],base[j,],v)
      }
      data2.wrap[i,] <- EXP.half_plane(base[t,],v)
    }
  }
  return(data2.wrap)
}

#' Normal distribution plot
#' 
#' conf.level does not necessarily mean confidence level. It is used for convenience for deciding the size of circle drawn.
#' @param mu location param. of the norm. distn.
#' @param sigma scale param. of the norm. distn.
#' @param conf.level how big should the circle be?
#' @return No return. Draw a circle on the existing plot.
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100, -4,4)
#' x <- sort(x)[seq(1,100,by=10)]
#' y <- runif(100, 0, 4)[seq(1,100,by=10)]
#' data <- half_plane(x, y)
#' 
#' plot(data, ylim = c(-0.1,4)); 
#' for(i in 1:10) normal_plot.half_plane(mu = data[i,], sigma = 1/4, col="red")
normal_plot.half_plane <- function(mu, sigma, conf.level = 0.95,...){
  k <- sqrt(qchisq(conf.level, df = 2))
  theta <- seq(0,2*pi, by = pi/360)
  e <- matrix(nrow = length(theta),ncol=2)
  for(t in seq_along(theta)){
    e[t,] <- EXP.half_plane(mu, k*sigma*c(cos(theta[t]),sin(theta[t])))
  }
  # points(mu[1], mu[2], pch = "+", ...)
  lines(e,...)
}

#' Normalizing constant of half-plane
#' 
#' @param sigma scale param. of the norm. distn.
#' @return Normalizing constant of the Riemannian normal distribution of the given scale parameter.
#' @export
normal_const.half_plane <- function(sigma){
  formula <- function(x) exp(x-x^2/2/sigma^2)/2 - exp(-x-x^2/2/sigma^2)/2
  integrate(formula,0,Inf)$value * 2 * pi
}

#' Frechet mean (Intrinsic mean)
#' 
#' @param data each row is data point in half-plane
#' @param alpha step size
#' @return frechet mean
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100, -4,4)
#' x <- sort(x)[seq(1,100,by=10)]
#' y <- runif(100, 0, 4)[seq(1,100,by=10)]
#' data <- half_plane(x, y)
#' 
#' f.mu <- FRECHET(data)
#' plot(data)
#' points(f.mu[1], f.mu[2], pch = "f")
#' points(colMeans(data)[1], colMeans(data)[2], pch = "m")
#' title("f.mean = 'f', e.mean  = 'm'")
FRECHET.half_plane <- function(data, alpha = 1, weight = NULL, vis = F){
  n <- nrow(data)
  if(is.null(weight)) weight = rep(1,n)
  mu <- colMeans(data)
  delta.magnitude <- 1
  
  while(delta.magnitude > 1e-3){
    delta <- c(0,0)
    for(j in 1:n){
      if(weight[j] == 0) next()
      delta <- delta + LOG.half_plane(mu, data[j,])*weight[j]
    }
    if(abs(sqrt(sum(delta^2)) - delta.magnitude)<1e-4 && delta.magnitude < 1e-2) alpha <- alpha*2
    delta.magnitude <- sqrt(sum(delta^2))
    mu <- EXP.half_plane(mu, alpha * delta/n)
    if(vis) cat(paste0(delta.magnitude,"\t"))
  }
  return(mu)
}

#' Principal Geodesic Analysis (PGA)
#' 
#' @param data each row is data point in half-plane
#' @return geodesic parametrized by 2-vector
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100, -4,4)
#' x <- sort(x)[seq(1,100,by=10)]
#' y <- runif(100, 0, 4)[seq(1,100,by=10)]
#' data <- half_plane(x, y)
#' 
#' geo.pga <- PGA(data)
#' fm <- FRECHET(data)
#' plot(data)
#' geodesic_full.half_plane(geo.pga, col = "red")
#' points(fm[1], fm[2], pch = "+")
PGA.half_plane <- function(data){
  n <- nrow(data)
  fm <- FRECHET.half_plane(data)
  dir <- seq(-pi,0,length.out = 102)[2:101]
  z <- tan(dir)*fm[2]+fm[1]
  r <- 0
  z.max.val <- 0
  z.max.index <- 1
  
  for(j in 1:100){
    if(abs(z[j]-fm[1]) < 1e-5){
      r = Inf
    }else{
      r = sqrt(sum((fm-c(z[j],0))^2))
    }
    temp <- 0
    for(i in 1:n){ 
      temp <- temp + DIST.half_plane(project.half_plane(data[i,],c(z[j],r)),fm)^2
    }
    if(temp > z.max.val){
      z.max.val <- temp
      z.max.index <- j
    }
  }
  z <- z[z.max.index]
  if(abs(z-fm[1]) < 1e-5){
    r = Inf
  }else{
    r = sqrt(sum((fm-c(z,0))^2))
  } 
  return(c(z,r))
}


#' Sigma: rule-of-thumb
#' 
#' 0.1n-nearest neighborhood, median*qchisq(0.8,2)
#' @param data each row is data point in half-plane
#' @return rule-of-thumb sigma value
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100, -4,4)
#' x <- sort(x)[seq(1,100,by=10)]
#' y <- runif(100, 0, 4)[seq(1,100,by=10)]
#' data <- half_plane(x, y)
#' 
#' sigma_thumb.half_plane(data)
sigma_thumb.half_plane <- function(data){
  n <- nrow(data)
  dist.mat <- matrix(ncol = n, nrow = n)
  for(i in 1:n){
    for(j in 1:n){
      dist.mat[i,j] <- DIST.half_plane(data[i,], data[j,])
    }
  }
  
  # mean 0.1-nearest neigh distance
  nn.cov <- sapply(1:n, function(r) sqrt(sum(sort(dist.mat[r,])[2:(n/10+1)]^2)/(n/10)))
  sigma.thumb <- median(nn.cov) * sqrt(qchisq(0.8,2))
  sigma.thumb
}

#' Fit a probabilistic principal curve
#' 
#' @param data A nx2 matrix of each row on H^2.
#' @param lambda A nonnegative smoothing parameter. 
#' @max.iter maximum number of iterations to be done.
#' @vis If \code{TRUE}. show the fitted PPC of each iteration in 1*2 plot.
#' @geo.pga  Manually set the initial curve other than the geodesic.
#' @return Fitted PPC consisting of list(f, sigma)
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100, -4,4)
#' x <- sort(x)[seq(1,100,by=10)]
#' y <- runif(100, 0, 4)[seq(1,100,by=10)]
#' data <- half_plane(x, y)
#' 
#' PPC(data, lambda = 0.1)
PPC.half_plane <- function(data, lambda, max.iter = 20, vis = FALSE, geo.pga = NULL){
  n <- nrow(data)
  if(is.null(geo.pga)) geo.pga <- PGA.half_plane(data)
  
  pbit <- progress_bar$new(
    format = paste("ITER :its / Gradient Ascent at :what in :elapsed"),
    clear = FALSE, total = 1e7, width = 60)
  pbit$tick(tokens = list(its = 1,what="waiting..."))
  
  f <- t(apply(data, 1, function(row) project.half_plane(row, geo.pga)))
  f <- f[order(f[,1]),]
  f.init <- f
  sigma.max <- Inf
  a.const <- 1
  a <- vector(length = n-1)
  for(i in seq_along(a)){
    a[i] <- DIST.half_plane(f[i,],f[i+1,])
  }
  a <- c(0,cumsum(a))/sum(a) *a.const
  sigma.thumb <- sigma_thumb.half_plane(data)
  sigma <- rep(sigma.thumb,n)
  v <- rep(1/n,n)
  w <- matrix(1, nrow = n, ncol = n)
  for(i in 1:n){
    for(k in 1:n){
      d <- DIST.half_plane(data[i,],f[k,])
      g <- exp(-d^2/2/sigma[k]^2)/normal_const.half_plane(sigma[k])
      w[i,k] <- g * v[k]
    }
  }
  w <- w/rowSums(w)
  timestamp <- apply(w, 1, which.max)
  
  ## iteration 1
  data.unwrap <- UNWRAP.half_plane(data, timestamp, f)
  iter1 <- ppc.iter.e(data.unwrap, w, a, v, sigma, sigma.max = sigma.max, lambda)
  
  f.new <- WRAP.half_plane(iter1, 1:n, f)
  ## SGD
  tau <- 1/sigma^2
  sgd.const2 <-  sapply(1:n, function(k) 
    sum(sapply(1:n, function(i) -1/2*w[i,k]*DIST.half_plane(data[i,],f.new[k,])^2)))
  sgd.const1 <- vector(length = n)
  alpha <- rep(1, n) #step-size
  grad <- rep(1,n)
  iters <- 0
  
  while(sum(grad[alpha!=0]^2)>1e-7){
    pbit$tick(tokens = list(its=1,what = sum(grad[alpha!=0]^2)))
    iters <- iters + 1
    # cat(sum(grad[alpha!=0]^2), "\t")
    for(k in 1:n){
      formula <- function(r) r^2/2*(exp(r-r^2/2/sigma[k]^2)/2 - exp(-r-r^2/2/sigma[k]^2)/2)
      sgd.const1[k] <- 2*pi*integrate(formula, 0, Inf)$value/normal_const.half_plane(sigma[k])
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
  par(mfrow = c(1,2), oma=c(0,0,3,0))
  while(TRUE){
    if(max(sapply(1:n, function(i) DIST.half_plane(f[i,],f.new[i,]))) < 1e-2 || kk > max.iter) break()
    pbit$tick(tokens = list(its = kk,what="waiting..."))
    
    f <- f.new
    a <- vector(length = n-1)
    for(i in seq_along(a)){
      a[i] <- DIST.half_plane(f[i,],f[i+1,])
    }
    a <- c(0,cumsum(a))/sum(a) *a.const
    v <- colSums(w)/n
    for(i in 1:n){
      for(k in 1:n){
        d <- DIST.half_plane(data[i,],f[k,])
        g <- exp(-d^2/2/sigma[k]^2)/normal_const.half_plane(sigma[k])
        w[i,k] <- g * v[k]
      }
    }
    w <- w/rowSums(w)
    if(!all(!is.nan(w)) || !all(colSums(w) != 0)){
      message("stopped because some latent variables are useless")
      break()
    }
    timestamp <- apply(w, 1, which.max)
    
    data.unwrap <- UNWRAP.half_plane(data, timestamp, f)
    iter2 <- ppc.iter.e(data.unwrap, w, a, v, sigma, sigma.max = sigma.max, lambda)
    
    f.new <- WRAP.half_plane(iter2, 1:n, f)
    ## SGD
    tau <- 1/sigma^2
    sgd.const2 <-  sapply(1:n, function(k) 
      sum(sapply(1:n, function(i) -1/2*w[i,k]*DIST.half_plane(data[i,],f.new[k,])^2)))
    sgd.const1 <- vector(length = n)
    alpha <- rep(0.5, n) #step-size
    grad <- rep(1,n)
    iters <- 0
    
    while(sum(grad[alpha!=0]^2)>1e-7){
      pbit$tick(tokens = list(its=kk,what = sum(grad[alpha!=0]^2)))
      iters <- iters + 1
      for(k in 1:n){
        formula <- function(r) r^2/2*(exp(r-r^2/2/sigma[k]^2)/2 - exp(-r-r^2/2/sigma[k]^2)/2)
        sgd.const1[k] <- 2*pi*integrate(formula, 0, Inf)$value/normal_const.half_plane(sigma[k])
      }
      sgd.const1 <- colSums(w) * sgd.const1
      grad.old <- grad
      grad <- sgd.const1 + sgd.const2
      if(sum((grad.old-grad)^2)<1e-7)  alpha <- alpha * 1.1
      if(iters >= 500) alpha <- alpha * 0.9
      alpha[which(abs(grad) < 1e-7)] <- 0
      # alpha <- alpha*0.99
      tau <- tau + alpha*grad
      tau[tau < 0] <- min(abs(tau[tau!=0]))
      sigma <- sqrt(1/tau)
    }
    
    if(vis){
      plot(data)
      title(paste("H^2 : ",kk," / diff=", max(sapply(1:n, function(i) DIST.half_plane(f[i,],f.new[i,]))) ) )
      geodesics.half_plane(f.new, col = "blue"); geodesics.half_plane(f, col = "red");
      legend("topright", legend = c("old","new"), lty=1, lwd=1, col=c("red","blue"))
      for(i in 1:n){
        normal_plot.half_plane(f.new[i,],sigma[i],col = "blue",lwd=0.2)
      }
    }
    kk <- kk + 1
  }
  
  return(list(f.new, sigma))
}

#' biweight kernel
quartic_kernel <- function(d,q){
  if(d/q > 1) return(0)
  return((1-(d/q)^2)^2)
}

#' Fit a principal curve of Hauberg (2016)
#' 
#' Note that plot is drawn at each iteration step.
#' @param data A nx2 matrix of each row on H^2.
#' @param q A bandwidth of the quartic (biweight) kernel. 
#' @max.iter maximum number of iterations to be done.
#' @return Fitted PC consisting of matrix of knot points
#' @export
PC_Hauberg.half_plane <- function(data, q, max.iter = 20){
  n <- nrow(data)
  f.new <- matrix(0,nrow=n,ncol=2)
  geo.pga <- PGA.half_plane(data)
  f <- t(apply(data, 1, function(row) project.half_plane(row, geo.pga)))
  f <- f[order(f[,1]),]
  iter <- 0
  
  pbit <- progress_bar$new(
    format = paste("ITER :its / point :what in :elapsed"),
    clear = FALSE, total = 1e7, width = 60)
  
  while(TRUE){
    iter <- iter + 1
    pbit$tick(tokens = list(its = iter,what="waiting..."))
    if(iter > max.iter) break()
    #projection index
    t <- sapply(1:n, function(j) which.min(sapply(1:n, function(i) DIST.half_plane(f[i,],data[j,]))))
    f_proj <- f[t,]
    for(t in 1:n){
      pbit$tick(tokens = list(its = iter,what=t))
      wt <- sapply(1:n, function(i) quartic_kernel(DIST.half_plane(f[t,],f_proj[i,]),q))
      f.new[t,] <- FRECHET.half_plane(data,weight=wt, vis = F)
    }
    if(max(rowSums(f-f.new)^2)<1e-3) break()
    f <- f.new
    plot(data); lines(f)
  }
  # plot(data); lines(f,pch=20)
  return(f.new)
}

