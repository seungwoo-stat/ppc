source("../ppc_euc.R")
source("../generics.R")

library(rgl)
library(sphereplot)
library(MASS)
library(progress)

#' define class `sphere2`
new_sphere2 <- function(x = double(), y = double(), z = double()){
  stopifnot(is.double(c(x,y,z)) && all(abs(x^2 + y^2 + z^2 - 1) < 1e-5))
  structure(cbind(x, y, z), class = "sphere2")
}

#' define class `sphere3`
new_sphere3 <- function(x = double(), y = double(), z = double(), w = double()){
  stopifnot(is.double(c(x,y,z,w)) && all(abs(x^2 + y^2 + z^2 + w^2 - 1) < 1e-5))
  structure(cbind(x, y, z, w), class = "sphere3")
}

#' Constructor of class \code{sphere2} 
#' 
#' @param x A vector, or a matrix with 3 columns.
#' @param y A vector. May be missing when x is a matrix.
#' @param z A vector. May be missing when x is a matrix.
#' @return n x 3 matrix with coordinates \code{x}, \code{y}, and \code{z}.
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100)
#' y <- sqrt((1-x^2)/3)
#' z <- sqrt((1-x^2)/3*2)
#' s1 <- sphere2(x,y,z)
#' s2 <- sphere2(cbind(x,y,z))
sphere2 <- function(x, y=NULL, z=NULL){
  if(is.null(y) && is.null(z)){
    if(is.vector(x)) return(new_sphere2(x[1],x[2],x[3]))
    return(new_sphere2(as.double(x[,1]), 
                       as.double(x[,2]),
                       as.double(x[,3])))
  }
  y <- as.double(y)
  z <- as.double(z)
  new_sphere2(x, y, z)
}

#' Constructor of class \code{sphere3} 
#' 
#' @param x A vector, or a matrix with 4 columns.
#' @param y A vector. May be missing when x is a matrix.
#' @param z A vector. May be missing when x is a matrix.
#' @param w A vector. May be missing when x is a matrix.
#' @return n x 4 matrix with coordinates \code{x}, \code{y}, \code{z}, and \code{w}.
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100)
#' y <- sqrt((1-x^2)/6)
#' z <- sqrt((1-x^2)/6*2)
#' w <- sqrt((1-x^2)/6*3)
#' s1 <- sphere3(x,y,z,w)
#' s2 <- sphere3(cbind(x,y,z,w))
sphere3 <- function(x, y=NULL, z=NULL, w=NULL){
  if(is.null(y) && is.null(z) && is.null(w)){
    if(is.vector(x)) return(new_sphere3(x[1],x[2],x[3],x[4]))
    return(new_sphere3(as.double(x[,1]), 
                       as.double(x[,2]),
                       as.double(x[,3]),
                       as.double(x[,4])))
  }
  y <- as.double(y)
  z <- as.double(z)
  w <- as.double(w)
  new_sphere3(x, y, z, w)
}


`[.sphere2` <- function(x, i, j, ...) {
  out <- unclass(x)
  if(is.vector(out)){
    out <- `[`(out,i)
    if(length(out) == 3){
      return(sphere2(out))
    }
    return(out)
  }
  if(missing(j)){
    if(missing(i)) i <- 1:nrow(out)
    out <- out[i,,...]
    return(sphere2(out))
  }
  if(length(j) != 3){
    out <- out[i,j,...]
    return(out)
  }
  return(sphere2(out[i,j,...]))
}

`[.sphere3` <- function(x, i, j, ...) {
  out <- unclass(x)
  if(is.vector(out)){
    out <- `[`(out,i)
    if(length(out) == 4){
      return(sphere3(out))
    }
    return(out)
  }
  if(missing(j)){
    if(missing(i)) i <- 1:nrow(out)
    out <- out[i,,...]
    return(sphere3(out))
  }
  if(length(j) != 4){
    out <- out[i,j,...]
    return(out)
  }
  return(sphere3(out[i,j,...]))
}

## plot methods for sphere 2 ---------------------------------------------------

#' Method of \code{plot} generic for \code{sphere2}
#' 
#' @param data object of class \code{sphere2}
#' @return No return, draws a plot.
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100)
#' y <- sqrt((1-x^2)/3)
#' z <- sqrt((1-x^2)/3*2)
#' s1 <- sphere2(x,y,z)
#' plot(s1, col="red")
plot.sphere2 <- function(data,...){
  rgl.sphgrid(col.long='white', col.lat='white', radaxis = FALSE)
  rgl.sphpoints(car2sph(data[,1],data[,2],data[,3]),...)
}

#' Draws a geodesic between two points.
#' 
#' There should be an existing plot to draw the geodesic onto.
#' Vectorized operations are not supported. Use \code{geodesics.sphere2}.
#' @param a A 3-vector of cartesian coord.
#' @param b A 3-vector of cartesian coord.
#' @return No return. Draws a geodesic. 
#' @export
#' @examples
#' plot(sphere2(rbind(c(1,0,0),c(0,0,1))))
#' geodesic.sphere2(c(1,0,0),c(0,0,1))
geodesic.sphere2 <- function(a,b, ...){
  d <- DIST.sphere2(a,b)
  if(abs(d-pi)<1e-5) message("unique geodesic does not exist!")
  
  logmap <- LOG.sphere2(a, b)
  theta <- seq(0, 1, length.out = 100)
  expmap <- sapply(1:100, function(i) EXP.sphere2(a, logmap*theta[i]))
  lines3d(t(expmap), ...)
}

#' Draw geodesics between data
#' 
#' Draw geodesics connecting each row of the given data onto the existing plot.
#' @param data An n x 3 matrix.
#' @return No return. Draws a geodesic. 
#' @export
#' @examples
#' plot(sphere2(rbind(c(1,0,0),c(1/sqrt(2),1/sqrt(2),0),c(1/sqrt(3),-1/sqrt(3),1/sqrt(3)))))
#' geodesics.sphere2(rbind(c(1,0,0),c(1/sqrt(2),1/sqrt(2),0),c(1/sqrt(3),-1/sqrt(3),1/sqrt(3))), col = "red")
geodesics.sphere2 <- function(data, ...){
  n <- nrow(data)
  for(i in 1:(n-1)){
    geodesic.sphere2(data[i,],data[i+1,], ...)
  }
}

#' Draws a full geodesic
#' 
#' @param geo A 2*3 matrix consisting of parameters of a geodesic. A point on S^2 and a direction of the geodesic at that point.
#' @return No return. Just draws a geodesic onto the existing plot.
#' @export
#' @examples
#' geo1 <- rbind(c(1/sqrt(2),1/sqrt(2),0),c(0,0,1))
#' plot(sphere2(c(1/sqrt(2),1/sqrt(2),0)))
#' geodesic_full.sphere2(geo1, col = "red", lwd = 3)
geodesic_full.sphere2 <- function(geo, ...){
  a <- geo[1,]
  dir <- geo[2,]/sqrt(sum(geo[2,]^2))
  theta <- seq(0, 2*pi, length.out = 100)
  expmap <- sapply(1:100, function(i) EXP.sphere2(a, dir*theta[i]))
  lines3d(t(expmap), ...)
}

## -----------------------------------------------------------------------------

#' Geodesic distance between two points. 
#' 
#' Vectorized operations are not supported.
#' @param a A 3-vector of length 1.
#' @param b A 3-vector of length 1.
#' @return The geodesic distance between \code{a} and \code{b}. 
#' @export
#' @examples
#' DIST.sphere2(c(1,0,0), c(0,1,0))
DIST.sphere2 <- function(a, b){ # a and b are given by 4-vectors (x,y,z,w)
  s <- sum(a*b)
  if(s>=1) return(0)
  if(s<=-1) return(pi)
  acos(sum(a*b))
}

#' Geodesic distance between two points. 
#' 
#' Vectorized operations are not supported.
#' @param a A 4-vector of length 1.
#' @param b A 4-vector of length 1.
#' @return The geodesic distance between \code{a} and \code{b}. 
#' @export
#' @examples
#' DIST.sphere3(c(1,0,0,0), c(0,1,0,0))
DIST.sphere3 <- function(a, b){ # a and b are given by 4-vectors (x,y,z,w)
  s <- sum(a*b)
  if(s>=1) return(0)
  if(s<=-1) return(pi)
  acos(sum(a*b))
}

#' Returns a logarithmic map
#' 
#' @param p A 3-vector on ambient space of S^2.
#' @param v A 3-vector on ambient space of S^2. 
#' @return Log_p(v) in R^3. 
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(5)
#' y <- -sqrt((1-x^2)/3)
#' z <- sqrt((1-x^2)/3*2)
#' s1 <- sphere2(x,y,z)
#' LOG(s1[1,],s1[2,])
LOG.sphere2 <- function(p,v){
  d <- DIST.sphere2(p,v)
  if(d < 1e-10) return(c(0,0,0))
  if(abs(d-pi)<1e-10) stop("log mapping not defined")
  dir <- v - p * cos(d)
  map <- d * dir/sqrt(sum(dir^2))
  return(unclass(map))
}

#' Returns a logarithmic map
#' 
#' @param p A 4-vector on ambient space of S^3.
#' @param v A 4-vector on ambient space of S^3. 
#' @return Log_p(v) in R^4. 
#' @export
LOG.sphere3 <- function(p,v){
  d <- DIST.sphere3(p,v)
  if(d < 1e-10) return(c(0,0,0,0))
  if(abs(d-pi)<1e-10) stop("log mapping not defined")
  dir <- v - p * cos(d)
  map <- d * dir/sqrt(sum(dir^2))
  return(unclass(map))
}

#' Returns an exponential map
#' 
#' @param p A 3-vector on S^2.
#' @param v A 3-vector on T_p(S^2). Thus, p must be orthogonal to v. 
#' @return Exp_p(v). 
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(5)
#' y <- -sqrt((1-x^2)/3)
#' z <- sqrt((1-x^2)/3*2)
#' s1 <- sphere2(x,y,z)
#' EXP(s1[1,],c(0,0,0))
#' EXP(s1[1,],c(-s1[1,2],s1[1,1],0))
EXP.sphere2 <- function(p,v){
  d <- sqrt(sum(v^2))
  if(d < 1e-10) return(p)
  res <- as.vector(p*cos(d) + v/d*sin(d))
  return(res)
}

#' Returns an exponential map
#' 
#' @param p A 4-vector on S^3.
#' @param v A 4-vector on T_p(S^3). Thus, p must be orthogonal to v.
#' @return Exp_p(v). 
#' @export
EXP.sphere3 <- function(p,v){
  d <- sqrt(sum(v^2))
  if(d < 1e-10) return(p)
  res <- as.vector(p*cos(d) + v/d*sin(d))
  return(res)
}

#' Parallel transport a vector on the tangent space via geodesic
#' 
#' @param start A 3-vector on S^2.
#' @param end A 3-vector on S^2. 
#' @param vec A 3-vector on T_start(S^2), orthogonal to start vector
#' @return Parallel transport \code{vec} from \code{start} to \code{end} via a geodesic. 
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(5)
#' y <- -sqrt((1-x^2)/3)
#' z <- sqrt((1-x^2)/3*2)
#' s1 <- sphere2(x,y,z)
#' PARALLEL(s1[1,],s1[2,],c(0,s1[1,3],-s1[1,2]))
PARALLEL.sphere2 <- function(start, end, vec){
  if(sum((start+end)^2) < 1e-7) stop("parallel transport to antipodals not supported")
  theta <- DIST.sphere2(start, end)
  if(sum(theta^2)< 1e-7) return(vec)
  e4 <- LOG.sphere2(start, end)
  e4 <- e4 / sqrt(sum(e4^2))
  a <- sum(e4*vec) # projection of vec to e4 direction
  stationary <- vec - a*e4
  res <- stationary + a*(e4*cos(theta) - start*sin(theta))
  return(res)
}



#' Parallel transport a vector on the tangent space via geodesic
#' 
#' @param start A 4-vector on S^3.
#' @param end A 4-vector on S^3. 
#' @param vec A 4-vector on T_start(S^3), orthogonal to start vector
#' @return Parallel transport \code{vec} from \code{start} to \code{end} via a geodesic. 
#' @export
PARALLEL.sphere3 <- function(start, end, vec){
  if(sum((start+end)^2) < 1e-7) stop("parallel transport to antipodals not supported")
  theta <- DIST.sphere3(start, end)
  if(sum(theta^2)< 1e-7) return(vec)
  e4 <- LOG.sphere3(start, end)
  e4 <- e4 / sqrt(sum(e4^2))
  a <- sum(e4*vec) # projection of vec to e4 direction
  stationary <- vec - a*e4
  res <- stationary + a*(e4*cos(theta) - start*sin(theta))
  return(res)
}

#' Project a point onto a geodesic
#' 
#' @param a A 3-vector on S^2.
#' @param geo A 2*3 matrix, first row is the point on the geodesic, and the second row is the direction of the geodesic at that point. Thus, two rows are orthogonal.
#' @return Projected point of \code{a} to \code{geo}.
#' @export
#' @examples
#' plot(sphere2(rbind(c(1,0,0),c(0,1/sqrt(2),1/sqrt(2)))))
#' geo <- rbind(c(0,1/sqrt(2),1/sqrt(2)),c(1/sqrt(2),-2,2))
#' geodesic_full.sphere2(geo,lwd=3)
#' p <- project.sphere2(c(1,0,0),geo)
#' geodesic.sphere2(p, c(1,0,0))
project.sphere2 <- function(a, geo){
  e1 <- geo[1,]
  dir <- geo[2,]
  e2 <- dir/sqrt(sum(dir^2))
  proj <- sum(e1*a)*e1 + sum(e2*a)*e2
  if(sum(proj^2) < 1e-7) stop("projection not defined")
  proj <- proj/sqrt(sum(proj^2))
  return(sphere2(proj))
}

#' Project a point onto a geodesic
#' 
#' @param a A 4-vector on S^3.
#' @param geo A 2*4 matrix, first row is the point on the geodesic, and the second row is the direction of the geodesic at that point. Thus, two rows are orthogonal.
#' @return Projected point of \code{a} to \code{geo}.
#' @export
project.sphere3 <- function(a, geo){
  e1 <- geo[1,]
  dir <- geo[2,]
  e2 <- dir/sqrt(sum(dir^2))
  proj <- sum(e1*a)*e1 + sum(e2*a)*e2
  if(sum(proj^2) < 1e-7) stop("projection not defined")
  proj <- proj/sqrt(sum(proj^2))
  return(sphere3(proj))
}

#' Unroll a piecewise geodesic curve
#' 
#' @param data knots of the curve.
#' @param ret.type dimension of the returned curve. If 2, return the parallel transported curve to the north pole, if 3, return the curve in the ambient space.
#' @return Unrolled curve of the data on T_{data[1,]}(S^2)
#' @export
#' @examples
#' dt <- sphere2(rbind(c(1,0,0),c(0,1,0),c(0,0,1)))
#' plot(dt)
#' geodesics.sphere2(dt, lwd = 3)
#' uroll <- UNROLL(dt, ret.type=3)
#' lines3d(uroll+matrix(rep(dt[1,],3),byrow=T,ncol=3))
UNROLL.sphere2 <- function(data, ret.type = 2){
  data <- as.matrix(data)
  n <- nrow(data)
  unroll.mat.3d <- matrix(ncol = 3, nrow = n)
  unroll.mat.3d[1,] <- c(0,0,0)
  for(i in 2:n){
    pp <- LOG.sphere2(data[i-1,],data[i,])
    if(i >= 3){
      for(j in (i-1):2){
        pp <- PARALLEL.sphere2(data[j,],data[j-1,],pp)
      }
    }
    unroll.mat.3d[i,] <- pp + unroll.mat.3d[i-1,]
  }
  if(ret.type == 3L) return(unroll.mat.3d) 
  pt <- t(sapply(1:n, function(i) PARALLEL.sphere2(data[1,],c(0,0,1),unroll.mat.3d[i,]) ))
  return(pt[,1:3])
}


#' Unroll a piecewise geodesic curve
#' 
#' @param data knots of the curve.
#' @param ret.type dimension of the returned curve. If 2, return the parallel transported curve to the north pole, if 3, return the curve in the ambient space.
#' @return Unrolled curve of the data on T_{data[1,]}(S^2)
#' @export
UNROLL.sphere3 <- function(data, ret.type = 3){
  data <- as.matrix(data)
  n <- nrow(data)
  unroll.mat.4d <- matrix(ncol = 4, nrow = n)
  unroll.mat.4d[1,] <- c(0,0,0,0)
  for(i in 2:n){ #i <- 2
    pp <- LOG.sphere3(data[i-1,],data[i,])
    if(i >= 3){
      for(j in (i-1):2){
        pp <- PARALLEL.sphere3(data[j,],data[j-1,],pp)
      }
    }
    unroll.mat.4d[i,] <- pp + unroll.mat.4d[i-1,]
  }
  if(ret.type == 4L) return(unroll.mat.4d) 
  pt <- t(sapply(1:n, function(i) PARALLEL.sphere3(data[1,],c(0,0,0,1),unroll.mat.4d[i,]) ))
  return(pt[,1:3])
}


#' Unwrap points with timestamp
#' 
#' @param data2 points on the 2-sphere
#' @param timestamp timestamp of each data, corresponding to each row of \code{data2}.
#' @param base base curve to unwrap w.r.t.
#' @return Unwrapped data on T_{base[1,]}(S^2). (Parallel transported to north pole)
#' @export
#' @examples
#' set.seed(0)
#' theta <- runif(10,0,360)
#' phi <- runif(10,0,180)
#' base <- sphere2(sph2car(theta, phi, radius=1))
#' data <- sphere2(sph2car(theta+rnorm(10,0,10), phi+rnorm(10,0,10), radius=1))
#' UNWRAP(data, 1:10, base)
UNWRAP.sphere2 <- function(data2, timestamp, base){
  data2 <- as.matrix(data2)
  n <- nrow(data2)
  base.unroll <- UNROLL.sphere2(base, ret.type = 3)
  data2.unwrap <- matrix(ncol=3, nrow=n)
  
  for(i in 1:n){
    if(timestamp[i] == 1){
      data2.unwrap[i,] <- LOG.sphere2(base[1,],data2[i,])
    }else{
      t <- timestamp[i]
      data2.log <- LOG.sphere2(base[t,],data2[i,])
      for(j in t:2){
        data2.log <- PARALLEL.sphere2(base[j,],base[j-1,],data2.log)
      }
      data2.unwrap[i,] <- data2.log + base.unroll[t,]
    }
  }
  pt <- t(sapply(1:n, function(i) PARALLEL.sphere2(base[1,],c(0,0,1),data2.unwrap[i,]) ))
  return(pt[,1:2])
}

#' Unwrap points with timestamp
#' 
#' @param data2 points on the 3-sphere
#' @param timestamp timestamp of each data, corresponding to each row of \code{data2}.
#' @param base base curve to unwrap w.r.t.
#' @return Unwrapped data on T_{base[1,]}(S^3). (Parallel transported to north pole)
#' @export
UNWRAP.sphere3 <- function(data2, timestamp, base){
  data2 <- as.matrix(data2)
  n <- nrow(data2)
  base.unroll <- UNROLL.sphere3(base, ret.type = 4)
  data2.unwrap <- matrix(ncol=4,nrow=n)
  
  for(i in 1:n){
    if(timestamp[i] == 1){
      data2.unwrap[i,] <- LOG.sphere3(base[1,],data2[i,])
    }else{
      t <- timestamp[i]
      data2.log <- LOG.sphere3(base[t,],data2[i,])
      for(j in t:2){
        data2.log <- PARALLEL.sphere3(base[j,],base[j-1,],data2.log)
      }
      data2.unwrap[i,] <- data2.log + base.unroll[t,]
    }
  }
  pt <- t(sapply(1:n, function(i) PARALLEL.sphere3(base[1,],c(0,0,0,1),data2.unwrap[i,]) ))
  return(pt[,1:3])
}

#' Wrap points with timestamp
#' 
#' @param data2 points on the 2-sphere (in 2-dim)
#' @param timestamp timestamp of each data, corresponding to each row of \code{data2}.
#' @param base base curve to wrap w.r.t.
#' @return Wrapped data on S^2
#' @export
#' @examples
#' set.seed(0)
#' theta <- runif(10,0,360)
#' phi <- runif(10,0,180)
#' base <- sphere2(sph2car(theta, phi, radius=1))
#' data <- sphere2(sph2car(theta+rnorm(10,0,10), phi+rnorm(10,0,10), radius=1))
#' uw <- UNWRAP(data, 1:10, base)
#' w <- WRAP.sphere2(uw, 1:10, base)
#' max(abs(data - w))
WRAP.sphere2 <- function(data2, timestamp, base){ 
  data2 <- as.matrix(data2)
  n <- nrow(data2)
  base.unroll <- UNROLL.sphere2(base, ret.type = 3)
  data2.wrap <- matrix(ncol=3,nrow=n)
  data2 <- t(sapply(1:n, function(i) PARALLEL.sphere2(c(0,0,1),base[1,],c(data2[i,],0)) ))
  ## now the data2 start at (0,0,0) and is in the space embedded in R^3
  ## the space is orthogonal to the base[1,]
  
  for(i in 1:n){
    if(timestamp[i] == 1){
      data2.wrap[i,] <- EXP.sphere2(base[1,],data2[i,])
    }else{
      t <- timestamp[i]
      v <- data2[i,] - base.unroll[t,]
      for(j in 2:t){
        v <- PARALLEL.sphere2(base[j-1,],base[j,],v)
      }
      data2.wrap[i,] <- EXP.sphere2(base[t,],v)
    }
  }
  return(sphere2(data2.wrap/sqrt(rowSums(data2.wrap^2))))
}

#' Wrap points with timestamp
#' 
#' @param data2 points on the 3-sphere (in 3-dim)
#' @param timestamp timestamp of each data, corresponding to each row of \code{data2}.
#' @param base base curve to wrap w.r.t.
#' @return Wrapped data on S^3
#' @export
WRAP.sphere3 <- function(data2, timestamp, base){ 
  data2 <- as.matrix(data2)
  n <- nrow(data2)
  base.unroll <- UNROLL.sphere3(base, ret.type = 4)
  data2.wrap <- matrix(ncol=4,nrow=n)
  data2 <- t(sapply(1:n, function(i) PARALLEL.sphere3(c(0,0,0,1),base[1,],c(data2[i,],0)) ))
  ## now the data2 start at (0,0,0,0) and is in the space embedded in R^4
  ## the space is orthogonal to the base[1,]
  
  for(i in 1:n){
    if(timestamp[i] == 1){
      data2.wrap[i,] <- EXP.sphere3(base[1,],data2[i,])
    }else{
      t <- timestamp[i]
      v <- data2[i,] - base.unroll[t,]
      for(j in 2:t){
        v <- PARALLEL.sphere3(base[j-1,],base[j,],v)
      }
      data2.wrap[i,] <- EXP.sphere3(base[t,],v)
    }
  }
  return(sphere3(data2.wrap/sqrt(rowSums(data2.wrap^2))))
}

#' Normalizing constant of 2-sphere
#' 
#' @param sigma scale param. of the norm. distn.
#' @return Normalizing constant of the Riemannian normal distribution of the given scale parameter.
#' @export
normal_const.sphere2 <- function(sigma){
  formula <- function(x) exp(-x^2/2/sigma^2)*sin(x)
  res <- integrate(formula, 0, pi)$value
  return(res * 2 * pi)
}

#' Normalizing constant of 3-sphere
#' 
#' @param sigma scale param. of the norm. distn.
#' @return Normalizing constant of the Riemannian normal distribution of the given scale parameter.
#' @export
normal_const.sphere3 <- function(sigma){
  formula <- function(x) exp(-x^2/2/sigma^2)*sin(x)^2
  res <- integrate(formula, 0, pi)$value
  return(res * 4 * pi)
}

#' Frechet mean (Intrinsic mean)
#' 
#' @param data each row is data point in 2-sphere
#' @param alpha step size
#' @return frechet mean
#' @export
#' @examples
#' set.seed(0)
#' theta <- runif(10,0,360)
#' phi <- runif(10,0,180)
#' base <- sphere2(sph2car(theta, phi, radius=1))
#' FRECHET(base)
FRECHET.sphere2 <- function (data, alpha = 1) 
{
  n <- nrow(data)
  mu <- data[1,,drop=FALSE]
  delta.magnitude <- 1
  
  while(delta.magnitude > 1e-3){
    delta <- c(0,0,0)
    for(j in 1:n){
      delta <- delta + LOG.sphere2(mu, data[j,])
    }
    delta.magnitude <- sqrt(sum(delta^2))
    mu <- EXP.sphere2(mu, alpha * delta/n)
  }
  return(mu)
}

#' Frechet mean (Intrinsic mean)
#' 
#' @param data each row is data point in 3-sphere
#' @param alpha step size
#' @return frechet mean
#' @export
FRECHET.sphere3 <- function (data, alpha = 1) 
{
  n <- nrow(data)
  mu <- data[1,,drop=FALSE]
  delta.magnitude <- 1
  
  while(delta.magnitude > 1e-3){
    delta <- c(0,0,0,0)
    for(j in 1:n){
      delta <- delta + LOG.sphere3(mu, data[j,])
    }
    delta.magnitude <- sqrt(sum(delta^2))
    mu <- EXP.sphere3(mu, alpha * delta/n)
  }
  return(mu)
}

#' Principal Geodesic Analysis (PGA)
#' 
#' @param data each row is data point in 2-sphere
#' @param fmean to reduce computation time, plug in the precomputed frechet mean
#' @return geodesic parametrized by 2*3 matrix
#' @export
PGA.sphere2 <- function(data, f.mean=NULL){
  n <- nrow(data)
  if(is.null(f.mean)) f.mean <- FRECHET.sphere2(data)
  theta <- seq(0,2*pi,length.out = 100)
  dir <- c(0,0,0)
  res.sum <- 0
  for(i in seq_along(theta)){
    xy <- c(cos(theta[i]),sin(theta[i]))
    xyz <- PARALLEL.sphere2(c(0,0,1),f.mean,c(xy,0))
    pv <- t(sapply(1:n, function(k) project.sphere2(data[k,],rbind(f.mean,xyz))))
    pv.sum <- sum(sapply(1:n, function(k) LOG.sphere2(f.mean,pv[k,]))^2)
    if(pv.sum > res.sum){
      res.sum <- pv.sum
      dir <- xyz
    }
    # cat(i, "\t")
  }
  return(rbind(f.mean,dir/sqrt(sum(dir^2))))
}

#' Principal Geodesic Analysis (PGA)
#' 
#' @param data each row is data point in 3-sphere
#' @param fmean to reduce computation time, plug in the precomputed frechet mean
#' @return geodesic parametrized by 2*4 matrix
#' @export
PGA.sphere3 <- function(data, f.mean=NULL){
  n <- nrow(data)
  if(is.null(f.mean)) f.mean <- FRECHET.sphere3(data)
  theta <- seq(0,2*pi,length.out = 30)
  phi <- seq(0,pi,length.out=30)
  dir <- c(0,0,0,0)
  res.sum <- 0
  for(i in seq_along(theta)){
    for(j in seq_along(phi)){
      xyz <- sph2car(theta[i], phi[j], radius = 1, deg = FALSE)
      xyzw <- PARALLEL.sphere3(c(0,0,0,1),f.mean,c(xyz,0))
      pv <- t(sapply(1:n, function(k) project.sphere3(data[k,],rbind(f.mean,xyzw))))
      pv.sum <- sum(sapply(1:n, function(k) LOG.sphere3(f.mean,pv[k,]))^2)
      if(pv.sum > res.sum){
        res.sum <- pv.sum
        dir <- xyzw
      }
    }
    cat(theta[i], "\t")
  }
  return(rbind(f.mean,dir/sqrt(sum(dir^2))))
}

#' Sigma: rule-of-thumb
#' 
#' 0.1n-nearest neighborhood, median*qchisq(0.8,2)
#' @param data each row is data point in 2-sphere
#' @return rule-of-thumb sigma value
#' @export
sigma_thumb.sphere2 <- function(data){
  n <- nrow(data)
  dist.mat <- matrix(ncol = n, nrow = n)
  for(i in 1:n){
    for(j in 1:n){
      dist.mat[i,j] <- DIST.sphere2(data[i,], data[j,])
    }
  }
  
  # mean 0.1-nearest neigh distance
  nn.cov <- sapply(1:n, function(r) sqrt(sum(sort(dist.mat[r,])[2:(n/10+1)]^2)/(n/10)))
  sigma.thumb <- median(nn.cov) * sqrt(qchisq(0.8,2))
  sigma.thumb
}

#' Sigma: rule-of-thumb
#' 
#' 0.1n-nearest neighborhood, median*qchisq(0.8,3)
#' @param data each row is data point in 2-sphere
#' @return rule-of-thumb sigma value
#' @export
sigma_thumb.sphere3 <- function(data){
  n <- nrow(data)
  dist.mat <- matrix(ncol = n, nrow = n)
  for(i in 1:n){
    for(j in 1:n){
      dist.mat[i,j] <- DIST.sphere3(data[i,], data[j,])
    }
  }
  
  # mean 0.1-nearest neigh distance
  nn.cov <- sapply(1:n, function(r) sqrt(sum(sort(dist.mat[r,])[2:(n/10+1)]^2)/(n/10)))
  sigma.thumb <- median(nn.cov) * sqrt(qchisq(0.8,3))
  sigma.thumb
}


PPC.sphere2 <- function(data, lambda, max.iter = 20, vis = FALSE, geo.pga = NULL){
  n <- nrow(data)
  if(is.null(geo.pga)) geo.pga <- PGA.sphere2(data)
  f <- t(apply(data, 1, function(row) project.sphere2(row, geo.pga)))
  d <- matrix(0,ncol=n,nrow=n)
  for(i in 1:n){
    for(j in 1:n){
      d[i,j] <- DIST.sphere2(f[i,],f[j,])
    }
  }
  fit <- isoMDS(d, k=1)$points
  
  pbit <- progress_bar$new(
    format = paste("ITER :its / Gradient Ascent at :what in :elapsed"),
    clear = FALSE, total = 1e7, width = 60)
  pbit$tick(tokens = list(its = 1,what="waiting..."))
  # message("ITER 1")
  f <- f[order(fit),]
  sigma.max <- Inf
  a.const <- 1
  a <- vector(length = n-1)
  for(i in seq_along(a)){
    a[i] <- DIST.sphere2(f[i,],f[i+1,])
  }
  a <- c(0,cumsum(a))/sum(a) *a.const
  sigma.thumb <- sigma_thumb.sphere2(data)
  sigma <- rep(sigma.thumb,n)
  v <- rep(1/n,n)
  w <- matrix(ncol = n, nrow = n)
  for(i in 1:n){
    for(k in 1:n){
      d <- DIST.sphere2(data[i,],f[k,])
      g <- exp(-d^2/2/sigma[k]^2)/normal_const.sphere2(sigma[k])
      w[i,k] <- g * v[k]
    }
  }
  w <- w/rowSums(w)
  timestamp <- apply(w, 1, which.max)
  
  ## iteration 1
  data.unwrap <- UNWRAP.sphere2(data, timestamp, f)
  iter1 <- ppc.iter.e(data.unwrap, w, a, v, sigma, sigma.max = sigma.max, lambda)
  f.new <- WRAP.sphere2(iter1, 1:n, f)
  ## SGD
  tau <- 1/sigma^2
  sgd.const2 <-  sapply(1:n, function(k) 
    sum(sapply(1:n, function(i) -1/2*w[i,k]*DIST.sphere2(data[i,],f.new[k,])^2)))
  sgd.const1 <- vector(length = n)
  alpha <- rep(1, n) #step-size
  grad <- rep(1,n)
  iters <- 0
  # pb <- progress_bar$new(
  #   format = paste("  Gradient Ascent at :what, in :elapsed"),
  #   clear = FALSE, total = 1e7, width = 60)
  while(sum(grad[alpha!=0]^2)>1e-7){
    pbit$tick(tokens = list(its=1,what = sum(grad[alpha!=0]^2)))
    iters <- iters + 1
    # if(vis) cat(sum(grad[alpha!=0]^2), "\t")
    for(k in 1:n){
      formula <- function(r) r^2/2*exp(-r^2/2/sigma[k]^2)*sin(r)
      sgd.const1[k] <- 2*pi*integrate(formula, 0, pi)$value/normal_const.sphere2(sigma[k])
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
  
  if(vis){
    plot(data)
    geodesics.sphere2(f,col="red",lwd=3)
    geodesics.sphere2(f.new,col="blue",lwd=3)
    text3d(0,0,1.5,"Iter 1 / blue: new, red: old")
  }
  
  ## iteration >= 2
  kk <- 2
  while(TRUE){
    if(max(sapply(1:n, function(i) DIST.sphere2(f[i,],f.new[i,]))) < 1e-2 || kk > max.iter) break()
    pbit$tick(tokens = list(its = kk,what="waiting..."))
    # message(paste("\nITER ",kk))
    f <- f.new
    a <- vector(length = n-1)
    for(i in seq_along(a)){
      a[i] <- DIST.sphere2(f[i,],f[i+1,])
    }
    a <- c(0,cumsum(a))/sum(a) *a.const
    v <- colSums(w)/n
    w <- matrix(ncol = n, nrow = n)
    for(i in 1:n){
      for(k in 1:n){
        d <- DIST.sphere2(data[i,],f[k,])
        g <- exp(-d^2/2/sigma[k]^2)/normal_const.sphere2(sigma[k])
        w[i,k] <- g * v[k]
      }
    }
    w <- w/rowSums(w)
    if(!all(!is.nan(w)) || !all(colSums(w) != 0)){
      message("stopped because some latent variables are useless")
      break()
    }
    timestamp <- apply(w, 1, which.max)
    data.unwrap <- UNWRAP.sphere2(data, timestamp, f)
    iter2 <- ppc.iter.e(data.unwrap, w, a, v, sigma, sigma.max = sigma.max, lambda)
    
    if(vis){
      plot(data.unwrap,pch=".")
      text(data.unwrap[,1],data.unwrap[,2],labels=as.character(timestamp))
      points(iter2,pch="+")
      for(i in 1:n) normal.e.plot(iter2[i,],sigma[i], col="gray")
    }
    
    f.new <- WRAP.sphere2(iter2, 1:n, f)
    tau <- 1/sigma^2
    tau.old <- tau
    sigma.old <- sigma
    sgd.const2 <-  sapply(1:n, function(k) 
      sum(sapply(1:n, function(i) -1/2*w[i,k]*DIST.sphere2(data[i,],f.new[k,])^2)))
    sgd.const1 <- vector(length = n)
    alpha <- rep(1, n) #step-size
    grad <- rep(1,n)
    iters <- 0
    # pb <- progress_bar$new(
    #   format = paste("  Gradient Ascent at :what, in :elapsed"),
    #   clear = FALSE, total = 1e7, width = 60)
    while(sum(grad[alpha!=0]^2)>1e-7){
      pbit$tick(tokens = list(its=kk,what = sum(grad[alpha!=0]^2)))
      iters <- iters + 1
      # if(vis) cat(sum(grad[alpha!=0]^2), "\t")
      for(k in 1:n){
        formula <- function(r) r^2/2*exp(-r^2/2/sigma[k]^2)*sin(r)
        sgd.const1[k] <- 2*pi*integrate(formula, 0, pi)$value/normal_const.sphere2(sigma[k])
      }
      sgd.const1 <- colSums(w) * sgd.const1
      grad.old <- grad
      grad <- sgd.const1 + sgd.const2
      if(!all(!is.na(grad))) return(list(f,sigma.old))
      if(sum((grad.old-grad)^2)<1e-7) alpha <- alpha * 1.1
      alpha[which(abs(grad) < 1e-7)] <- 0
      if(iters > 1000) return(list(f,sigma.old))
      tau <- tau + alpha*grad
      tau[tau < 0] <- min(abs(tau[tau!=0]))
      sigma <- sqrt(1/tau)
    }
    
    if(is.na(sum(f.new))){
      message("f is not numeric. NA")
      return(list(f, sigma))
    }
    
    if(vis){
      plot(data)
      lines3d(f,col="red",lwd=3)
      lines3d(f.new,col="blue",lwd=3)
      text3d(0,0,1.5,paste("Iter ",kk," / blue: new, red: old / ",max(sapply(1:n, function(i) DIST.sphere2(f[i,],f.new[i,]))) ))
    }
    kk <- kk + 1
  }
  return(list(f.new, sigma))
}

PPC.sphere3 <- function(data, lambda, max.iter = 20, geo.pga = NULL){
  n <- nrow(data)
  if(is.null(geo.pga)) geo.pga <- PGA.sphere3(data)
  f <- t(apply(data, 1, function(row) project.sphere3(row, geo.pga)))
  d <- matrix(0,ncol=n,nrow=n)
  for(i in 1:n){
    for(j in 1:n){
      d[i,j] <- DIST.sphere3(f[i,],f[j,])
    }
  }
  fit <- isoMDS(d, k=1)$points
  
  message("ITER 1")
  f <- f[order(fit),]
  sigma.max <- Inf
  a.const <- 1
  a <- vector(length = n-1)
  for(i in seq_along(a)){
    a[i] <- DIST.sphere3(f[i,],f[i+1,])
  }
  a <- c(0,cumsum(a))/sum(a) *a.const
  sigma.thumb <- sigma_thumb.sphere3(data)
  sigma <- rep(sigma.thumb,n)
  v <- rep(1/n,n)
  w <- matrix(ncol = n, nrow = n)
  for(i in 1:n){
    for(k in 1:n){
      d <- DIST.sphere3(data[i,],f[k,])
      g <- exp(-d^2/2/sigma[k]^2)/normal_const.sphere3(sigma[k])
      w[i,k] <- g * v[k]
    }
  }
  w <- w/rowSums(w)
  timestamp <- apply(w, 1, which.max)
  
  ## iteration 1
  data.unwrap <- UNWRAP.sphere3(data, timestamp, f)
  iter1 <- ppc.iter.e.3d(data.unwrap, w, a, v, sigma, sigma.max = sigma.max, lambda)
  f.new <- WRAP.sphere3(iter1, 1:n, f)
  ## SGD
  tau <- 1/sigma^2
  sgd.const2 <-  sapply(1:n, function(k) 
    sum(sapply(1:n, function(i) -1/2*w[i,k]*DIST.sphere3(data[i,],f.new[k,])^2)))
  sgd.const1 <- vector(length = n)
  alpha <- rep(1, n) #step-size
  grad <- rep(1,n)
  iters <- 0
  pb <- progress_bar$new(
    format = paste("  Gradient Ascent at :what, in :elapsed"),
    clear = FALSE, total = 1e7, width = 60)
  while(sum(grad[alpha!=0]^2)>1e-7){
    pb$tick(tokens = list(what = sum(grad[alpha!=0]^2)))
    iters <- iters + 1
    # if(vis) cat(sum(grad[alpha!=0]^2), "\t")
    for(k in 1:n){
      formula <- function(r) r^2/2*exp(-r^2/2/sigma[k]^2)*sin(r)^2
      sgd.const1[k] <- 4*pi*integrate(formula, 0, pi)$value/normal_const.sphere3(sigma[k])
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
    if(max(sapply(1:n, function(i) DIST.sphere3(f[i,],f.new[i,]))) < 1e-2 || kk > max.iter) break()
    message(paste("\nITER ",kk))
    f <- f.new
    a <- vector(length = n-1)
    for(i in seq_along(a)){
      a[i] <- DIST.sphere3(f[i,],f[i+1,])
    }
    a <- c(0,cumsum(a))/sum(a) *a.const
    v <- colSums(w)/n
    w <- matrix(ncol = n, nrow = n)
    for(i in 1:n){
      for(k in 1:n){
        d <- DIST.sphere3(data[i,],f[k,])
        g <- exp(-d^2/2/sigma[k]^2)/normal_const.sphere3(sigma[k])
        w[i,k] <- g * v[k]
      }
    }
    w <- w/rowSums(w)
    if(!all(!is.nan(w)) || !all(colSums(w) != 0)){
      message("stopped because some latent variables are useless")
      break()
    }
    timestamp <- apply(w, 1, which.max)
    data.unwrap <- UNWRAP.sphere3(data, timestamp, f)
    iter2 <- ppc.iter.e.3d(data.unwrap, w, a, v, sigma, sigma.max = sigma.max, lambda)
    
    f.new <- WRAP.sphere3(iter2, 1:n, f)
    tau <- 1/sigma^2
    tau.old <- tau
    sigma.old <- sigma
    sgd.const2 <-  sapply(1:n, function(k) 
      sum(sapply(1:n, function(i) -1/2*w[i,k]*DIST.sphere3(data[i,],f.new[k,])^2)))
    sgd.const1 <- vector(length = n)
    alpha <- rep(1, n) #step-size
    grad <- rep(1,n)
    iters <- 0
    pb <- progress_bar$new(
      format = paste("  Gradient Ascent at :what, in :elapsed"),
      clear = FALSE, total = 1e7, width = 60)
    while(sum(grad[alpha!=0]^2)>1e-7){
      pb$tick(tokens = list(what = sum(grad[alpha!=0]^2)))
      iters <- iters + 1
      # if(vis) cat(sum(grad[alpha!=0]^2), "\t")
      for(k in 1:n){
        formula <- function(r) r^2/2*exp(-r^2/2/sigma[k]^2)*sin(r)^2
        sgd.const1[k] <- 4*pi*integrate(formula, 0, pi)$value/normal_const.sphere3(sigma[k])
      }
      sgd.const1 <- colSums(w) * sgd.const1
      grad.old <- grad
      grad <- sgd.const1 + sgd.const2
      if(sum((grad.old-grad)^2)<1e-7) alpha <- alpha * 1.1
      alpha[which(abs(grad) < 1e-7)] <- 0
      if(iters > 1000) alpha <- alpha/2
      tau <- tau + alpha*grad
      tau[tau < 0] <- min(abs(tau[tau!=0]))
      sigma <- sqrt(1/tau)
    }
    
    if(is.na(sum(f.new))){
      message("f is not numeric. NA")
      return(list(f, sigma))
    }
    kk <- kk + 1
  }
  return(list(f.new, sigma,w))
}

