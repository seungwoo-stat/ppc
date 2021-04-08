source("half_plane_funs.R")

## Simulation 1 ----------------------------------------------------------------
## data generation
n <- 100
set.seed(0)
t <- runif(n,-3,3)
data1 <- matrix(cbind(t,2+3/4*abs(t)), nrow=n, ncol=2)
err <- rnorm(2*n,0,1/8)
for(i in 1:n){
  data1[i,] <- EXP.half_plane(data1[i,],c(err[i],err[n+i]))
}
data1 <- half_plane(data1)

## true curve
xt1 <- seq(-3,3,by = 0.01)
yt1 <- 2 + 3/4*abs(xt1)

## PGA
fm1 <- FRECHET(data1)
geo.pga1 <- PGA(data1)
plot(data1, xlim = c(-4,4), ylim = c(0,6))
lines(xt1, yt1)
points(fm1[1], fm1[2], pch = "f", col = "red")
geodesic_full.half_plane(geo.pga1, col = "red")
title("PGA and the Frechet mean")

## PPC
PPC1 <- PPC(data1, lambda = 0.03, max.iter = 20, vis = TRUE)


## Simulation 2 ----------------------------------------------------------------
n <- 100
set.seed(1)
t <- runif(n,-3,3)
data2 <- matrix(cbind(t,10), nrow=n, ncol=2)
err <- rnorm(2*n,0,1/32)
for(i in 1:n){
  data2[i,] <- EXP.half_plane(data2[i,],c(err[i],err[n+i]))
}
data2 <- half_plane(data2)

## true curve
xt2 <- seq(-3,3,by = 0.01)
yt2 <- rep(10,length(xt2))

## PGA
fm2 <- FRECHET(data2)
geo.pga2 <- PGA(data2)
par(mfrow=c(1,1))
plot(data2, xlim = c(-4,4), ylim = c(8,12))
lines(xt2, yt2)
points(fm2[1], fm2[2], pch = "f", col = "red")
geodesic_full.half_plane(geo.pga2, col = "red")
title("PGA and the Frechet mean")

## PPC
PPC2 <- PPC(data2, lambda = 2, max.iter = 20, vis = TRUE)


## Simulation 3 ----------------------------------------------------------------
n <- 100
set.seed(1)
t <- runif(n,-3,3)
data3 <- matrix(cbind(t,5-sin(t)), nrow=n, ncol=2)
err <- rnorm(2*n,0,1/16)
for(i in 1:n){
  data3[i,] <- EXP.half_plane(data3[i,],c(err[i],err[n+i]))
}
data3 <- half_plane(data3)

## true curve
xt3 <- seq(-3,3,by = 0.01)
yt3 <- 5-sin(xt3)

## PGA
fm3 <- FRECHET(data3)
geo.pga3 <- PGA(data3)
par(mfrow=c(1,1))
plot(data3, xlim = c(-4,4), ylim = c(2,7))
lines(xt3, yt3)
points(fm3[1], fm3[2], pch = "f", col = "red")
geodesic_full.half_plane(geo.pga3, col = "red")
title("PGA and the Frechet mean")

## PPC
PPC3 <- PPC(data3, lambda = 0.05, max.iter = 25, vis = TRUE)


## Simulation 4 ----------------------------------------------------------------
n <- 100
set.seed(0)
t <- runif(n,-3,3)
data4 <- matrix(cbind(t,10+1/2*t^2), nrow=n, ncol=2)
err <- rnorm(2*n,0,1/16)
for(i in 1:n){
  data4[i,] <- EXP.half_plane(data4[i,],c(err[i],err[n+i]))
}
data4 <- half_plane(data4)

## true curve
xt4 <- seq(-3,3,by = 0.01)
yt4 <- 10+1/2*xt4^2

## PGA
fm4 <- FRECHET(data4)
geo.pga4 <- PGA(data4)
par(mfrow=c(1,1))
plot(data4, xlim = c(-4,4), ylim = c(5,20))
lines(xt4, yt4)
points(fm4[1], fm4[2], pch = "f", col = "red")
geodesic_full.half_plane(geo.pga4, col = "red")
title("PGA and the Frechet mean")

## PPC
PPC4 <- PPC(data4, lambda = 0.1, max.iter = 25, vis = T)


## PLOTS -----------------------------------------------------------------------
par(mfrow=c(1,4))
plot(data1, pch = 20)
lines(xt1, yt1, lty = 2, lwd=2)
geodesic_full.half_plane(geo.pga1, col = "gray", lwd=2)
lines(PPC1[[1]], col = "red", lwd=2)
# title("Simulation Set 1")
# legend("topright", legend=c("Data Point","True Curve","Initial Curve (PGA)","PPC"), 
       # col=c("black","black","gray","red"), lty=c(NA,2,1,1), pch=c(20,NA,NA,NA), lwd=2)

plot(data2, pch = 20)
lines(xt2, yt2, lty=2, lwd=2)
geodesic_full.half_plane(geo.pga2, col = "gray", lwd=2)
lines(PPC2[[1]], col = "red", lwd=2)
# title("Simulation Set 2")
# legend("topright", legend=c("Data Point","True Curve","Initial Curve (PGA)","PPC"), 
#        col=c("black","black","gray","red"), lty=c(NA,2,1,1), pch=c(20,NA,NA,NA), lwd=2)

plot(data3, pch = 20)
lines(xt3, yt3, lty=2, lwd=2)
geodesic_full.half_plane(geo.pga3, col = "gray", lwd=2)
lines(PPC3[[1]], col = "red", lwd=2)
# title("Simulation Set 3")
# legend("topright", legend=c("Data Point","True Curve","Initial Curve (PGA)","PPC"), 
#        col=c("black","black","gray","red"), lty=c(NA,2,1,1), pch=c(20,NA,NA,NA), lwd=2)

plot(data4, pch = 20)
lines(xt4, yt4, lty=2, lwd=2)
geodesic_full.half_plane(geo.pga4, col = "gray", lwd=2)
lines(PPC4[[1]], col = "red", lwd=2)
# title("Simulation Set 4")
# legend("topright", legend=c("Data Point","True Curve","Initial Curve (PGA)","PPC"), 
#        col=c("black","black","gray","red"), lty=c(NA,2,1,1), pch=c(20,NA,NA,NA), lwd=2)


## Comparision -----------------------------------------------------------------
## (1) PPC on data4
## (2) HPC on data4
## (3) PGA on data4
library(RColorBrewer)
display.brewer.all()

HPC4_1 <- PC_Hauberg.half_plane(data4, q=0.12, max.iter = 20)
HPC4_2 <- PC_Hauberg.half_plane(data4, q=0.13, max.iter = 20)
HPC4_3 <- PC_Hauberg.half_plane(data4, q=0.14, max.iter = 20)
HPC4_4 <- PC_Hauberg.half_plane(data4, q=0.15, max.iter = 20)
HPC4_5 <- PC_Hauberg.half_plane(data4, q=0.16, max.iter = 20)
PPC4_1 <- PPC(data4, lambda = 0.005, max.iter = 25, vis = T)
PPC4_2 <- PPC(data4, lambda = 0.1, max.iter = 25, vis = F)
PPC4_3 <- PPC(data4, lambda = 2, max.iter = 25, vis = T)
geo.pga4 <- PGA(data4)

## PLOT ------------------------------------------------------------------------
par(mfrow=c(2,3))
plot(data1, pch = 20); lines(xt1, yt1, lty=2, lwd=2)
geodesic_full.half_plane(geo.pga1, col = "gray", lwd=2)
lines(PPC1[[1]], col = "red", lwd=2); title("(a) V-shape")
plot(data2, pch = 20); lines(xt2, yt2, lty=2, lwd=2)
geodesic_full.half_plane(geo.pga2, col = "gray", lwd=2)
lines(PPC2[[1]], col = "red", lwd=2); title("(b) Straight")
plot(data3, pch = 20); lines(xt3, yt3, lty=2, lwd=2)
geodesic_full.half_plane(geo.pga3, col = "gray", lwd=2)
lines(PPC3[[1]], col = "red", lwd=2); title("(c) Wave")

plot(data4, pch = 20); lines(xt4, yt4, lty=2, lwd=2)
lines(PPC4_1[[1]], lwd=3, col = brewer.pal(8, "RdYlGn")[1])
lines(PPC4_2[[1]], lwd=3, col = brewer.pal(8, "RdYlGn")[2])
lines(PPC4_3[[1]], lwd=3, col = brewer.pal(8, "RdYlGn")[3])
legend("topright", c(expression(paste(lambda,"=",0.005)),
                     expression(paste(lambda,"=",0.1)),
                     expression(paste(lambda,"=",2))), 
       col=brewer.pal(8, "RdYlGn")[1:3], lty=1,lwd=3)
title("(d) Parabola, PPC")

plot(data4, pch = 20); lines(xt4, yt4, lty=2, lwd=2)
geodesics.half_plane(HPC4_1, col = brewer.pal(11, "RdYlGn")[1], lwd = 3)
geodesics.half_plane(HPC4_2, col = brewer.pal(11, "RdYlGn")[2], lwd = 3)
geodesics.half_plane(HPC4_3, col = brewer.pal(11, "RdYlGn")[3], lwd = 3)
geodesics.half_plane(HPC4_4, col = brewer.pal(11, "RdYlGn")[4], lwd = 3)
geodesics.half_plane(HPC4_5, col = brewer.pal(11, "RdYlGn")[5], lwd = 3)
legend("topright", c(paste0("q=0.1",2:6)), col=brewer.pal(11, "RdYlGn")[1:5], lty=1,lwd=3)
title("(e) Parabola, HPC")

plot(data4, pch = 20); lines(xt4, yt4, lty=2, lwd=2)
geodesic_full.half_plane(geo.pga4, col = "red", lwd=3)
title("(f) Parabola, PGA")




