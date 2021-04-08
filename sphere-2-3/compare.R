source("sphere_funs.R")

## wave data
n <- 100
to <- 2*pi*3/4
t <- seq(0,to,length.out = 1000)
true <- sph2car(t,sin(3*t)/3,radius = 1,deg=FALSE)
plot(sphere2(true))
dist_true <- sapply(1:999, function(i) DIST.sphere2(true[i,],true[i+1,]))
dist_true_cum <- c(0,cumsum(dist_true))
true_dist_seq <- seq(0,sum(dist_true),length = 100)
true_points <- true[sapply(1:100, function(i) which.min(dist_true_cum < true_dist_seq[i])),]

## change this error level, and test the performance of several methods related to PPC.
sigma <- 0.1

lambda0 <- 0.003
q0 <- 0.04

## PPC -------------------------------------------------------------------------
## Repeat 50 times

CRITERIA <- vector(length = 50)
for(rep in 1:50){
  message(rep)
  set.seed(rep)
  
  theta2 <- runif(n,0,to)
  phi2 <- sin(3*theta2)/3
  data2 <- sphere2(sph2car(theta2,phi2,radius=1,deg=FALSE))
  for(i in 1:n){
    pp <- PARALLEL(sphere2(0,0,1),data2[i,],c(rnorm(2,0,sigma),0))
    data2[i,] <- EXP(data2[i,],pp)
  }
  
  ## PPC
  PPC2 <- PPC(data2, lambda=lambda0, max.iter = 20, vis = F, geo.pga = rbind(c(1,0,0),c(0,1,0)))
  PPC2 <- PPC2[[1]]
  ppc2 <- lapply(1:(n-1), function(i) t(sapply(0:9, function(j) (9-j)/9*PPC2[i,]+j/9*PPC2[i+1,])))
  ppc2 <- do.call('rbind',ppc2)
  ppc2 <- ppc2/sqrt(rowSums(ppc2^2)) # dense points in fitted curve
  RE <- vector(length = 100)
  for(i in 1:100){
    RE[i] <- min(sapply(1:990, function(j) DIST.sphere2(true_points[i,],ppc2[j,])))
  }
  CRITERIA[rep] <- sum(RE)/100
  message(CRITERIA[rep])
}
CRITERIA
mean(CRITERIA); sd(CRITERIA)


## HPC -------------------------------------------------------------------------

CRITERIA_HPC <- vector(length = 50)
for(rep in 1:50){
  message(rep)
  set.seed(rep)
  
  theta2 <- runif(n,0,to)
  phi2 <- sin(3*theta2)/3
  data2 <- sphere2(sph2car(theta2,phi2,radius=1,deg=FALSE))
  for(i in 1:n){
    pp <- PARALLEL(sphere2(0,0,1),data2[i,],c(rnorm(2,0,sigma),0))
    data2[i,] <- EXP(data2[i,],pp)
  }
  
  HPC2 <- SPC.Hauberg(car2sph(data2[,1],data2[,2],data2[,3])[,1:2], q=q0, deletePoints = TRUE)$prin.curves
  geodesics.sphere2(true_points, lwd = 5)
  HPC2 <- unique(HPC2)
  HPC2 <- sph2car(HPC2[,1],HPC2[,2])
  HPC2 <- lapply(1:(dim(HPC2)[1]-1), function(i) t(sapply(0:9, function(j) (9-j)/9*HPC2[i,]+j/9*HPC2[i+1,])))
  HPC2 <- do.call('rbind',HPC2)
  HPC2 <- HPC2/sqrt(rowSums(HPC2^2)) # dense points in fitted curve
  RE <- vector(length = 100)
  for(i in 1:100){
    RE[i] <- min(sapply(1:(dim(HPC2)[1]), function(j) DIST.sphere2(true_points[i,],HPC2[j,])))
  }
  CRITERIA_HPC[rep] <- sum(RE)/100
  message(CRITERIA_HPC[rep])
}
CRITERIA_HPC
mean(CRITERIA_HPC); sd(CRITERIA_HPC)


## PGA -------------------------------------------------------------------------

CRITERIA_PGA <- vector(length = 50)
for(rep in 1:50){
  message(rep)
  set.seed(rep)
  
  theta2 <- runif(n,0,to)
  phi2 <- sin(3*theta2)/3
  data2 <- sphere2(sph2car(theta2,phi2,radius=1,deg=FALSE))
  for(i in 1:n){
    pp <- PARALLEL(sphere2(0,0,1),data2[i,],c(rnorm(2,0,sigma),0))
    data2[i,] <- EXP(data2[i,],pp)
  }
  
  PGA2 <- PGA.sphere2(data2)
  pga_t <- seq(0,2*pi, length.out = 1000)
  PGA2 <- t(sapply(1:1000,function(i) EXP.sphere2(PGA2[1,],pga_t[i]*PGA2[2,])))
  
  RE <- vector(length = 100)
  for(i in 1:100){
    RE[i] <- min(sapply(1:(dim(PGA2)[1]), function(j) DIST.sphere2(true_points[i,],PGA2[j,])))
  }
  CRITERIA_PGA[rep] <- sum(RE)/100
  message(CRITERIA_PGA[rep])
}
CRITERIA_PGA
mean(CRITERIA_PGA); sd(CRITERIA_PGA)


## PNS -------------------------------------------------------------------------

CRITERIA_PNS <- vector(length = 50)
for(rep in 1:50){
  message(rep)
  set.seed(rep)
  
  theta2 <- runif(n,0,to)
  phi2 <- sin(3*theta2)/3
  data2 <- sphere2(sph2car(theta2,phi2,radius=1,deg=FALSE))
  for(i in 1:n){
    pp <- PARALLEL(sphere2(0,0,1),data2[i,],c(rnorm(2,0,sigma),0))
    data2[i,] <- EXP(data2[i,],pp)
  }
  
  PNS2 <- PrincipalCircle(car2sph(data2[,1],data2[,2],data2[,3])[,1:2])
  PNS2 <- GenerateCircle(PNS2[1:2], radius = PNS2[3])
  PNS2 <- sph2car(PNS2[,1],PNS2[,2])
  
  RE <- vector(length = 100)
  for(i in 1:100){
    RE[i] <- min(sapply(1:(dim(PNS2)[1]), function(j) DIST.sphere2(true_points[i,],PNS2[j,])))
  }
  CRITERIA_PNS[rep] <- sum(RE)/100
  message(CRITERIA_PNS[rep])
}
CRITERIA_PNS
mean(CRITERIA_PNS); sd(CRITERIA_PNS)


## SPC -------------------------------------------------------------------------

CRITERIA_SPC <- vector(length = 50)
for(rep in 1:50){
  message(rep)
  set.seed(rep)
  
  theta2 <- runif(n,0,to)
  phi2 <- sin(3*theta2)/3
  data2 <- sphere2(sph2car(theta2,phi2,radius=1,deg=FALSE))
  for(i in 1:n){
    pp <- PARALLEL(sphere2(0,0,1),data2[i,],c(rnorm(2,0,sigma),0))
    data2[i,] <- EXP(data2[i,],pp)
  }
  
  SPC2 <- SPC(car2sph(data2[,1],data2[,2],data2[,3])[,1:2], q=q0, deletePoints = TRUE)$prin.curves
  geodesics.sphere2(true_points, lwd=5)
  SPC2 <- unique(SPC2)
  SPC2 <- sph2car(SPC2[,1],SPC2[,2])
  SPC2 <- SPC2[c(rep(TRUE,69),!(SPC2[70:dim(SPC2)[1],1] > 0)),]
  SPC2 <- lapply(1:(dim(SPC2)[1]-1), function(i) t(sapply(0:9, function(j) (9-j)/9*SPC2[i,]+j/9*SPC2[i+1,])))
  SPC2 <- do.call('rbind',SPC2)
  SPC2 <- SPC2/sqrt(rowSums(SPC2^2)) # dense points in fitted curve
  RE <- vector(length = 100)
  for(i in 1:100){
    RE[i] <- min(sapply(1:(dim(SPC2)[1]), function(j) DIST.sphere2(true_points[i,],SPC2[j,])))
  }
  CRITERIA_SPC[rep] <- sum(RE)/100
  message(CRITERIA_SPC[rep])
}
CRITERIA_SPC
mean(CRITERIA_SPC); sd(CRITERIA_SPC)


